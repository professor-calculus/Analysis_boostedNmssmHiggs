// CPP headers
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

// ROOT headers
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

// CMSSW headers
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

// Headers from this package
#include "Analysis/Analysis_boostedNmssmHiggs/interface/Kinematics.h"
// #include "Analysis/Analysis_boostedNmssmHiggs/interface/PlottingInvestigatePreSelection.h"


/* TODO

no track info for lepton veto atm


*/

// PAT TUPLE FORMAT ONLY !!!

// preliminary running, compile with scram b and then
// $ investigatePreSelection inputfiles=XYZ outputfile=ABC orderedsecondaryfiles=0
// ($CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/InvestigatePreSelection/InvestigatePreSelection)
// nb. if you are running it on DICE, include the word runOnDice at the end of the arguments of the executable

/* Notes on runOnDice mode
watchout...with this toggle the executable can now overwrite outputs!!!
hacks the original script abit so it can do some things slightly differently
it does so using the quantity: bool runOnDice
*/

// NB: *known bug* if running on just plotting mode you still need to give it an inputfiles argument so that the script can get a title without seg faulting!!!

void CreateHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&);
void WriteHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::string);
void WriteHistogramsDICE(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::string);
std::string getOutputDirFromOutputFile(std::string);
bool getCascadeParticles(edm::Handle<std::vector<reco::GenParticle>>,int,std::string,unsigned int&,std::vector<reco::GenParticle>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&);


int main(int argc, char* argv[]) 
{
	gSystem->Load("libFWCoreFWLite.so");
	FWLiteEnabler::enable();
	
	////////////////////---------------------------------------------------------------------
	////////////////////----------------------------------------------
	// SET PARAMETERS //-----------------------
	////////////////////
	// first element must always be 'off'
	// make sure that the string and parameter inputs match up
	// nb: -1 on the cut parameter means no criteria to be met
	
	// step1: lepton/photon veto stage
	std::vector<std::string> step1labels = {"STEP1off", "STEP1on"};
	std::vector<int> step1_logicCut = {-1, 1};

	// step2: ht cut stage
	std::vector<std::string> step2labels = {"STEP2off", "STEP2on_HT800", "STEP2on_HT1000", "STEP2on_HT1200"};
	std::vector<int> step2_htCut = {-1, 800, 1000, 1200};

	// step3: ak8 jet pt cut stage (with |eta| < 2.4) 
	std::vector<std::string> step3labels = {"STEP3off", "STEP3on_AK8Leading170_AK8Secondary170", "STEP3AK8Leading200_AK8Secondary170"}
	std::vector<std::vector<int>> step3_ak8Cut = {{-1}, {170,170}, {200,170}}; // the sub-vectors correspond to the leading and secondary jet transverse momenta

	// step4: ak4 ISO jet cut stage (with |eta| < 3.0)
	std::vector<std::string> step4labels = {"STEP4off", "STEP4on_AK4Leading100_AK4Secondary100", "STEP4on_AK4Leading200_AK4Secondary100"}
	std::vector<std::vector<int>> step4_ak4Cut = {{-1}, {170,170}, {200,170}}; // the sub-vectors correspond to the leading and secondary jet transverse momenta

	////////////////////-----------------------
	////////////////////----------------------------------------------
	////////////////////---------------------------------------------------------------------

	// see if it is in runOnDice mode
	bool runOnDice = false;
    for (int i = 0; i < argc; ++i) {
        std::string argvString(argv[i]);
        if (argvString == "runOnDice") runOnDice = true;
    }

	// Create histograms, they are accessed by eg: h_["fatJetMass_loose"]->Fill(125.0);
	std::map<std::string, TH1F*> h_;
	std::map<std::string, TH2F*> h2_;
	CreateHistograms(h_, h2_);

	// Initialize command line parser
	optutl::CommandLineParser parser ("InvestigatePreSelection ");
	//////////////////
	// Set defaults //
	parser.integerValue ("maxevents"      ) = -1; // -1 for all events
	parser.integerValue ("outputevery"    ) = 100;
	parser.stringVector ("inputfiles"   ) = {"/hdfs/user/jt15104/Analysis_boostedNmssmHiggs/patTuples/CMSSW_8_0_21/signalSamples/nmssmSignalCascadeV05_13TeV_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/nmssmSignalCascadeV05_13TeV_patTupleAddBTag_ed8021v1_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/bTagPatTuple_10.root"};
	parser.stringValue  ("outputfile"     ) = "output_InvestigatePreSelection/histos.root";
	parser.boolValue    ("orderedsecondaryfiles") = false;
	//////////////////

	// Parse arguments_
	if (runOnDice) parser.parseArguments (argc-1, argv);
	else parser.parseArguments (argc, argv);
	int maxEvents_ = parser.integerValue("maxevents");
	unsigned int outputEvery_ = parser.integerValue("outputevery");
	std::vector<std::string> inputFiles_ = parser.stringVector("inputfiles");
	std::string outputFile_ = parser.stringValue("outputfile");
	bool justDoPlotting_ = parser.boolValue("orderedsecondaryfiles");

	if (runOnDice==false){ // don't do this if on DICE as it cannot do this command locally

		std::string outputDirectory_ = getOutputDirFromOutputFile(outputFile_);
		bool makeDir = !(std::system(Form("mkdir %s",outputDirectory_.c_str())));
		if (justDoPlotting_ == false && makeDir == false){
			std::cout << "The chosen output directory already exists, or the parent directory does not exist:" << std::endl;
			std::cout << "Do not wish to overwrite ROOT file: Exiting..." << std::endl;
			return 1;
		}
		
		// copy the code used to make the histogram ROOT file into the same directory (this could be out of sync if you edit after compilation)
		// also copy the parser values that were used to a .txt file (ie the input data used)
		if (justDoPlotting_ == false){
			
			std::system(Form("cp $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/bin/InvestigatePreSelection.cpp %s",outputDirectory_.c_str()));
			
			std::ofstream parserRecord;
			parserRecord.open(Form("%soptionsUsed.txt",outputDirectory_.c_str()));
			parserRecord << "MaxEvents: " << maxEvents_ << "\n";
			parserRecord << "InputFiles: " << "\n";
			for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
				parserRecord << inputFiles_[iFile] << "\n";
			}
			parserRecord.close();
		}
	} // closes 'if' runOnDice

	///////////////////////////////
	// Loop through the input files
	int ievt=0;
	// but we are just doing the plotting goto the end
	if (justDoPlotting_) goto plottingLabel;

	for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
	    
	    // Open input file
	    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
	    if( inFile ){
		
			fwlite::Event ev(inFile);

			// Loop through the events for this file
			for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
				edm::EventBase const & event = ev;

				// break loop if maximal number of events is reached 
				if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
				
				// event counter
				if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false){
					std::cout << "   File " << iFile+1 << " of " << inputFiles_.size() << ":";
					std::cout << "  processing event: " << ievt << std::endl;
				}

				// Handle to the ak4Jet collection
				edm::Handle<std::vector<pat::Jet>> ak4Jets;
				event.getByLabel(std::string("selectedPatJets"), ak4Jets); // old version

				// Handle to the fatJet collection
				edm::Handle<std::vector<pat::Jet>> fatJets;
				event.getByLabel(std::string("selectedPatJetsAK8PFCHS"), fatJets);

				// Handle to the Electron collection
				edm::Handle<std::vector<pat::Electron>> electrons;
				event.getByLabel(std::string("selectedPatElectrons"), electrons);

				// Handle to the Muon collection
				edm::Handle<std::vector<pat::Muon>> muons;
				event.getByLabel(std::string("selectedPatMuons"), muons);

				// Handle to the Photon collection
				edm::Handle<std::vector<pat::Photon>> photons;
				event.getByLabel(std::string("selectedPatPhotons"), photons);

				// ----------------//
				// ------------------------------------------------------------------------------------------------------------//
				// A N A L Y S I S
				// ------------------------------------------------------------------------------------------------------------//
				// ----------------//
				// ----------------//

				// calculate all the quantites that we need, eg ht
				for (size_t iJ = 0; iJ < ak4Jets->size(); ++iJ){
					const pat::Jet & ak4Jet = (*ak4Jets)[iJ];

					if (ak4Jet.pt() > 40.0 && fabs(ak4Jet.eta())<3.0){ // corresponds to the alpha_t HT definition
						ht = ht + ak4Jet.pt();
						//function: adds two pt vectors together in the (2d) transverse plane
						//output: first element is magnitude, second element is phi
						mhtVec = addTwoPtVectors(mhtVec[0], mhtVec[1], ak4Jet.pt(), ak4Jet.phi());
						
						indexBiasDeltaPhiJets.push_back(iJ); // use to index jets valid for biasedDeltaPhi and alphaT calculations
						epsilon = epsilon + ak4Jet.et(); // a component of the alphaT calculation
					}

					if (ak4Jet.pt() > leadingAk4JetPt){
						secondaryAk4JetPt = leadingAk4JetPt;
						leadingAk4JetPt = ak4Jet.pt();
						secondaryAk4JetEta = leadingAk4JetEta;
						leadingAk4JetEta = ak4Jet.eta();
					}
					else if (ak4Jet.pt() > secondaryAk4JetPt){
						secondaryAk4JetPt = ak4Jet.pt();
						secondaryAk4JetEta = ak4Jet.eta();
					}
				} // closes loop through the ak4 jets



				// loop through Step 1 states
				for (size_t iStep1 = 0; iStep1 < step1labels.size(); ++iStep1){
					
					// loop through Step 2 states if Step 1 is applied
					for (size_t iStep2 = 0; iStep1!=0 && iStep2 < step2labels.size(); ++iStep2){

						// loop through Step 3 states if Step 2 is applied
						for (size_t iStep3 = 0; iStep2!=0 && iStep3 < step3labels.size(); ++iStep3){

							// loop through Step 4 states if Step 3 is applied
							for (size_t iStep4 = 0; iStep3!=0 && iStep4 < step4labels.size(); ++iStep4){


								h_[Form("%s__%s__%s__%s__ht",step1labels[iStep1],step2labels[iStep2],step3labels[iStep3],step4labels[iStep4])]->Fill();


							} // closes loop through Step 4 states								
						} // closes loop through Step 3 states
					} // closes loop through Step 2 states
				} // closes loop through Step 1 states



				// -----------------------------
				// -----------------------------
				// E N D _ O F _ A N A L Y S I S
				// ------------------------------------------------------------------------------------------------------------//
				// ------------------------------------------------------------------------------------------------------------//

			} // closes loop through events for this file
			inFile->Close();
	    } // closes 'if' the file exists

	    // break loop if maximal number of events is reached:
	    // this has to be done twice to stop the file loop as well
	    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;

	} // closes loop through files

	// don't do the .pdf plotting on DICE (we will want to hadd outputs first)
	if (runOnDice==true){ 
		WriteHistogramsDICE(h_, h2_, outputFile_.c_str());
		return 0;
	}  

	WriteHistograms(h_, h2_, outputFile_.c_str());
	
	plottingLabel:
	// sketchy hack to get the higgsMass, squarkMass, ratio and mass splitting from the root file path name
	size_t endNameContainer = 8888;
	size_t beginNameContainer = 8888;
	std::string forwardSlash = "/";
	for (size_t c = inputFiles_[0].size()-2; c > 0; --c){
		
		if (inputFiles_[0][c] == forwardSlash[0]) endNameContainer = c; 
		if (inputFiles_[0].substr(c,2) == "mH"){
			beginNameContainer = c;
			break;
		}
	}
	std::string titleName = inputFiles_[0].substr(beginNameContainer,endNameContainer-beginNameContainer);
	std::cout << titleName << std::endl;

	if (std::system(Form("test -e %s",outputFile_.c_str())) == 0) // check the file exists
		PlottingMcSignalStudiesCMSSW createPlots(outputFile_.c_str(), titleName.c_str());
	else{
		std::cout << "the following output root file does not exist to make plots from:" << std::endl;
		std::cout << outputFile_ << std::endl;
	}

return 0;
} // closes the 'main' function







void CreateHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_)
{

	// fill the efficiency denominators
	for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinningD.size()-1; ++iEtaBin){

		if( fabs(higssBbGenParticleD.eta()) >= etaBinningD[iEtaBin] && fabs(higssBbGenParticleD.eta()) < etaBinningD[iEtaBin+1]){
			h_[Form("effDenominator_eta%.2f-%.2f", etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )]->Fill(higssBbGenParticleD.pt());
			h_[Form("effDenominator_eta%.2f-%.2f_fcnDR", etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )]->Fill(dRbb);

		} // closes 'if' eta within the set bin
	}




	// loop through Step 1 states
	for (size_t iStep1 = 0; iStep1 < step1labels.size(); ++iStep1){
		
		// loop through Step 2 states if Step 1 is applied
		for (size_t iStep2 = 0; iStep1!=0 && iStep2 < step2labels.size(); ++iStep2){

			// loop through Step 3 states if Step 2 is applied
			for (size_t iStep3 = 0; iStep2!=0 && iStep3 < step3labels.size(); ++iStep3){

				// loop through Step 4 states if Step 3 is applied
				for (size_t iStep4 = 0; iStep3!=0 && iStep4 < step4labels.size(); ++iStep4){


					h_[Form("%s__%s__%s__%s__ht",step1labels[iStep1],step2labels[iStep2],step3labels[iStep3],step4labels[iStep4])] = 				
					new TH1F(Form("%s__%s__%s__%s__ht",step1labels[iStep1],step2labels[iStep2],step3labels[iStep3],step4labels[iStep4]), ";detector H_{T} (GeV);a.u.", 100, 0, 7000);

					// h_[Form("%s__%s__%s__%s__ht",step1labels[iStep1],step2labels[iStep2],step3labels[iStep3],step4labels[iStep4])] = 				
					// new TH1F(Form("%s__%s__%s__%s__ht",step1labels[iStep1],step2labels[iStep2],step3labels[iStep3],step4labels[iStep4]), ";detector H_{T} (GeV);a.u.", 100, 0, 7000);

				} // closes loop through Step 4 states								
			} // closes loop through Step 3 states
		} // closes loop through Step 2 states
	} // closes loop through Step 1 states







} //closes the function 'CreateHistograms'







void WriteHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_, std::string outputFileD)
{
	TFile * outFile = new TFile(outputFileD.c_str(),"RECREATE");
	for ( auto & h : h_ ){
		h.second->Write();
	}
	for ( auto & h : h2_ ){
		h.second->Write();
	}
	outFile -> Close();
	delete outFile;
}







void WriteHistogramsDICE(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_, std::string outputFileD)
{
	TFile outFile(outputFileD.c_str(),"RECREATE");
	for ( auto & h : h_ ){
		h.second->Write();
	}
	for ( auto & h : h2_ ){
		h.second->Write();
	}
	outFile.Close();
	std::system("cp	*.root	../../.	"); // get the output back on soolin
}







std::string getOutputDirFromOutputFile(std::string outputFile)
{
	std::string forwardSlash = "/";
	std::string outputDirectory = outputFile;
	// strip the directory from the outputfile name
	for (size_t c = outputFile.size()-1; c > 0; --c){
		if (outputFile[c] == forwardSlash[0]){
			outputDirectory = outputFile.substr(0, c+1);
			break;
		}
	}
	return outputDirectory;
}