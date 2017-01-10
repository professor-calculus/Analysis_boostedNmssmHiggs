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


// do i need more of the above for met and stuff???
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

// Headers from this package
#include "Analysis/Analysis_boostedNmssmHiggs/interface/Kinematics.h"
// #include "Analysis/Analysis_boostedNmssmHiggs/interface/PlottingDoubleBTaggerEfficiencyStudies.h"


// TODO:
// TChain instead of TTree?
// New headers for MET and Gen and stuff?

				// the first element is leading arm (decay from highest pt squark)
				// the second element is for the secondary arm (decay from secondary pt squark)


// preliminary running, compile with scram b and then
// $ McSignalStudiesCMSSW inputfiles=XYZ outputfile=ABC orderedsecondaryfiles=0
// nb. if you are running it on DICE, include the word runOnDice at the end of the arguments of the executable

/* Notes on runOnDice mode
watchout...with this toggle the executable can now overwrite outputs!!!
hacks the original script abit so it can do some things slightly differently
it does so using the quantity: bool runOnDice
*/


void CreateHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&);
void WriteHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::string);
void WriteHistogramsDICE(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::string);
std::string getOutputDirFromOutputFile(std::string);
bool getCascadeParticles(edm::Handle<std::vector<reco::GenParticle>>,int,std::string,unsigned int&,std::vector<reco::GenParticle>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&);


int main(int argc, char* argv[]) 
{
	gSystem->Load("libFWCoreFWLite.so");
	FWLiteEnabler::enable();
	
	////////////////////
	////////////////////
	// SET PARAMETERS //
	std::vector<double> doubleBtagWP = {0.3, 0.6, 0.8, 0.9}; // these WP vectors must correspond to one-another
	std::vector<std::string> doubleBtagWPname = {"loose", "medium", "tight", "veryTight"};
	/////////////////////
	/////////////////////

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
	optutl::CommandLineParser parser ("McSignalStudiesCMSSW ");
	//////////////////
	// Set defaults //
	parser.integerValue ("maxevents"      ) = -1; // -1 for all events
	parser.integerValue ("outputevery"    ) = 100;
	parser.stringVector  ("inputfiles"    ) = {"/hdfs/user/jt15104/Analysis_boostedNmssmHiggs/patTuples/CMSSW_8_0_20/signalSamples/nmssmSignalCascadeV05_13TeV_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/nmssmSignalCascadeV05_13TeV_patTupleAddBTag_ed12_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/bTagPatTuple_888.root"};
	parser.stringValue  ("outputfile"     ) = "output_McSignalStudiesCMSSW/output_McSignalStudiesCMSSW.root";
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
			
			std::system(Form("cp $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/bin/McSignalStudiesCMSSW.cpp %s",outputDirectory_.c_str()));
			
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
				event.getByLabel(std::string("selectedPatJets"), ak4Jets);

				// Handle to the fatJet collection
				edm::Handle<std::vector<pat::Jet>> fatJets;
				event.getByLabel(std::string("selectedPatJetsAK8PFCHS"), fatJets);

				// Handle to the MET collection
				edm::Handle<std::vector<pat::MET>> METvec;
				event.getByLabel(std::string("patMETs"), METvec);

				// Handle to the Gen Particle collection
 				edm::Handle<std::vector<reco::GenParticle>> genParticles;
				event.getByLabel(std::string("genParticles"), genParticles);	


				// A N A L Y S I S
				// ------------------------------------------------------------------------------------------------------------//
				// first get all the gen particles in the cascade

				// for (size_t iJ = 0; iJ < ak4Jets->size(); ++iJ){
				// 	const pat::Jet & ak4Jet = (*ak4Jets)[iJ];
				// }

				// unsigned int numberOfFatJets = fatJets->size();
				// for (size_t iFJ = 0; iFJ < fatJets->size(); ++iFJ){
				// 	const pat::Jet & fatJet = (*fatJets)[iFJ];
				// }

				// for (size_t iM = 0; iM < METvec->size(); ++iM){
				// 	const pat::MET & MET = (*METvec)[iM];
				// }				


				// count how many gluinos are involved, we know we expect 2 squarks
				unsigned int gluinoCount = 0;
				// the first element is leading arm version of particle(decay from highest pt squark)
				// the second element is for the secondary arm version of particle (decay from secondary pt squark)
				std::vector<reco::GenParticle> squarkVec; // note that this object is slightly different to the others
				std::vector<const reco::Candidate*> qjetVec;
				std::vector<const reco::Candidate*> nlspVec;
				std::vector<const reco::Candidate*> lspVec;
				std::vector<const reco::Candidate*> higgsVec;
				std::vector<const reco::Candidate*> bVec;
				std::vector<const reco::Candidate*> bbarVec;

				bool gotCascadeParticles = getCascadeParticles(genParticles,ievt,inputFiles_[iFile],gluinoCount,squarkVec,qjetVec,nlspVec,lspVec,higgsVec,bVec,bbarVec);
				if (gotCascadeParticles == false) continue;

				std::cout << "leading squark pt: " << squarkVec[0].pt() << "   secondary squark pt: " << squarkVec[1].pt() << std::endl;
				std::cout << "leading lsp pt: " << lspVec[0]->pt() << "   second lsp pt: " << lspVec[1]->pt();
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
	// if (std::system(Form("test -e %s",outputFile_.c_str())) == 0) // check the file exists
	// 	PlottingDoubleBTaggerEfficiencyStudies createPlots(outputFile_.c_str(), doubleBtagWPname, etaBinning, massCut);
	// else{
	// 	std::cout << "the following output root file does not exist to make plots from:" << std::endl;
	// 	std::cout << outputFile_ << std::endl;
	// }

return 0;
} // closes the 'main' function








void CreateHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_)
{
	// create the debugging histograms
	h_["DEBUG_higgsBbDRpreMatching"] = new TH1F("DEBUG_higgsBbDRpreMatching", ";dR_bb;a.u.", 50, 0, 2.50);
	h_["DEBUG_higgsBbDRpreMatching_pt400to415_eta0to2p4"] = new TH1F("DEBUG_higgsBbDRpreMatching_pt400to415_eta0to2p4", ";dR_bb;a.u.", 100, 0, 2.50);
	h_["DEBUG_higgsBbDRpreMatching_pt450to465_eta0to2p4"] = new TH1F("DEBUG_higgsBbDRpreMatching_pt450to465_eta0to2p4", ";dR_bb;a.u.", 100, 0, 2.50);
	h_["DEBUG_higgsBbDRpreMatching_pt500to515_eta0to2p4"] = new TH1F("DEBUG_higgsBbDRpreMatching_pt500to515_eta0to2p4", ";dR_bb;a.u.", 100, 0, 2.50);
	h_["DEBUG_higgsBbDRpreMatching_3p400to415_eta0to2p4"] = new TH1F("DEBUG_higgsBbDRpreMatching_3p400to415_eta0to2p4", ";dR_bb;a.u.", 100, 0, 2.50);
	h_["DEBUG_higgsBbDRpreMatching_3p450to465_eta0to2p4"] = new TH1F("DEBUG_higgsBbDRpreMatching_3p450to465_eta0to2p4", ";dR_bb;a.u.", 100, 0, 2.50);
	h_["DEBUG_higgsBbDRpreMatching_3p500to515_eta0to2p4"] = new TH1F("DEBUG_higgsBbDRpreMatching_3p500to515_eta0to2p4", ";dR_bb;a.u.", 100, 0, 2.50);
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






// function to get the nmssm cascade particles, returns true if does so sucessfully:)
bool getCascadeParticles(edm::Handle<std::vector<reco::GenParticle>> genParticles, int ievt, std::string filename, unsigned int & gluinoCount, std::vector<reco::GenParticle> & squarkVec, std::vector<const reco::Candidate*> & qjetVec, std::vector<const reco::Candidate*> & nlspVec, std::vector<const reco::Candidate*> & lspVec, std::vector<const reco::Candidate*> & higgsVec, std::vector<const reco::Candidate*> & bVec, std::vector<const reco::Candidate*> & bbarVec)
{
	for (size_t iGen=0; iGen<genParticles->size(); ++iGen){
		const reco::GenParticle & genParticle = (*genParticles)[iGen];
		// get the squarks
		if ( abs(genParticle.pdgId()) == 1000021 ) gluinoCount++;
		if ( (abs(genParticle.pdgId()) == 1000001
		     || abs(genParticle.pdgId()) == 2000001
		     || abs(genParticle.pdgId()) == 1000002
		     || abs(genParticle.pdgId()) == 2000002
		     || abs(genParticle.pdgId()) == 1000003
		     || abs(genParticle.pdgId()) == 2000003
		     || abs(genParticle.pdgId()) == 1000004
		     || abs(genParticle.pdgId()) == 2000004) 
			 && genParticle.numberOfDaughters() == 2){

			int pdgId_0 = abs(genParticle.daughter(0)->pdgId());
			int pdgId_1 = abs(genParticle.daughter(1)->pdgId());

			if ( (pdgId_1 == 1 || pdgId_1 == 2 || pdgId_1 == 3 || pdgId_1 ==4) && pdgId_0 == 1000023 ){
				squarkVec.push_back(genParticle);
				qjetVec.push_back(genParticle.daughter(1));
				nlspVec.push_back(genParticle.daughter(0));		
			} 
			else if ( (pdgId_0 == 1 || pdgId_0 == 2 || pdgId_0 == 3 || pdgId_0 ==4) && pdgId_1 == 1000023 ){
				squarkVec.push_back(genParticle);					
				qjetVec.push_back(genParticle.daughter(0));
				nlspVec.push_back(genParticle.daughter(1));
			}
			else continue; // go look at the next candidate genParticle

			// nlsp decays
			// loop to ensure we get the decaying version of the nlsp
			while (nlspVec.back()->numberOfDaughters() != 0){
				if (nlspVec.back()->numberOfDaughters() == 2) break;
				else nlspVec.back() = nlspVec.back()->daughter(0);
			}

			pdgId_0 = abs(nlspVec.back()->daughter(0)->pdgId());
			pdgId_1 = abs(nlspVec.back()->daughter(1)->pdgId());
			if ( pdgId_1 == 35 && pdgId_0 == 1000022 ){
				higgsVec.push_back(nlspVec.back()->daughter(1));
				lspVec.push_back(nlspVec.back()->daughter(0));			
			} 
			else if ( pdgId_0 == 35 && pdgId_1 == 1000022 ){
				higgsVec.push_back(nlspVec.back()->daughter(0));
				lspVec.push_back(nlspVec.back()->daughter(1));
			}
			else {
				std::cerr << "ERROR: nlsp does not have the correct daughters" << std::endl;
				std::cerr << "File: " << filename << std::endl;
				std::cerr << "Event: " << ievt << std::endl;
				std::cerr << std::endl;
				return false;
			}

			// higgs decays
			if (higgsVec.back()->numberOfDaughters() != 2 ){
				std::cerr << "ERROR: higgs does not have two daughters" << std::endl;
				std::cerr << "File: " << filename << std::endl;
				std::cerr << "Event: " << ievt << std::endl;
				std::cerr << std::endl;
				return false;
			}
			pdgId_0 = higgsVec.back()->daughter(0)->pdgId();
			pdgId_1 = higgsVec.back()->daughter(1)->pdgId();
			if ( pdgId_0 == 5 && pdgId_1 == -5 ){
				bVec.push_back(higgsVec.back()->daughter(0));
				bbarVec.push_back(higgsVec.back()->daughter(1));
			}
			else if ( pdgId_1 == 5 && pdgId_0 == -5 ){
				bVec.push_back(higgsVec.back()->daughter(1));
				bbarVec.push_back(higgsVec.back()->daughter(0));			
			} 
			else {
				std::cerr << "ERROR: higgs does not have the correct daughters" << std::endl;
				std::cerr << "File: " << filename << std::endl;
				std::cerr << "Event: " << ievt << std::endl;
				std::cerr << std::endl;
				return false;
			}

		} // closes 'if' genParticle is a squark with two daughters
	} // closes loop through genParticles

	// check we have the right number of correctly decaying squarks
	if (squarkVec.size() != 2){
		std::cerr << "ERROR: haven't got two correctly decaying squarks" << std::endl;
		std::cerr << "File: " << filename << std::endl;
		std::cerr << "Event: " << ievt << std::endl;
		std::cerr << std::endl;
		return false;	
	}

	// ensure squarkVec entries are in pt order, and that decay arms follow this rule also
	if (squarkVec[1].pt() > squarkVec[0].pt()){
		std::swap(squarkVec[0],squarkVec[1]);
		std::swap(qjetVec[0],qjetVec[1]);
		std::swap(nlspVec[0],nlspVec[1]);
		std::swap(lspVec[0],lspVec[1]);
		std::swap(higgsVec[0],higgsVec[1]);
		std::swap(bVec[0],bVec[1]);
		std::swap(bbarVec[0],bbarVec[1]);
	}

	return true;
} // closes the funciton getCascadeParticles