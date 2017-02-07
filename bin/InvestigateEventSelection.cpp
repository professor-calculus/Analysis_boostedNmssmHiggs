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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

// Headers from this package
#include "Analysis/Analysis_boostedNmssmHiggs/interface/Kinematics.h"

// PAT TUPLE FORMAT ONLY !!!

/* TODO
no tau(track?) info for the lepton veto
*/

// preliminary running, compile with scram b and then
// $ InvestigateEventSelection inputfiles=XYZ outputfile=ABC orderedsecondaryfiles=0
// ($CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/InvestigateEventSelection/InvestigateEventSelection)
// nb. if you are running it on DICE, include the word runOnDice at the end of the arguments of the executable

/* Notes on runOnDice mode
watchout...with this toggle the executable can now overwrite outputs!!!
hacks the original script abit so it can do some things slightly differently
it does so using the quantity: bool runOnDice
*/

// NB: *known bug* if running on just plotting mode you still need to give it an inputfiles argument so that the script can get a title without seg faulting!!!

void CreateHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::vector<std::string>, std::vector<std::string>, std::vector<std::string>, std::vector<std::string>);
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
	// make sure that the string and parameter inputs match up <---DONT' FORGET
	// nb: -1 on the cut parameter means no criteria to be met
	
	// step1: lepton/photon veto stage
	std::vector<std::string> step1labels = {"STEP1off", "STEP1on"};

	// step2: ht cut stage
	std::vector<std::string> step2labels = {"STEP2off", "STEP2on_HT1000", "STEP2on_HT2000", "STEP2on_HT3000"};
	std::vector<int> step2_htCut = {-1, 1000, 2000, 3000};

	// step3: ak8 jet pt cut stage (with |eta| < eta_centralFatJet) 
	std::vector<std::string> step3labels = {"STEP3off", "STEP3on_AK8Leading200_AK8Secondary200", "STEP3on_AK8Leading300_AK8Secondary300", "STEP3on_AK8Leading400_AK8Secondary400"};
	std::vector<std::vector<int>> step3_ak8Cut = {{-1}, {200,200}, {300,300}, {400,400}}; // the sub-vectors correspond to the leading and secondary jet transverse momenta. 170 should always be absolute minimum.

	// step4: double b tag requirements on ak8 jets
	std::vector<std::string> step4labels = {"STEP4off", "STEP4on_looseDoubleBTags", "STEP4on_mediumDoubleBTags", "STEP4on_tightDoubleBTags"};
	std::vector<double> step4_doubleBTagWPs = {-999.99, 0.3, 0.6, 0.9}; // notOn, loose, medium and tight repsectively

	// additional parameters
	double eta_centralFatJet = 2.4;
	
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
	CreateHistograms(h_, h2_, step1labels, step2labels, step3labels, step4labels);

	// Initialize command line parser
	optutl::CommandLineParser parser ("InvestigateEventSelection ");
	//////////////////
	// Set defaults // (command line doesn't seem to override these options, in which case make sure they are commented out)
	parser.integerValue ("maxevents"      ) = -1; // -1 for all events
	parser.integerValue ("outputevery"    ) = 100;
	// parser.stringVector ("inputfiles"   ) = {"/hdfs/user/jt15104/Analysis_boostedNmssmHiggs/patTuples/CMSSW_8_0_21/signalSamples/nmssmSignalCascadeV05_13TeV_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/nmssmSignalCascadeV05_13TeV_patTupleAddBTag_ed8021v1_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/bTagPatTuple_10.root"};
	// parser.stringValue  ("outputfile"     ) = "output_InvestigateEventSelection/histos.root";
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
			
			std::system(Form("cp $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/bin/InvestigateEventSelection.cpp %s",outputDirectory_.c_str()));
			
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

				// Handle to the Tau collection
				edm::Handle<std::vector<pat::Tau>> taus;
				event.getByLabel(std::string("selectedPatTaus"), taus);

				// Handle to the Photon collection
				edm::Handle<std::vector<pat::Photon>> photons;
				event.getByLabel(std::string("selectedPatPhotons"), photons);

				// ----------------//
				// ------------------------------------------------------------------------------------------------------------//
				// A N A L Y S I S
				// ------------------------------------------------------------------------------------------------------------//
				// ----------------//
				// *************************************************
				// *************************************************
				// *************************************************
				// first calculate all the quantites that we need...
				
				// ak8 fatJets
				std::vector<double> leadingCentralFatJetPt(step4_doubleBTagWPs.size(), 0.0); // doubleBTag: no, loose, medium, tight 
				std::vector<double> secondaryCentralFatJetPt(step4_doubleBTagWPs.size(), 0.0);  // doubleBTag: no, loose, medium, tight

				for (size_t iFJ = 0; iFJ < fatJets->size(); ++iFJ){
					const pat::Jet & fatJet = (*fatJets)[iFJ];
					
					for (size_t iWP = 0; iWP < step4_doubleBTagWPs.size(); ++iWP){
						if  (fatJet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > step4_doubleBTagWPs[iWP] && fabs(fatJet.eta())<eta_centralFatJet){

							if (fatJet.pt() > leadingCentralFatJetPt[iWP]){
							
									secondaryCentralFatJetPt[iWP] = leadingCentralFatJetPt[iWP];
									leadingCentralFatJetPt[iWP] = fatJet.pt();
							}

							else if (fatJet.pt() > secondaryCentralFatJetPt[iWP]) secondaryCentralFatJetPt[iWP] = fatJet.pt();
						
						} // closes 'if' fat jet passes doubleBTag cut and is central
					} // closes loop through the doubleBTag WP's
				} // closes loop through the fat jets


				// ak4 jets and ht
				double ht = 0.0;
				for (size_t iJ = 0; iJ < ak4Jets->size(); ++iJ){
					const pat::Jet & ak4Jet = (*ak4Jets)[iJ];

					if (ak4Jet.pt() > 40.0 && fabs(ak4Jet.eta())<3.0) ht = ht + ak4Jet.pt();
				} // closes loop through the ak4 jets


				// leptons / photons
				bool eventPassLeptonPhotonVeto = true;
		
				for (size_t iE = 0; iE < electrons->size(); ++iE){
					
					const pat::Electron & electron = (*electrons)[iE];
					int electronIdAndIsolationValue = int(electron.electronID("eidLoose")) % 4;
					if (electronIdAndIsolationValue == 3 && electron.pt() > 10.0 && fabs(electron.eta()) < 2.5){
						eventPassLeptonPhotonVeto = false;
						break;
					}
				} // closes loop through the electrons

				for (size_t iM = 0; eventPassLeptonPhotonVeto && iM < muons->size(); ++iM){
					
					const pat::Muon & muon = (*muons)[iM];
					if (muon.isLooseMuon() && muon.pt() > 10.0 && fabs(muon.eta()) < 2.5 && muon.trackIso() < 3.0){
						eventPassLeptonPhotonVeto = false;
						break;
					}
				} // closes loop through muons

				// for (size_t iT = 0; eventPassLeptonPhotonVeto && iT < taus->size(); ++iT){
					
				// 	const pat::Tau & tau = (*taus)[iT];
				// } // closes loop through taus
				
				for (size_t iP = 0; eventPassLeptonPhotonVeto && iP < photons->size(); ++iP){
					
					const pat::Photon & photon = (*photons)[iP];
					if (photon.photonID("PhotonCutBasedIDLoose") && photon.pt() > 25.0 && fabs(photon.eta()) < 2.5 && photon.trackIso() < 3.0){
						eventPassLeptonPhotonVeto = false;
						break;
					}
				} // closes loop through photons			

				// *************************************************
				// *************************************************
				// *************************************************

				// loop through Step 1 states
				for (size_t iStep1 = 0; iStep1 < step1labels.size(); ++iStep1){
					
					// loop through Step 2 states if Step 1 is applied
					for (size_t iStep2 = 0; iStep2 < step2labels.size(); ++iStep2){

						// loop through Step 3 states if Step 2 is applied
						for (size_t iStep3 = 0; iStep3 < step3labels.size(); ++iStep3){

							// loop through Step 4 states if Step 3 is applied
							for (size_t iStep4 = 0; iStep4 < step4labels.size(); ++iStep4){

								if (iStep1 == 0 || eventPassLeptonPhotonVeto){
									if (iStep2 == 0 || ht > step2_htCut[iStep2]){
										if (iStep3 == 0 || (leadingCentralFatJetPt[0] > step3_ak8Cut[iStep3][0] && secondaryCentralFatJetPt[0] > step3_ak8Cut[iStep3][1]) ){
											if (iStep4 == 0 || (leadingCentralFatJetPt[iStep4] > step3_ak8Cut[iStep3][0] && secondaryCentralFatJetPt[iStep4] > step3_ak8Cut[iStep3][1]) ){

												h_[Form("ht__%s__%s__%s__%s",step1labels[iStep1].c_str(),step2labels[iStep2].c_str(),step3labels[iStep3].c_str(),step4labels[iStep4].c_str())]->Fill(ht);
												h_[Form("leadingCentralFatJetPt__%s__%s__%s__%s",step1labels[iStep1].c_str(),step2labels[iStep2].c_str(),step3labels[iStep3].c_str(),step4labels[iStep4].c_str())]->Fill(leadingCentralFatJetPt[iStep4]);
												h_[Form("secondaryCentralFatJetPt__%s__%s__%s__%s",step1labels[iStep1].c_str(),step2labels[iStep2].c_str(),step3labels[iStep3].c_str(),step4labels[iStep4].c_str())]->Fill(secondaryCentralFatJetPt[iStep4]);
												h2_[Form("leadingVsSecondaryCentralFatJetPt__%s__%s__%s__%s",step1labels[iStep1].c_str(),step2labels[iStep2].c_str(),step3labels[iStep3].c_str(),step4labels[iStep4].c_str())]->Fill(secondaryCentralFatJetPt[iStep4],leadingCentralFatJetPt[iStep4]);

											} // closes 'if' event passes Step4
										} // closes 'if' event passes Step3
									} // closes 'if' event passes Step2
								} // closes 'if' event passes Step1

								if (iStep3 == 0) break;
							} // closes loop through Step 4 states
							if (iStep2 == 0) break;
						} // closes loop through Step 3 states
						if (iStep1 == 0) break;
					} // closes loop through Step 2 states
				} // closes loop through Step 1 states

				// *************************************************
				// *************************************************
				// *************************************************
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

	// // // // // // // // // // // // // // // // // // //
	// // // // // // // // // // // // // // // // // // //
	// // // // // // // // // // // // // // // // // // //
	// make the table
	std::string outputDirectoryTable = "";
	if (justDoPlotting_){
	 	// get the name of directory holding the .root file, so we can save .pdfs here also
		for (size_t c = outputFile_.size()-1; c > 0; --c){
			std::string forwardSlash = "/";
			if (outputFile_[c] == forwardSlash[0]){
				outputDirectoryTable = outputFile_.substr(0, c+1);
				break;
			}
		}
	}
	else std::string outputDirectoryTable = getOutputDirFromOutputFile(outputFile_);
	// std::cout << outputDirectoryTable << std::endl;

	TFile * f = TFile::Open(outputFile_.c_str());

	int nCol = (step1labels.size()-1) * (step2labels.size()-1) * (step3labels.size()-1) * (step4labels.size()-1); // number of additional columns
	int divFactorStep2 = nCol / (step2labels.size()-1);
	int divFactorStep3 = nCol / ((step2labels.size()-1)*(step3labels.size()-1));

	std::ofstream table;
	table.open(Form("%sselectionTable.csv",outputDirectoryTable.c_str()));
	table << titleName << "\n";
	table << "Number of events surviving," << "\n";
	
	// Pre-Cuts
	table << "Pre-Cuts:," << "\n";
	TH1F * h0 = (TH1F*)f->Get(Form("ht__%s__%s__%s__%s",step1labels[0].c_str(),step2labels[0].c_str(),step3labels[0].c_str(),step4labels[0].c_str()));
	table << "," << h0->GetEntries() << "\n";
	
	// Step 1
	table << "Step 1 (lepton/photon veto):," << "\n";
	TH1F * h1 = (TH1F*)f->Get(Form("ht__%s__%s__%s__%s",step1labels[1].c_str(),step2labels[0].c_str(),step3labels[0].c_str(),step4labels[0].c_str()));
	table << "," << h1->GetEntries() << "\n";

	// Step 2
	table << "Step 2 (ht cut):,";
	for (int i = 0; i < nCol; ++i){
		if (i % divFactorStep2 == 0){
			int indexStep2 = i / divFactorStep2 + 1;
			table << step2_htCut[indexStep2] << " (GeV),";
		}
		else table << ",";
	}
	table << "\n";
	// step 2 values loop
	for (int i = 0; i < nCol; ++i){

		if (i % divFactorStep2 == 0){

			int indexStep2 = i / divFactorStep2 + 1;
			TH1F * h2 = (TH1F*)f->Get(Form("ht__%s__%s__%s__%s",step1labels[1].c_str(),step2labels[indexStep2].c_str(),step3labels[0].c_str(),step4labels[0].c_str()));
			table << "," << h2->GetEntries();
		}
		else table << ",";
	}
	table << ",\n";

	// Step 3
	table << "Step 3 (ak8 cut):,";
	for (int i = 0; i < nCol; ++i){
		if (i % divFactorStep3 == 0){
			int indexStep3 = ( (i/divFactorStep3) % (step3labels.size()-1) ) + 1;
			table << step3_ak8Cut[indexStep3][0] << ":" << step3_ak8Cut[indexStep3][1] << " (GeV),";
		}
		else table << ",";
	}
	table << "\n";
	// step 3 values loop
	for (int i = 0; i < nCol; ++i){

		if (i % divFactorStep3 == 0){

		int indexStep2 = i/divFactorStep2 + 1;
		int indexStep3 = ( (i/divFactorStep3) % (step3labels.size()-1) ) + 1;

		TH1F * h3 = (TH1F*)f->Get(Form("ht__%s__%s__%s__%s",step1labels[1].c_str(),step2labels[indexStep2].c_str(),step3labels[indexStep3].c_str(),step4labels[0].c_str()));
		table << "," << h3->GetEntries();
		}
		else table << ",";
	}
	table << ",\n";

	// Step 4
	table << "Step 4 (doubleBTag cut):,";
	for (int i = 0; i < nCol; ++i){

			int indexStep4 = (i%(step4labels.size()-1)) + 1;
			table << step4_doubleBTagWPs[indexStep4] << " WP,";
	}
	table << "\n";
	// step 4 values loop
	for (int i = 0; i < nCol; ++i){

		int indexStep2 = i/divFactorStep2 + 1;
		int indexStep3 = ( (i/divFactorStep3) % (step3labels.size()-1) ) + 1;
		int indexStep4 = (i%(step4labels.size()-1)) + 1;

		TH1F * h4 = (TH1F*)f->Get(Form("ht__%s__%s__%s__%s",step1labels[1].c_str(),step2labels[indexStep2].c_str(),step3labels[indexStep3].c_str(),step4labels[indexStep4].c_str()));
		table << "," << h4->GetEntries();		
	}
	table << ",\n";

	table.close();
	// // // // // // // // // // // // // // // // // // //
	// // // // // // // // // // // // // // // // // // //
	// // // // // // // // // // // // // // // // // // //
	return 0;
} // closes the 'main' function







void CreateHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_, std::vector<std::string> step1labels, std::vector<std::string> step2labels, std::vector<std::string> step3labels, std::vector<std::string> step4labels)
{
	// loop through Step 1 states
	for (size_t iStep1 = 0; iStep1 < step1labels.size(); ++iStep1){
		
		// loop through Step 2 states if Step 1 is applied
		for (size_t iStep2 = 0; iStep2 < step2labels.size(); ++iStep2){

			// loop through Step 3 states if Step 2 is applied
			for (size_t iStep3 = 0; iStep3 < step3labels.size(); ++iStep3){

				// loop through Step 4 states if Step 3 is applied
				for (size_t iStep4 = 0; iStep4 < step4labels.size(); ++iStep4){

					h_[Form("ht__%s__%s__%s__%s",step1labels[iStep1].c_str(),step2labels[iStep2].c_str(),step3labels[iStep3].c_str(),step4labels[iStep4].c_str())] = 				
					new TH1F(Form("ht__%s__%s__%s__%s",step1labels[iStep1].c_str(),step2labels[iStep2].c_str(),step3labels[iStep3].c_str(),step4labels[iStep4].c_str()), ";detector H_{T} (GeV);a.u.", 100, 0, 7000);

					h_[Form("leadingCentralFatJetPt__%s__%s__%s__%s",step1labels[iStep1].c_str(),step2labels[iStep2].c_str(),step3labels[iStep3].c_str(),step4labels[iStep4].c_str())] = 				
					new TH1F(Form("leadingCentralFatJetPt__%s__%s__%s__%s",step1labels[iStep1].c_str(),step2labels[iStep2].c_str(),step3labels[iStep3].c_str(),step4labels[iStep4].c_str()), ";leadingCentralFatJetPt (GeV);a.u.", 100, 0, 2500);

					h_[Form("secondaryCentralFatJetPt__%s__%s__%s__%s",step1labels[iStep1].c_str(),step2labels[iStep2].c_str(),step3labels[iStep3].c_str(),step4labels[iStep4].c_str())] = 				
					new TH1F(Form("secondaryCentralFatJetPt__%s__%s__%s__%s",step1labels[iStep1].c_str(),step2labels[iStep2].c_str(),step3labels[iStep3].c_str(),step4labels[iStep4].c_str()), ";secondaryCentralFatJetPt (GeV);a.u.", 100, 0, 2500);

					h2_[Form("leadingVsSecondaryCentralFatJetPt__%s__%s__%s__%s",step1labels[iStep1].c_str(),step2labels[iStep2].c_str(),step3labels[iStep3].c_str(),step4labels[iStep4].c_str())] = 				
					new TH2F(Form("leadingVsSecondaryCentralFatJetPt__%s__%s__%s__%s",step1labels[iStep1].c_str(),step2labels[iStep2].c_str(),step3labels[iStep3].c_str(),step4labels[iStep4].c_str()), ";secondaryCentralFatJetPt (GeV);leadingCentralDoubleBTagFatJetPt (GeV)", 200, 0, 2500, 200, 0, 2500);

					if (iStep3 == 0) break;
				} // closes loop through Step 4 states
				if (iStep2 == 0) break;
			} // closes loop through Step 3 states
			if (iStep1 == 0) break;
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