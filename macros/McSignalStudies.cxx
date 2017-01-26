// CPP headers
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h> 

// ROOT headers
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>
#include <TBranch.h>
#include <TClonesArray.h>

// DELPHES headers
#include "Analysis/Analysis_boostedNmssmHiggs/interface/DelphesClasses.h"
#include "Analysis/Analysis_boostedNmssmHiggs/interface/ExRootTreeReader.h"

// Headers from this package
#include "Analysis/Analysis_boostedNmssmHiggs/interface/Kinematics.h"
#include "Analysis/Analysis_boostedNmssmHiggs/interface/PlottingMcSignalStudies.h"

// *** NOTE ***
// need to have already loaded in a root session:
// gSystem->Load("libTreePlayer");
// gSystem->Load("$HOME/MG5_aMC_v2_3_3/Delphes/libDelphes.so");
// ************

// run with
// $ root -l -b -q 'McSignalStudies.cxx(std::vector<std::string> inputFiles_, std::string outputFile_)'
// eg.
// $ root -l -b -q 'McSignalStudies.cxx({"/storage/jt15104/madGraphProjects/testing/mH125p0_mSusy1000p0_ratio0p96_splitting2p0_10000events/dirA/dirB/dirC/tag_1_delphes_events.root"}, "/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/test/output.root")'
// 
// or use submit_McSignalStudies.cxx
// NOTE THAT THE INPUT MUST CONFORM TO THE mhXXXpX_mSusyXXXXpX_ratioXpXX_splittingXpX_XXXXXevents format

void CreateHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&);
void WriteHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::string);
bool indexAllCascadeParticles(std::vector<GenParticle*>,int,std::string,std::ofstream&,unsigned int&,std::vector<int>&,std::vector<int>&,std::vector<int>&,std::vector<int>&,std::vector<int>&,std::vector<int>&,std::vector<int>&);

void McSignalStudies(std::vector<std::string> inputFiles_, std::string outputFile_) 
{
	///////////////////////
	///////////////////////
	///////////////////////
	// Running Options ////
	int maxEvents_ = -1; // -1 for all events
	unsigned int outputEvery_ = 2000;
	// std::vector<std::string> inputFiles_ = {"/storage/jt15104/madGraphProjects/testing/mH125p0_mSusy1000p0_ratio0p96_splitting2p0_10000events/dirA/dirB/dirC/tag_1_delphes_events.root"};
	// std::string outputFile_ = "/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/test/output.root";
	bool justDoPlotting_ = false;
	// bool justDoPlotting_ = true;
	///////////////////////
	///////////////////////
	///////////////////////

	// Create histograms, they are accessed by eg: h_["fatJetMass_loose"]->Fill(125.0);
	std::map<std::string, TH1F*> h_;
	std::map<std::string, TH2F*> h2_;
	CreateHistograms(h_, h2_);

	std::string outputDirectory_;
	std::string forwardSlash = "/";
	// strip the directory from the outputfile name
	for (size_t c = outputFile_.size()-1; c > 0; --c){
		if (outputFile_[c] == forwardSlash[0]){
			outputDirectory_ = outputFile_.substr(0, c+1);
			break;
		}
	}
	bool makeDir = !(std::system(Form("mkdir %s",outputDirectory_.c_str())));
	if (justDoPlotting_ == false && makeDir == false){
		std::cout << "The chosen output directory already exists, or the parent directory does not exist:" << std::endl;
		std::cout << "Do not wish to overwrite ROOT file: Exiting..." << std::endl;
		return;
	}

	// copy the code used to make the histogram ROOT file into the same directory (this could be out of sync if you edit after compilation)
	if (justDoPlotting_ == false) std::system(Form("cp $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/macros/McSignalStudies.cxx %s",outputDirectory_.c_str()));		

	// create an error log, for events that don't have the particles we expect... 
	std::ofstream errorRecord;
	if (justDoPlotting_ == false) errorRecord.open(Form("%serrors.txt",outputDirectory_.c_str()));
	
	///////////////////////////////
	// Loop through the input files
	///////////////////////////////
	int ievt=0;
	/////////////////////////////////////////////////
	// if we are just doing the plotting goto the end
	if (justDoPlotting_) goto plottingLabel;
	/////////////////////////////////////////////////

	for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
	    
	    // Open input file
	    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
	    if( inFile ){
		
			TTree * t = (TTree*)inFile->Get("Delphes");
			ExRootTreeReader *reader = new ExRootTreeReader(t);
			TClonesArray *branchParticle = reader->UseBranch("Particle");
			TClonesArray *branchMET = reader->UseBranch("MissingET");
			TClonesArray *branchHT = reader->UseBranch("ScalarHT");
			TClonesArray *branchJet = reader->UseBranch("Jet");

		    Int_t nevents = (Int_t)t->GetEntries();
		    for (Int_t ievt=0; ievt<nevents; ievt++){  

				// break loop if maximal number of events is reached 
				if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
					
				// event counter
				if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false){
					std::cout << "   File " << iFile+1 << " of " << inputFiles_.size() << ":";
					std::cout << "  processing event: " << ievt << " of " << nevents << std::endl;
				}

				// load info for this event
		        reader->ReadEntry(ievt);

				// load the gen particle info
				std::vector<GenParticle*> particleVec;
				for (int iPar=0; iPar<branchParticle->GetEntries(); ++iPar){
				    GenParticle * particle = (GenParticle*) branchParticle->At(iPar);
				    particleVec.push_back(particle);
				    // std::cout << particleVec[iPar]->PT << std::endl; // an example of how to read info from this object we created
				}

				// load the jet info
				std::vector<Jet*> jetVec;
				for (int iJet=0; iJet<branchJet->GetEntries(); ++iJet){
				    Jet * jet = (Jet*) branchJet->At(iJet);
				    jetVec.push_back(jet);
				    // std::cout << particleVec[iPar]->PT << std::endl; // an example of how to read info from this object we created
				}

				// load the met info
				if (branchMET->GetEntries() != 1){
					std::cout << "We do not have ONE met quantity: Exiting..." << std::endl;
					return;
				}
				MissingET* metDetector = (MissingET*) branchMET->At(0);

				// load the ht info NOW WE WORK IT OUT FOR OURSELVES
				// if (branchHT->GetEntries() != 1){
				// 	std::cout << "We do not have ONE ht quantity: Exiting..." << std::endl;
				// 	return;
				// }
				// ScalarHT* htDetector = (ScalarHT*) branchHT->At(0);

				// ------------------------------------------------------------------------------------------------------------//
				// A N A L Y S I S
				// ------------------------------------------------------------------------------------------------------------//
				// ------------------------------------------------------------------------------------------------------------//

				// count how many gluinos are involved, we know we expect 2 squarks
				unsigned int gluinoCount = 0;
				// the first element is leading arm (decay from highest pt squark)
				// the second element is for the secondary arm (decay from secondary pt squark)
				// access leading lsp PT like so: particleVec[lspIndices[0]]->PT
				// access secondary lsp PT like so: particleVec[lspIndices[1]]->PT
				std::vector<int> squarkIndices;
				std::vector<int> qjetIndices;
				std::vector<int> nlspIndices;
				std::vector<int> lspIndices;
				std::vector<int> higgsIndices;
				std::vector<int> bIndices;
				std::vector<int> bbarIndices;
				bool allParticlesPresent = indexAllCascadeParticles(particleVec,ievt,inputFiles_[iFile],errorRecord,gluinoCount,squarkIndices,qjetIndices,nlspIndices,lspIndices,higgsIndices,bIndices,bbarIndices);
				if (allParticlesPresent==false) continue;


				// Fill histograms with properties of these particles:)
				h_["numberOfGluinos"]->Fill(gluinoCount);
				
				h_["leadingSquarkPt"]->Fill(particleVec[squarkIndices[0]]->PT);
				if (gluinoCount == 0) h_["leadingSquarkPt_zeroGluinos"]->Fill(particleVec[squarkIndices[0]]->PT);
				if (gluinoCount == 1) h_["leadingSquarkPt_oneGluinos"]->Fill(particleVec[squarkIndices[0]]->PT);
				if (gluinoCount == 2) h_["leadingSquarkPt_twoGluinos"]->Fill(particleVec[squarkIndices[0]]->PT);

				h_["leadingSquarkEta"]->Fill(particleVec[squarkIndices[0]]->Eta);
				if (gluinoCount == 0) h_["leadingSquarkEta_zeroGluinos"]->Fill(particleVec[squarkIndices[0]]->Eta);
				if (gluinoCount == 1) h_["leadingSquarkEta_oneGluinos"]->Fill(particleVec[squarkIndices[0]]->Eta);
				if (gluinoCount == 2) h_["leadingSquarkEta_twoGluinos"]->Fill(particleVec[squarkIndices[0]]->Eta);

				h_["secondarySquarkPt"]->Fill(particleVec[squarkIndices[1]]->PT);
				if (gluinoCount == 0) h_["secondarySquarkPt_zeroGluinos"]->Fill(particleVec[squarkIndices[1]]->PT);
				if (gluinoCount == 1) h_["secondarySquarkPt_oneGluinos"]->Fill(particleVec[squarkIndices[1]]->PT);
				if (gluinoCount == 2) h_["secondarySquarkPt_twoGluinos"]->Fill(particleVec[squarkIndices[1]]->PT);

				h_["secondarySquarkEta"]->Fill(particleVec[squarkIndices[1]]->Eta);
				if (gluinoCount == 0) h_["secondarySquarkEta_zeroGluinos"]->Fill(particleVec[squarkIndices[1]]->Eta);
				if (gluinoCount == 1) h_["secondarySquarkEta_oneGluinos"]->Fill(particleVec[squarkIndices[1]]->Eta);
				if (gluinoCount == 2) h_["secondarySquarkEta_twoGluinos"]->Fill(particleVec[squarkIndices[1]]->Eta);

				h2_["leadingSquarkPt_SecondarySquarkPt"]->Fill(particleVec[squarkIndices[1]]->PT,particleVec[squarkIndices[0]]->PT);
				h2_["leadingSquarkEta_SecondarySquarkEta"]->Fill(particleVec[squarkIndices[1]]->Eta,particleVec[squarkIndices[0]]->Eta);
				h2_["leadingSquarkPhi_SecondarySquarkPhi"]->Fill(particleVec[squarkIndices[1]]->Phi,particleVec[squarkIndices[0]]->Phi);

				h_["leadingQjetPt"]->Fill(particleVec[qjetIndices[0]]->PT);
				h_["leadingQjetEta"]->Fill(particleVec[qjetIndices[0]]->Eta);
				h_["secondaryQjetPt"]->Fill(particleVec[qjetIndices[1]]->PT);
				h_["secondaryQjetEta"]->Fill(particleVec[qjetIndices[1]]->Eta);
				h2_["leadingQjetPt_secondaryQjetPt"]->Fill(particleVec[qjetIndices[1]]->PT, particleVec[qjetIndices[0]]->PT);
				h2_["leadingQjetEta_secondaryQjetEta"]->Fill(particleVec[qjetIndices[1]]->Eta, particleVec[qjetIndices[0]]->Eta);
				h2_["leadingQjetPhi_secondaryQjetPhi"]->Fill(particleVec[qjetIndices[1]]->Phi, particleVec[qjetIndices[0]]->Phi);

				h_["leadingNlspPt"]->Fill(particleVec[nlspIndices[0]]->PT);
				h_["leadingNlspEta"]->Fill(particleVec[nlspIndices[0]]->Eta);
				h_["secondaryNlspPt"]->Fill(particleVec[nlspIndices[1]]->PT);
				h_["secondaryNlspEta"]->Fill(particleVec[nlspIndices[1]]->Eta);

				h2_["leadingQjetPt_leadingHiggsPt"]->Fill(particleVec[higgsIndices[0]]->PT, particleVec[qjetIndices[0]]->PT);
				h2_["leadingQjetEta_leadingHiggsEta"]->Fill(particleVec[higgsIndices[0]]->Eta, particleVec[qjetIndices[0]]->Eta);
				h2_["leadingQjetPhi_leadingHiggsPhi"]->Fill(particleVec[higgsIndices[0]]->Phi, particleVec[qjetIndices[0]]->Phi);
				h2_["secondaryQjetPt_secondaryHiggsPt"]->Fill(particleVec[higgsIndices[1]]->PT, particleVec[qjetIndices[1]]->PT);
				h2_["secondaryQjetEta_secondaryHiggsEta"]->Fill(particleVec[higgsIndices[1]]->Eta, particleVec[qjetIndices[1]]->Eta);
				h2_["secondaryQjetPhi_secondaryHiggsPhi"]->Fill(particleVec[higgsIndices[1]]->Phi, particleVec[qjetIndices[1]]->Phi);

				h_["leadingHiggsPt"]->Fill(particleVec[higgsIndices[0]]->PT);
				h_["leadingHiggsEta"]->Fill(particleVec[higgsIndices[0]]->Eta);
				h_["secondaryHiggsPt"]->Fill(particleVec[higgsIndices[1]]->PT);
				h_["secondaryHiggsEta"]->Fill(particleVec[higgsIndices[1]]->Eta);
				h2_["leadingHiggsPt_secondaryHiggsPt"]->Fill(particleVec[higgsIndices[1]]->PT, particleVec[higgsIndices[0]]->PT);
				h2_["leadingHiggsEta_secondaryHiggsEta"]->Fill(particleVec[higgsIndices[1]]->Eta, particleVec[higgsIndices[0]]->Eta);
				h2_["leadingHiggsPhi_secondaryHiggsPhi"]->Fill(particleVec[higgsIndices[1]]->Phi, particleVec[higgsIndices[0]]->Phi);
				double dphiLeadingHiggsQjet = delPhi(particleVec[higgsIndices[0]]->Phi,particleVec[qjetIndices[0]]->Phi);
				if (dphiLeadingHiggsQjet < 0) dphiLeadingHiggsQjet += 2*M_PI;
				h_["leadingHiggsQjetDphi"]->Fill(dphiLeadingHiggsQjet);
				double dphiSecondaryHiggsQjet = delPhi(particleVec[higgsIndices[1]]->Phi,particleVec[qjetIndices[1]]->Phi);
				if (dphiSecondaryHiggsQjet < 0) dphiSecondaryHiggsQjet += 2*M_PI;
				h_["secondaryHiggsQjetDphi"]->Fill(dphiSecondaryHiggsQjet);
				double dphiLeadingHiggsSecondaryHiggs = delPhi(particleVec[higgsIndices[0]]->Phi,particleVec[higgsIndices[1]]->Phi);
				if (dphiLeadingHiggsSecondaryHiggs < 0) dphiLeadingHiggsSecondaryHiggs += 2*M_PI;
				h_["leadingHiggsSecondaryHiggsDphi"]->Fill(dphiLeadingHiggsSecondaryHiggs);
				double dphiLeadingQjetSecondaryQjet = delPhi(particleVec[qjetIndices[0]]->Phi,particleVec[qjetIndices[1]]->Phi);
				if (dphiLeadingQjetSecondaryQjet < 0) dphiLeadingQjetSecondaryQjet += 2*M_PI;
				h_["leadingQjetSecondaryQjetDphi"]->Fill(dphiLeadingQjetSecondaryQjet);

				h_["leadingLspPt"]->Fill(particleVec[lspIndices[0]]->PT);
				h_["leadingLspEta"]->Fill(particleVec[lspIndices[0]]->Eta);
				h_["secondaryLspPt"]->Fill(particleVec[lspIndices[1]]->PT);
				h_["secondaryLspEta"]->Fill(particleVec[lspIndices[1]]->Eta);
				//function: adds two pt vectors together in the (2d) transverse plane
				//output: first element is magnitude, second element is phi
				std::vector<double> lspMet = addTwoPtVectors(particleVec[lspIndices[0]]->PT, particleVec[lspIndices[0]]->Phi, particleVec[lspIndices[1]]->PT, particleVec[lspIndices[1]]->Phi);
				h_["lspMET"]->Fill(lspMet[0]);

				double leadingBBbarInvmass = invMass_v1(particleVec[bIndices[0]]->E, particleVec[bbarIndices[0]]->E, particleVec[bIndices[0]]->Px, particleVec[bbarIndices[0]]->Px, particleVec[bIndices[0]]->Py, particleVec[bbarIndices[0]]->Py, particleVec[bIndices[0]]->Pz, particleVec[bbarIndices[0]]->Pz, particleVec[bIndices[0]]->Mass, particleVec[bbarIndices[0]]->Mass);
				h_["leadingBBbarInvmass"]->Fill(leadingBBbarInvmass);
				double leadingBBbarSeperation = delR( delPhi(particleVec[bIndices[0]]->Phi, particleVec[bbarIndices[0]]->Phi), delEta(particleVec[bIndices[0]]->Eta, particleVec[bbarIndices[0]]->Eta) );
				h_["leadingBBbarSeperation"]->Fill(leadingBBbarSeperation);
				double secondaryBBbarInvmass = invMass_v1(particleVec[bIndices[1]]->E, particleVec[bbarIndices[1]]->E, particleVec[bIndices[1]]->Px, particleVec[bbarIndices[1]]->Px, particleVec[bIndices[1]]->Py, particleVec[bbarIndices[1]]->Py, particleVec[bIndices[1]]->Pz, particleVec[bbarIndices[1]]->Pz, particleVec[bIndices[1]]->Mass, particleVec[bbarIndices[1]]->Mass);
				h_["secondaryBBbarInvmass"]->Fill(secondaryBBbarInvmass);
				double secondaryBBbarSeperation = delR( delPhi(particleVec[bIndices[1]]->Phi, particleVec[bbarIndices[1]]->Phi), delEta(particleVec[bIndices[1]]->Eta, particleVec[bbarIndices[1]]->Eta) );
				h_["secondaryBBbarSeperation"]->Fill(secondaryBBbarSeperation);
				h2_["leadingBBbarSeperation_secondaryBBbarSeperation"]->Fill(secondaryBBbarSeperation,leadingBBbarSeperation);

				h2_["leadingBBbarSeperation_massHiggsOverPt"]->Fill(particleVec[higgsIndices[0]]->Mass / (particleVec[higgsIndices[0]]->PT) , leadingBBbarSeperation);
				h2_["secondaryBBbarSeperation_massHiggsOverPt"]->Fill(particleVec[higgsIndices[1]]->Mass / (particleVec[higgsIndices[1]]->PT), secondaryBBbarSeperation);
				double higgs3mom = particleVec[higgsIndices[0]]->PT * cosh(particleVec[higgsIndices[0]]->Eta);
				h2_["leadingBBbarSeperation_massHiggsOverP3"]->Fill(particleVec[higgsIndices[0]]->Mass / higgs3mom , leadingBBbarSeperation);
				higgs3mom = particleVec[higgsIndices[1]]->PT * cosh(particleVec[higgsIndices[1]]->Eta);
				h2_["secondaryBBbarSeperation_massHiggsOverP3"]->Fill(particleVec[higgsIndices[1]]->Mass / higgs3mom, secondaryBBbarSeperation);

				if (leadingBBbarSeperation <= 0.8 && secondaryBBbarSeperation <= 0.8) h2_["logicBbDr"]->Fill(0.1,0.1);
				if (leadingBBbarSeperation <= 0.8 && secondaryBBbarSeperation > 0.8) h2_["logicBbDr"]->Fill(1.1,0.1);
				if (leadingBBbarSeperation > 0.8 && secondaryBBbarSeperation <= 0.8) h2_["logicBbDr"]->Fill(0.1,1.1);
				if (leadingBBbarSeperation > 0.8 && secondaryBBbarSeperation > 0.8) h2_["logicBbDr"]->Fill(1.1,1.1);

				if (particleVec[higgsIndices[0]]->PT <= 170 && particleVec[higgsIndices[1]]->PT <= 170) h2_["logicHiggsPt"]->Fill(0.1,0.1);
				if (particleVec[higgsIndices[0]]->PT <= 170 && particleVec[higgsIndices[1]]->PT > 170) h2_["logicHiggsPt"]->Fill(1.1,0.1);
				if (particleVec[higgsIndices[0]]->PT > 170 && particleVec[higgsIndices[1]]->PT <= 170) h2_["logicHiggsPt"]->Fill(0.1,1.1);
				if (particleVec[higgsIndices[0]]->PT > 170 && particleVec[higgsIndices[1]]->PT > 170) h2_["logicHiggsPt"]->Fill(1.1,1.1);

				if (particleVec[higgsIndices[0]]->PT <= 170 && leadingBBbarSeperation <= 0.8) h2_["logicBbDrHiggsPt"]->Fill(0.1,0.1);
				if (particleVec[higgsIndices[0]]->PT > 170 && leadingBBbarSeperation <= 0.8) h2_["logicBbDrHiggsPt"]->Fill(1.1,0.1);
				if (particleVec[higgsIndices[0]]->PT <= 170 && leadingBBbarSeperation > 0.8) h2_["logicBbDrHiggsPt"]->Fill(0.1,1.1);
				if (particleVec[higgsIndices[0]]->PT > 170 && leadingBBbarSeperation > 0.8) h2_["logicBbDrHiggsPt"]->Fill(1.1,1.1);
				if (particleVec[higgsIndices[1]]->PT <= 170 && secondaryBBbarSeperation <= 0.8) h2_["logicBbDrHiggsPt"]->Fill(0.1,0.1);
				if (particleVec[higgsIndices[1]]->PT > 170 && secondaryBBbarSeperation <= 0.8) h2_["logicBbDrHiggsPt"]->Fill(1.1,0.1);
				if (particleVec[higgsIndices[1]]->PT <= 170 && secondaryBBbarSeperation > 0.8) h2_["logicBbDrHiggsPt"]->Fill(0.1,1.1);
				if (particleVec[higgsIndices[1]]->PT > 170 && secondaryBBbarSeperation > 0.8) h2_["logicBbDrHiggsPt"]->Fill(1.1,1.1);

				double numHiggsHighPt = 0.0;
				double numHighBbSep = 0.0;
				if (particleVec[higgsIndices[0]]->PT > 170) numHiggsHighPt++;
				if (particleVec[higgsIndices[1]]->PT > 170) numHiggsHighPt++;
				if (leadingBBbarSeperation > 0.8) numHighBbSep++;
				if (secondaryBBbarSeperation > 0.8) numHighBbSep++;
				numHiggsHighPt += 0.01;
				numHighBbSep += 0.01;
				h2_["logicNumberHiggsPtNumberBbDr"]->Fill(numHiggsHighPt,numHighBbSep);


				// Detector Plots //

				// MET
				h_["detectorMET"]->Fill(metDetector->MET);

				// ak4 Jets (and HT, MHT & SUSY variables) nb: they are not in pt order
				unsigned int numAk4JetsOver200GeV = 0;
				unsigned int numAk4JetsOver20GeV = 0;
				unsigned int numAk4JetsOver40GeV = 0;
				unsigned int numAk4JetsOver100GeV = 0;
				unsigned int numCentralAk4JetsOver200GeV = 0;
				unsigned int numCentralAk4JetsOver20GeV = 0;
				unsigned int numCentralAk4JetsOver40GeV = 0;
				unsigned int numCentralAk4JetsOver100GeV = 0;
				unsigned int numForwardAk4JetsOver200GeV = 0;
				unsigned int numForwardAk4JetsOver20GeV = 0;
				unsigned int numForwardAk4JetsOver40GeV = 0;
				unsigned int numForwardAk4JetsOver100GeV = 0;
				double ht = 0.0;
				std::vector<double> mhtVec = {0,0}; // first element magnitude, second element phi
				std::vector<int> indexBiasDeltaPhiJets; // use to index jets valid for biasedDeltaPhi and alphaT calculations
				double epsilon = 0.0; // a component of the alphaT calculation
				double leadingAk4JetPt = 0.0;
				double leadingAk4JetEta = 999.99;
				double secondaryAk4JetPt = 0.0;
				double secondaryAk4JetEta = 999.99;

				for (size_t iJ = 0; iJ < jetVec.size(); ++iJ){
					
					if (jetVec[iJ]->PT > 200.0) numAk4JetsOver200GeV++;
					if (jetVec[iJ]->PT > 20.0) numAk4JetsOver20GeV++;
					if (jetVec[iJ]->PT > 40.0) numAk4JetsOver40GeV++;
					if (jetVec[iJ]->PT > 100.0) numAk4JetsOver100GeV++;
					if (jetVec[iJ]->PT > 200.0 && fabs(jetVec[iJ]->Eta) <= 3.0) numCentralAk4JetsOver200GeV++;
					if (jetVec[iJ]->PT > 20.0 && fabs(jetVec[iJ]->Eta) <= 3.0) numCentralAk4JetsOver20GeV++;
					if (jetVec[iJ]->PT > 40.0 && fabs(jetVec[iJ]->Eta) <= 3.0) numCentralAk4JetsOver40GeV++;
					if (jetVec[iJ]->PT > 100.0 && fabs(jetVec[iJ]->Eta) <= 3.0) numCentralAk4JetsOver100GeV++;
					if (jetVec[iJ]->PT > 200.0 && fabs(jetVec[iJ]->Eta) > 3.0) numForwardAk4JetsOver200GeV++;
					if (jetVec[iJ]->PT > 20.0 && fabs(jetVec[iJ]->Eta) > 3.0) numForwardAk4JetsOver20GeV++;
					if (jetVec[iJ]->PT > 40.0 && fabs(jetVec[iJ]->Eta) > 3.0) numForwardAk4JetsOver40GeV++;
					if (jetVec[iJ]->PT > 100.0 && fabs(jetVec[iJ]->Eta) > 3.0) numForwardAk4JetsOver100GeV++;

					if (jetVec[iJ]->PT > 40.0 && abs(jetVec[iJ]->Eta)<3.0){ // corresponds to the alpha_t HT definition
						ht = ht + jetVec[iJ]->PT;
						//function: adds two pt vectors together in the (2d) transverse plane
						//output: first element is magnitude, second element is phi
						mhtVec = addTwoPtVectors(mhtVec[0], mhtVec[1], jetVec[iJ]->PT, jetVec[iJ]->Phi);
						
						indexBiasDeltaPhiJets.push_back(iJ); // use to index jets valid for biasedDeltaPhi and alphaT calculations
						
						double jetEt = sqrt( jetVec[iJ]->PT * jetVec[iJ]->PT + (jetVec[iJ]->Mass * jetVec[iJ]->Mass)/(1 + sinh(jetVec[iJ]->Eta) * sinh(jetVec[iJ]->Eta) ) );
						epsilon = epsilon + jetEt; // a component of the alphaT calculation
					}

					if (jetVec[iJ]->PT > leadingAk4JetPt){
						secondaryAk4JetPt = leadingAk4JetPt;
						leadingAk4JetPt = jetVec[iJ]->PT;
						secondaryAk4JetEta = leadingAk4JetEta;
						leadingAk4JetEta = jetVec[iJ]->Eta;
					}
					else if (jetVec[iJ]->PT > secondaryAk4JetPt){
						secondaryAk4JetPt = jetVec[iJ]->PT;
						secondaryAk4JetEta = jetVec[iJ]->Eta;
					}
				} // closes loop through the ak4 jets

				// *** biased delta phi calc ***
				std::vector<double> biasedDeltaPhiVec;
				for (size_t i = 0; i < indexBiasDeltaPhiJets.size(); ++i){ // i indexes the solo jet for a given configuration of biasedDeltaPhi

					std::vector<double> biasedMhtVec = {0,0}; // first element magnitude, second element phi
					for (size_t j = 0; j < indexBiasDeltaPhiJets.size(); ++j){ // j indexes all the non solo jets for a given configuration of biasedDeltaPhi
						if (j != i) biasedMhtVec = addTwoPtVectors(biasedMhtVec[0], biasedMhtVec[1], jetVec[indexBiasDeltaPhiJets[j]]->PT, -1*jetVec[indexBiasDeltaPhiJets[j]]->Phi);
					}
					double biasedDeltaPhi = delPhi(jetVec[indexBiasDeltaPhiJets[i]]->Phi,biasedMhtVec[1]);
					if (biasedDeltaPhi > M_PI) biasedDeltaPhi -= 2*M_PI;
					if (biasedDeltaPhi < -1 * M_PI) biasedDeltaPhi += 2*M_PI;
					biasedDeltaPhiVec.push_back(fabs(biasedDeltaPhi));
				}
				int minindex = min_element(biasedDeltaPhiVec.begin(), biasedDeltaPhiVec.end()) - biasedDeltaPhiVec.begin();
				double minBiasedDeltaPhi = biasedDeltaPhiVec[minindex];

				// *** alphaT calc ***
				double alphaT = 999.99; // set to a large value so the event default is to pass cut
				if (indexBiasDeltaPhiJets.size() > 1){
					// *** delta epsilon calc *** (part of the general case alpha_t calculation)
					std::vector<int> bitOnOrOff(indexBiasDeltaPhiJets.size(),1); 
					std::vector<int> sameBitCount(indexBiasDeltaPhiJets.size(),0);
					std::vector<double> deltaEpsilonVec;

					for (unsigned int i = 0; i < pow(2,indexBiasDeltaPhiJets.size()) / 2; ++i){

						double pseudoJetOneEnergy = 0.0;
						double pseudoJetTwoEnergy = 0.0;

						for (size_t j = 0; j < indexBiasDeltaPhiJets.size(); ++j){

							double jetEt = sqrt( jetVec[j]->PT * jetVec[j]->PT + (jetVec[j]->Mass * jetVec[j]->Mass)/(1 + sinh(jetVec[j]->Eta) * sinh(jetVec[j]->Eta) ) );
							if (bitOnOrOff[j] == 1) pseudoJetOneEnergy += jetEt;
							if (bitOnOrOff[j] == 0) pseudoJetTwoEnergy += jetEt;

							sameBitCount[j] = sameBitCount[j]+1;
							if (sameBitCount[j] == pow(2,j)){

								sameBitCount[j] = 0;
								if (bitOnOrOff[j] == 0) bitOnOrOff[j] = 1; 
								else if (bitOnOrOff[j] == 1) bitOnOrOff[j] = 0;
							}
						} // closes loop through jet number index

						deltaEpsilonVec.push_back(fabs(pseudoJetOneEnergy - pseudoJetTwoEnergy));
					} // closes loop through 2^(numJets) / 2 iterations for delta epsilon calc

					int minindex2 = min_element(deltaEpsilonVec.begin(), deltaEpsilonVec.end()) - deltaEpsilonVec.begin();
					double minDeltaEpsilon = deltaEpsilonVec[minindex2];

					alphaT = 0.5 * (epsilon - minDeltaEpsilon) / (sqrt(epsilon * epsilon - mhtVec[0] * mhtVec[0]));
				} // closes 'if' more than one valid jet for the calculation
				// *** end of alphaT calc ***
				
				h_["detectorNumAk4JetsOver200GeV"]->Fill(numAk4JetsOver200GeV);
				h_["detectorNumAk4JetsOver20GeV"]->Fill(numAk4JetsOver20GeV);
				h_["detectorNumAk4JetsOver40GeV"]->Fill(numAk4JetsOver40GeV);
				h_["detectorNumAk4JetsOver100GeV"]->Fill(numAk4JetsOver100GeV);
				h_["detectorNumCentralAk4JetsOver200GeV"]->Fill(numCentralAk4JetsOver200GeV);
				h_["detectorNumCentralAk4JetsOver20GeV"]->Fill(numCentralAk4JetsOver20GeV);
				h_["detectorNumCentralAk4JetsOver40GeV"]->Fill(numCentralAk4JetsOver40GeV);
				h_["detectorNumCentralAk4JetsOver100GeV"]->Fill(numCentralAk4JetsOver100GeV);
				h_["detectorNumForwardAk4JetsOver200GeV"]->Fill(numForwardAk4JetsOver200GeV);
				h_["detectorNumForwardAk4JetsOver20GeV"]->Fill(numForwardAk4JetsOver20GeV);
				h_["detectorNumForwardAk4JetsOver40GeV"]->Fill(numForwardAk4JetsOver40GeV);
				h_["detectorNumForwardAk4JetsOver100GeV"]->Fill(numForwardAk4JetsOver100GeV);
				h_["detectorLeadingAk4JetPt"]->Fill(leadingAk4JetPt);
				if (jetVec.size() > 0) h_["detectorLeadingAk4JetEta"]->Fill(leadingAk4JetEta);
				h_["detectorSecondaryAk4JetPt"]->Fill(secondaryAk4JetPt);
				if (jetVec.size() > 1) h_["detectorSecondaryAk4JetEta"]->Fill(secondaryAk4JetEta);
				h2_["detectorLeadingAk4JetPt_detectorSecondaryAk4JetPt"]->Fill(secondaryAk4JetPt,leadingAk4JetPt);
				h2_["detectorHT_detectorSecondaryAk4JetPt"]->Fill(secondaryAk4JetPt,ht);
				h2_["detectorHT_detectorLeadingAk4JetPt"]->Fill(leadingAk4JetPt,ht);
				h_["detectorHT"]->Fill(ht);
				h_["detectorMHT"]->Fill(mhtVec[0]);
				h2_["detectorMHT_detectorMET"]->Fill(metDetector->MET,mhtVec[0]);
				h2_["lspMET_detectorMET"]->Fill(metDetector->MET,lspMet[0]);
				h_["detectorMinBiasedDeltaPhi"]->Fill(minBiasedDeltaPhi);
				h_["detectorAlphaT"]->Fill(alphaT);
				double mhtOverMet = mhtVec[0] / metDetector->MET;
				h_["detectorMHToverMET"]->Fill(mhtOverMet);

				if (alphaT <= 0.52 && minBiasedDeltaPhi <= 0.50) h2_["detectorLogicAlphaTMinBiasedDeltaPhi"]->Fill(0.1, 0.1);
				if (alphaT > 0.52 && minBiasedDeltaPhi <= 0.50) h2_["detectorLogicAlphaTMinBiasedDeltaPhi"]->Fill(1.1, 0.1);
				if (alphaT <= 0.52 && minBiasedDeltaPhi > 0.50) h2_["detectorLogicAlphaTMinBiasedDeltaPhi"]->Fill(0.1, 1.1);
				if (alphaT > 0.52 && minBiasedDeltaPhi > 0.50) h2_["detectorLogicAlphaTMinBiasedDeltaPhi"]->Fill(1.1, 1.1);

				if (mhtOverMet <= 1.25 && mhtVec[0] <= 130) h2_["detectorLogicMHToverMETmht"]->Fill(0.1,0.1);
				if (mhtOverMet > 1.25 && mhtVec[0] <= 130) h2_["detectorLogicMHToverMETmht"]->Fill(1.1,0.1);
				if (mhtOverMet <= 1.25 && mhtVec[0] > 130) h2_["detectorLogicMHToverMETmht"]->Fill(0.1,1.1);
				if (mhtOverMet > 1.25 && mhtVec[0] > 130) h2_["detectorLogicMHToverMETmht"]->Fill(1.1,1.1);

				if (mhtVec[0] > 130 && mhtOverMet < 1.25 && minBiasedDeltaPhi > 0.50 && alphaT > 0.52) h2_["detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50andAlphaTmoreThan0p52"]->Fill(1.1, 0.1);
				else h2_["detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50andAlphaTmoreThan0p52"]->Fill(0.1,0.1);

				if (mhtVec[0] > 130 && mhtOverMet < 1.25 && minBiasedDeltaPhi > 0.50) h2_["detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50"]->Fill(1.1, 0.1);
				else h2_["detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50"]->Fill(0.1,0.1);
				// End of Analysis
				// ------------------------------------------------------------------------------------------------------------//
				// ------------------------------------------------------------------------------------------------------------//

			} // closes loop through events for this file
			inFile->Close();
	    } // closes 'if' the file exists

	    // break loop if maximal number of events is reached:
	    // this has to be done twice to stop the file loop as well
	    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;

	} // closes loop through files

	errorRecord.close();
	WriteHistograms(h_, h2_, outputFile_.c_str());

	plottingLabel:
	// sketchy hack to get the higgsMass, squarkMass, ratio and mass splitting from the root file path name
	size_t endNameContainer;
	size_t beginNameContainer;
	for (size_t c = inputFiles_[0].size()-2; c > 0; -- c){
		
		if (inputFiles_[0].substr(c-6,6) == "events") endNameContainer = c; 
		if (inputFiles_[0].substr(c,2) == "mH"){
			beginNameContainer = c;
			break;
		}
	}
	std::string titleName = inputFiles_[0].substr(beginNameContainer,endNameContainer-beginNameContainer);
	// std::cout << titleName << std::endl;
	std::string lowDash = "_";
	for (size_t c = titleName.size()-1; c > 0; --c){
		if (titleName[c] == lowDash[0]){
			endNameContainer = c;
			break;
		}
	}
	titleName = titleName.substr(0,endNameContainer);

	if (std::system(Form("test -e %s",outputFile_.c_str())) == 0) // check the file exists
		PlottingMcSignalStudies createPlots(outputFile_.c_str(), titleName.c_str());
	else{
		std::cout << "the following output root file does not exist to make plots from:" << std::endl;
		std::cout << outputFile_ << std::endl;
	}

	return;
} // closes the function "McSignalStudies"






void CreateHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_)
{
	// Gen Particle Histograms
    h_["numberOfGluinos"] = new TH1F("numberOfGluinos", ";Number of Gluinos;a.u.", 4, 0, 4);
	
	h_["leadingSquarkPt"] = new TH1F("leadingSquarkPt", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["leadingSquarkEta"] = new TH1F("leadingSquarkEta", ";#eta squark;a.u.", 50, -5, 5);
	h_["secondarySquarkPt"] = new TH1F("secondarySquarkPt", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["secondarySquarkEta"] = new TH1F("secondarySquarkEta", ";#eta squark;a.u.", 50, -5, 5);
	h2_["leadingSquarkPt_SecondarySquarkPt"] = new TH2F("leadingSquarkPt_SecondarySquarkPt", ";secondary squark p_{T} (GeV);leading squark p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);
	h2_["leadingSquarkEta_SecondarySquarkEta"] = new TH2F("leadingSquarkEta_SecondarySquarkEta", ";#eta secondary squark;#eta leading squark", 400, -5, 5, 400, -5, 5);
	h2_["leadingSquarkPhi_SecondarySquarkPhi"] = new TH2F("leadingSquarkPhi_SecondarySquarkPhi", ";secondary squark Phi;leading squark Phi", 400, -M_PI, M_PI, 400, -M_PI, M_PI);
	h_["leadingSquarkPt_zeroGluinos"] = new TH1F("leadingSquarkPt_zeroGluinos", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["leadingSquarkPt_oneGluinos"] = new TH1F("leadingSquarkPt_oneGluinos", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["leadingSquarkPt_twoGluinos"] = new TH1F("leadingSquarkPt_twoGluinos", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["leadingSquarkEta_zeroGluinos"] = new TH1F("leadingSquarkEta_zeroGluinos", ";#eta squark;a.u.", 50, -5, 5);
	h_["leadingSquarkEta_oneGluinos"] = new TH1F("leadingSquarkEta_oneGluinos", ";#eta squark;a.u.", 50, -5, 5);
	h_["leadingSquarkEta_twoGluinos"] = new TH1F("leadingSquarkEta_twoGluinos", ";#eta squark;a.u.", 50, -5, 5);
	h_["secondarySquarkPt_zeroGluinos"] = new TH1F("secondarySquarkPt_zeroGluinos", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["secondarySquarkPt_oneGluinos"] = new TH1F("secondarySquarkPt_oneGluinos", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["secondarySquarkPt_twoGluinos"] = new TH1F("secondarySquarkPt_twoGluinos", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["secondarySquarkEta_zeroGluinos"] = new TH1F("secondarySquarkEta_zeroGluinos", ";#eta squark;a.u.", 50, -5, 5);
	h_["secondarySquarkEta_oneGluinos"] = new TH1F("secondarySquarkEta_oneGluinos", ";#eta squark;a.u.", 50, -5, 5);
	h_["secondarySquarkEta_twoGluinos"] = new TH1F("secondarySquarkEta_twoGluinos", ";#eta squark;a.u.", 50, -5, 5);

	h_["leadingQjetPt"] = new TH1F("leadingQjetPt", ";quark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["leadingQjetEta"] = new TH1F("leadingQjetEta", ";#eta quark;a.u.", 50, -5, 5);
	h_["secondaryQjetPt"] = new TH1F("secondaryQjetPt", ";quark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["secondaryQjetEta"] = new TH1F("secondaryQjetEta", ";#eta quark;a.u.", 50, -5, 5);
	h2_["leadingQjetPt_secondaryQjetPt"] = new TH2F("leadingQjetPt_secondaryQjetPt", ";secondary quark p_{T} (GeV);leading quark p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);
	h2_["leadingQjetEta_secondaryQjetEta"] = new TH2F("leadingQjetEta_secondaryQjetEta", ";#eta secondary quark;#eta leading quark", 400, -5, 5, 400, -5, 5);
	h2_["leadingQjetPhi_secondaryQjetPhi"] = new TH2F("leadingQjetPhi_secondaryQjetPhi", ";secondary quark Phi;leading quark Phi", 400, -M_PI, M_PI, 400, -M_PI, M_PI);

	h_["leadingNlspPt"] = new TH1F("leadingNlspPt", ";NLSP p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["leadingNlspEta"] = new TH1F("leadingNlspEta", ";#eta NLSP;a.u.", 50, -5, 5);
	h_["secondaryNlspPt"] = new TH1F("secondaryNlspPt", ";NLSP p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["secondaryNlspEta"] = new TH1F("secondaryNlspEta", ";#eta NLSP;a.u.", 50, -5, 5);

	h_["leadingHiggsPt"] = new TH1F("leadingHiggsPt", ";higgs p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["leadingHiggsEta"] = new TH1F("leadingHiggsEta", "; #eta higgs;a.u.", 50, -5, 5);
	h_["secondaryHiggsPt"] = new TH1F("secondaryHiggsPt", ";higgs p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["secondaryHiggsEta"] = new TH1F("secondaryHiggsEta", ";#eta higgs;a.u.", 50, -5, 5);
	h2_["leadingQjetPt_leadingHiggsPt"] = new TH2F("leadingQjetPt_leadingHiggsPt", ";higgs p_{T} (GeV);quark p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);
	h2_["leadingQjetEta_leadingHiggsEta"] = new TH2F("leadingQjetEta_leadingHiggsEta", ";#eta higgs;#eta quark", 400, -5, 5, 400, -5, 5);
	h2_["leadingQjetPhi_leadingHiggsPhi"] = new TH2F("leadingQjetPhi_leadingHiggsPhi", ";higgs Phi;quark Phi", 400, -M_PI, M_PI, 400, -M_PI, M_PI);
	h2_["secondaryQjetPt_secondaryHiggsPt"] = new TH2F("secondaryQjetPt_secondaryHiggsPt", ";higgs p_{T};quark p_{T}", 400, 0, 2500, 400, 0, 2500);
	h2_["secondaryQjetEta_secondaryHiggsEta"] = new TH2F("secondaryQjetEta_secondaryHiggsEta", ";#eta higgs;#eta quark", 400, -5, 5, 400, -5, 5);
	h2_["secondaryQjetPhi_secondaryHiggsPhi"] = new TH2F("secondaryQjetPhi_secondaryHiggsPhi", ";higgs Phi;quark Phi", 400, -M_PI, M_PI, 400, -M_PI, M_PI);
	h2_["leadingHiggsPt_secondaryHiggsPt"] = new TH2F("leadingHiggsPt_secondaryHiggsPt", ";secondary higgs p_{T} (GeV);leading higgs p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);
	h2_["leadingHiggsEta_secondaryHiggsEta"] = new TH2F("leadingHiggsEta_secondaryHiggsEta", ";#eta secondary higgs;#eta leading higgs", 400, -5, 5, 400, -5, 5);
	h2_["leadingHiggsPhi_secondaryHiggsPhi"] = new TH2F("leadingHiggsPhi_secondaryHiggsPhi", ";secondary higgs Phi;leading higgs Phi", 400, -M_PI, M_PI, 400, -M_PI, M_PI);
	h_["leadingHiggsQjetDphi"] = new TH1F("leadingHiggsQjetDphi", ";higgs Phi - qjet Phi;a.u.", 50, 0, 2*M_PI);
	h_["secondaryHiggsQjetDphi"] = new TH1F("secondaryHiggsQjetDphi", ";higgs Phi - qjet Phi;a.u.", 50, 0, 2*M_PI);
	h_["leadingHiggsSecondaryHiggsDphi"] = new TH1F("leadingHiggsSecondaryHiggsDphi", ";leading higgs Phi - secondary higgs Phi;a.u.", 50, 0, 2*M_PI);
	h_["leadingQjetSecondaryQjetDphi"] = new TH1F("leadingQjetSecondaryQjetDphi", ";leading qjet Phi - secondary qjet Phi;a.u.", 50, 0, 2*M_PI);

	h2_["logicBbDr"] = new TH2F("logicBbDr", ";secondary arm dR_bb > 0.8;leading arm dR_bb > 0.8", 2, 0, 2, 2, 0, 2);
	h2_["logicHiggsPt"] = new TH2F("logicHiggsPt", ";secondary arm higgs p_{T} > 170 (GeV);leading arm higgs p_{T} > 170 (GeV)", 2, 0, 2, 2, 0, 2);
	h2_["logicBbDrHiggsPt"] = new TH2F("logicBbDrHiggsPt", ";higgs p_{T} > 170 (GeV); dR_bb > 0.8", 2, 0, 2, 2, 0, 2);
	h2_["logicNumberHiggsPtNumberBbDr"] = new TH2F("logicNumberHiggsPtNumberBbDr", ";number of higgs with p_{T} > 170 (GeV);number of dR_bb > 0.8", 3, 0, 3, 3, 0, 3);
	
	h_["leadingLspPt"] = new TH1F("leadingLspPt", ";LSP p_{T} (GeV);a.u.", 50, 0, 500);
	h_["leadingLspEta"] = new TH1F("leadingLspEta", ";#eta LSP;a.u.", 50, -5, 5);
	h_["secondaryLspPt"] = new TH1F("secondaryLspPt", ";LSP p_{T} (GeV);a.u.", 50, 0, 500);
	h_["secondaryLspEta"] = new TH1F("secondaryLspEta", ";#eta LSP;a.u.", 50, -5, 5);
	h_["lspMET"] = new TH1F("lspMET", ";LSP E_{T}^{miss} (GeV);a.u.", 50, 0, 500);

	h_["leadingBBbarSeperation"] = new TH1F("leadingBBbarSeperation", ";dR_bb;a.u.", 50, 0, 2.5);
	h_["secondaryBBbarSeperation"] = new TH1F("secondaryBBbarSeperation", ";dR_bb;a.u.", 50, 0, 2.5);
	h2_["leadingBBbarSeperation_secondaryBBbarSeperation"] = new TH2F("leadingBBbarSeperation_secondaryBBbarSeperation", ";secondary dR_bb;leading dR_bb", 400, 0, 2.5, 400, 0, 2.5);	
	h_["leadingBBbarInvmass"] = new TH1F("leadingBBbarInvmass", ";mass bb (GeV);a.u.", 150, 0, 150);
	h_["secondaryBBbarInvmass"] = new TH1F("secondaryBBbarInvmass", ";mass bb (GeV);a.u.", 150, 0, 150);
	h2_["leadingBBbarSeperation_massHiggsOverPt"] = new TH2F("leadingBBbarSeperation_massHiggsOverPt", ";(higgs) mass / p_{T};dR_bb", 400, 0, 1.5, 400, 0, 2.5);
	h2_["secondaryBBbarSeperation_massHiggsOverPt"] = new TH2F("secondaryBBbarSeperation_massHiggsOverPt",";(higgs) mass / p_{T};dR_bb", 400, 0, 1.5, 400, 0, 2.5);
	h2_["leadingBBbarSeperation_massHiggsOverP3"] = new TH2F("leadingBBbarSeperation_massHiggsOverP3", ";(higgs) mass / p3;dR_bb", 400, 0, 1.5, 400, 0, 2.5);
	h2_["secondaryBBbarSeperation_massHiggsOverP3"] = new TH2F("secondaryBBbarSeperation_massHiggsOverP3",";(higgs) mass / p3;dR_bb", 400, 0, 1.5, 400, 0, 2.5);

	// Detector Histograms
	h_["detectorMET"] = new TH1F("detectorMET", ";detector E_{T}^{miss} (GeV);a.u.", 50, 0, 800);
	h_["detectorHT"] = new TH1F("detectorHT", ";detector H_{T} (GeV);a.u.", 50, 0, 7000);
	h2_["lspMET_detectorMET"] = new TH2F("lspMET_detectorMET", ";detector E_{T}^{miss} (GeV);LSP E_{T}^{miss} (GeV)", 400, 0, 800, 400, 0, 800);
	h2_["detectorMHT_detectorMET"] = new TH2F("detectorMHT_detectorMET", ";detector E_{T}^{miss} (GeV);detector H_{T}^{miss} (GeV)", 400, 0, 800, 400, 0, 800);
	h_["detectorMHT"] = new TH1F("detectorMHT", ";detector H_{T}^{miss} (GeV);a.u.", 50, 0, 800);
	h_["detectorMinBiasedDeltaPhi"] = new TH1F("detectorMinBiasedDeltaPhi", ";Min Biased Delta Phi; a.u.", 50, 0, M_PI);
	h_["detectorAlphaT"] = new TH1F("detectorAlphaT", ";Alpha_T; a.u.", 100, 0, 1.0);
	h_["detectorMHToverMET"] = new TH1F("detectorMHToverMET", ";detector H_{T}^{miss} / E_{T}^{miss}; a.u.", 50, 0, 3.0);
	h_["detectorLeadingAk4JetPt"] = new TH1F("detectorLeadingAk4JetPt", ";AK4 Jet p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["detectorSecondaryAk4JetPt"] = new TH1F("detectorSecondaryAk4JetPt", ";AK4 Jet p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["detectorLeadingAk4JetEta"] = new TH1F("detectorLeadingAk4JetEta", ";#eta AK4 Jet;a.u.", 50, -5, 5);
	h_["detectorSecondaryAk4JetEta"] = new TH1F("detectorSecondaryAk4JetEta", ";#eta AK4 Jet;a.u.", 50, -5, 5);
	h2_["detectorLeadingAk4JetPt_detectorSecondaryAk4JetPt"] = new TH2F("detectorLeadingAk4JetPt_detectorSecondaryAk4JetPt", ";secondary AK4 Jet p_{T} (GeV);leading AK4 Jet p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);
	h2_["detectorHT_detectorSecondaryAk4JetPt"] = new TH2F("detectorHT_detectorSecondaryAk4JetPt", ";secondary AK4 Jet p_{T} (GeV);detector H_{T} (GeV)", 400, 0, 2500, 400, 0, 7000);
	h2_["detectorHT_detectorLeadingAk4JetPt"] = new TH2F("detectorHT_detectorLeadingAk4JetPt", ";leading AK4 Jet p_{T} (GeV);detector H_{T} (GeV)", 400, 0, 2500, 400, 0, 7000);
	h_["detectorNumAk4JetsOver200GeV"] = new TH1F("detectorNumAk4JetsOver200GeV", ";Number of AK4 Jets;a.u.", 15, 0, 15);
	h_["detectorNumAk4JetsOver20GeV"] = new TH1F("detectorNumAk4JetsOver20GeV", ";Number of AK4 Jets;a.u.", 15, 0, 15);
	h_["detectorNumAk4JetsOver40GeV"] = new TH1F("detectorNumAk4JetsOver40GeV", ";Number of AK4 Jets;a.u.", 15, 0, 15);
	h_["detectorNumAk4JetsOver100GeV"] = new TH1F("detectorNumAk4JetsOver100GeV", ";Number of AK4 Jets;a.u.", 15, 0, 15);
	h_["detectorNumCentralAk4JetsOver200GeV"] = new TH1F("detectorNumCentralAk4JetsOver200GeV", ";Number of AK4 Jets;a.u.", 15, 0, 15);
	h_["detectorNumCentralAk4JetsOver20GeV"] = new TH1F("detectorNumCentralAk4JetsOver20GeV", ";Number of AK4 Jets;a.u.", 15, 0, 15);
	h_["detectorNumCentralAk4JetsOver40GeV"] = new TH1F("detectorNumCentralAk4JetsOver40GeV", ";Number of AK4 Jets;a.u.", 15, 0, 15);
	h_["detectorNumCentralAk4JetsOver100GeV"] = new TH1F("detectorNumCentralAk4JetsOver100GeV", ";Number of AK4 Jets;a.u.", 15, 0, 15);
	h_["detectorNumForwardAk4JetsOver200GeV"] = new TH1F("detectorNumForwardAk4JetsOver200GeV", ";Number of AK4 Jets;a.u.", 15, 0, 15);
	h_["detectorNumForwardAk4JetsOver20GeV"] = new TH1F("detectorNumForwardAk4JetsOver20GeV", ";Number of AK4 Jets;a.u.", 15, 0, 15);
	h_["detectorNumForwardAk4JetsOver40GeV"] = new TH1F("detectorNumForwardAk4JetsOver40GeV", ";Number of AK4 Jets;a.u.", 15, 0, 15);
	h_["detectorNumForwardAk4JetsOver100GeV"] = new TH1F("detectorNumForwardAk4JetsOver100GeV", ";Number of AK4 Jets;a.u.", 15, 0, 15);
	h2_["detectorLogicAlphaTMinBiasedDeltaPhi"] = new TH2F("detectorLogicAlphaTMinBiasedDeltaPhi", ";Alpha T > 0.52;Biased Delta Phi > 0.50", 2, 0, 2, 2, 0, 2);
	h2_["detectorLogicMHToverMETmht"] = new TH2F("detectorLogicMHToverMETmht", ";H_{T}^{miss} / E_{T}^{miss} > 1.25;H_{T}^{miss} > 130 (GeV)", 2, 0, 2, 2, 0, 2);

	h2_["detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50andAlphaTmoreThan0p52"] = new TH2F("detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50andAlphaTmoreThan0p52", ";pass?;", 2, 0, 2, 1, 0, 1);
	h2_["detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50"] = new TH2F("detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50", ";pass?;", 2, 0, 2, 1, 0, 1);

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





// NOTE THAT THIS GETS THE PYTHIA HARD PROCESS PARTICLES. not the fully hadronization-showered particles
bool indexAllCascadeParticles(std::vector<GenParticle*> particleVec, int ievt, std::string filename, std::ofstream & errorRecord, unsigned int & gluinoCount, std::vector<int> & squarkIndices, std::vector<int> & qjetIndices, std::vector<int> & nlspIndices, std::vector<int> & lspIndices, std::vector<int> & higgsIndices, std::vector<int> & bIndices, std::vector<int> & bbarIndices)
{
// *** SQUARKS ***
// loop through gen particle entries to find the squarks (and count the gluinos)
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( abs(particleVec[iPar]->PID) == 1000021 ) gluinoCount++;
	if ( abs(particleVec[iPar]->PID) == 1000001
	|| abs(particleVec[iPar]->PID) == 2000001
	|| abs(particleVec[iPar]->PID) == 1000002
	|| abs(particleVec[iPar]->PID) == 2000002
	|| abs(particleVec[iPar]->PID) == 1000003
	|| abs(particleVec[iPar]->PID) == 2000003
	|| abs(particleVec[iPar]->PID) == 1000004
	|| abs(particleVec[iPar]->PID) == 2000004) squarkIndices.push_back(iPar);
} // closes loop through gen particle vector
if (squarkIndices.size() != 2){
	errorRecord << "ERROR - Not 2 Squarks" << std::endl;
	errorRecord << "File: " << filename << std::endl;
	errorRecord << "Event: " << ievt << std::endl;
	errorRecord << std::endl;
	return false;
} // closes 'if' we don't have the 2 squarks
// arrange so first element is the highest pt squark
if (particleVec[squarkIndices[1]]->PT > particleVec[squarkIndices[0]]->PT) std::swap(squarkIndices[0],squarkIndices[1]);



// NB: the squark object do not have daughters in root file
// so we have to find cascade particles in a backwards kind of way			

// *** QJET ***
// loop through gen particle entries to find the quarks (squark->QUARK+nlsp)
// nb: currently these are set in madgraph to only be u or d quarks (added s and c)
// leading arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( abs(particleVec[iPar]->PID) == 1 || abs(particleVec[iPar]->PID) == 2 || abs(particleVec[iPar]->PID) == 3 || abs(particleVec[iPar]->PID) == 4){
		if (particleVec[iPar]->M1 == squarkIndices[0]){
			qjetIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
// secondary arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( abs(particleVec[iPar]->PID) == 1 || abs(particleVec[iPar]->PID) == 2 || abs(particleVec[iPar]->PID) == 3 || abs(particleVec[iPar]->PID) == 4){
		if (particleVec[iPar]->M1 == squarkIndices[1]){
			qjetIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
if (qjetIndices.size() != 2){
	errorRecord << "ERROR - Not 2 qjets" << std::endl;
	errorRecord << "File: " << filename << std::endl;
	errorRecord << "Event: " << ievt << std::endl;
	errorRecord << std::endl;
	return false;
} // closes 'if' we don't have the 2 qjets



// *** NLSP ***
// loop through gen particle entries to find the nlsp's (squark->quark+NLSP)
// leading arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( abs(particleVec[iPar]->PID) == 1000023 ){
		if (particleVec[iPar]->M1 == squarkIndices[0]){
			nlspIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
// secondary arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( abs(particleVec[iPar]->PID) == 1000023 ){
		if (particleVec[iPar]->M1 == squarkIndices[1]){
			nlspIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
if (nlspIndices.size() != 2){
	errorRecord << "ERROR - Not 2 nlsp's" << std::endl;
	errorRecord << "File: " << filename << std::endl;
	errorRecord << "Event: " << ievt << std::endl;
	errorRecord << std::endl;
	return false;
} // closes 'if' we don't have the 2 nlsps's



// *** LSP ***
// loop through gen particle entries to find the lsp's (nlsp->LSP+higgs)
// leading arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( abs(particleVec[iPar]->PID) == 1000022 ){
		if (particleVec[iPar]->M1 == nlspIndices[0]){
			lspIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
// secondary arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( abs(particleVec[iPar]->PID) == 1000022 ){
		if (particleVec[iPar]->M1 == nlspIndices[1]){
			lspIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
if (lspIndices.size() != 2){
	errorRecord << "ERROR - Not 2 lsp's" << std::endl;
	errorRecord << "File: " << filename << std::endl;
	errorRecord << "Event: " << ievt << std::endl;
	errorRecord << std::endl;
	return false;
} // closes 'if' we don't have the 2 lsp's



// *** HIGGS ***
// loop through gen particle entries to find the lsp's (nlsp->lsp+HIGGS)
// leading arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( abs(particleVec[iPar]->PID) == 35 ){
		if (particleVec[iPar]->M1 == nlspIndices[0]){
			higgsIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
// secondary arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( abs(particleVec[iPar]->PID) == 35 ){
		if (particleVec[iPar]->M1 == nlspIndices[1]){
			higgsIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
if (higgsIndices.size() != 2){
	errorRecord << "ERROR - Not 2 higgs'" << std::endl;
	errorRecord << "File: " << filename << std::endl;
	errorRecord << "Event: " << ievt << std::endl;
	errorRecord << std::endl;
	return false;
} // closes 'if' we don't have the 2 higgs



//  *** b ***
// loop through gen particle entries to find the b's (higgs->B+bbar)
// leading arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( particleVec[iPar]->PID == 5 ){
		if (particleVec[iPar]->M1 == higgsIndices[0]){
			bIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
// secondary arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( particleVec[iPar]->PID == 5 ){
		if (particleVec[iPar]->M1 == higgsIndices[1]){
			bIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
if (bIndices.size() != 2){
	errorRecord << "ERROR - Not 2 b's" << std::endl;
	errorRecord << "File: " << filename << std::endl;
	errorRecord << "Event: " << ievt << std::endl;
	errorRecord << std::endl;
	return false;
} // closes 'if' we don't have the 2 b's				



//  *** bbar ***
// loop through gen particle entries to find the b's (higgs->b+BBAR)
// leading arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( particleVec[iPar]->PID == -5 ){
		if (particleVec[iPar]->M1 == higgsIndices[0]){
			bbarIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
// secondary arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( particleVec[iPar]->PID == -5 ){
		if (particleVec[iPar]->M1 == higgsIndices[1]){
			bbarIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
if (bIndices.size() != 2){
	errorRecord << "ERROR - Not 2 bbars's" << std::endl;
	errorRecord << "File: " << filename << std::endl;
	errorRecord << "Event: " << ievt << std::endl;
	errorRecord << std::endl;
	return false;
} // closes 'if' we don't have the 2 bbar's	

return true;
}