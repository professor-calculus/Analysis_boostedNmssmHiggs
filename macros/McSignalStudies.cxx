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
// $ root -l -b -q McSignalStudies.cxx

void CreateHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&);
void WriteHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::string);
bool indexAllCascadeParticles(std::vector<GenParticle*>,int,std::string,std::ofstream&,unsigned int&,std::vector<int>&,std::vector<int>&,std::vector<int>&,std::vector<int>&,std::vector<int>&,std::vector<int>&,std::vector<int>&);

void McSignalStudies() 
{
	// Running Options
	int maxEvents_ = -1; // -1 for all events
	unsigned int outputEvery_ = 1000;
	std::vector<std::string> inputFiles_ = {"/storage/jt15104/madGraphProjects/testing/mH125p0_mSusy1000p0_ratio0p96_splitting2p0_5000events/dirA/dirB/dirC/tag_1_delphes_events.root"};
	std::string outputFile_ = "ABC_v3/output.root";
	bool justDoPlotting_ = false;
	// bool justDoPlotting_ = true;

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

				// load the ht info
				if (branchHT->GetEntries() != 1){
					std::cout << "We do not have ONE ht quantity: Exiting..." << std::endl;
					return;
				}
				ScalarHT* htDetector = (ScalarHT*) branchHT->At(0);

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
				h2_["leadingNlspPt_secondaryNlspPt"]->Fill(particleVec[nlspIndices[1]]->PT, particleVec[nlspIndices[0]]->PT);
				h2_["leadingNlspEta_secondaryNlspEta"]->Fill(particleVec[nlspIndices[1]]->Eta, particleVec[nlspIndices[0]]->Eta);
				h2_["leadingNlspPhi_secondaryNlspPhi"]->Fill(particleVec[nlspIndices[1]]->Phi, particleVec[nlspIndices[0]]->Phi);

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

				h_["leadingLspPt"]->Fill(particleVec[lspIndices[0]]->PT);
				h_["leadingLspEta"]->Fill(particleVec[lspIndices[0]]->Eta);
				h_["secondaryLspPt"]->Fill(particleVec[lspIndices[1]]->PT);
				h_["secondaryLspEta"]->Fill(particleVec[lspIndices[1]]->Eta);
				//function: adds two pt vectors together in the (2d) transverse plane
				//output: first element is magnitude, second element is phi
				std::vector<double> lspMet = addTwoPtVectors(particleVec[lspIndices[0]]->PT, particleVec[lspIndices[1]]->PT, particleVec[lspIndices[0]]->Phi, particleVec[lspIndices[1]]->Phi);
				h_["lspMET"]->Fill(lspMet[0]);
				h_["detectorMET"]->Fill(metDetector->MET);
				h_["detectorHT"]->Fill(htDetector->HT);
				// get the leading and second leading detector jets
				double leadingDetJetPt = 0.0;
				double secondaryDetJetPt = 0.0;
				for (size_t iJet=0; iJet<jetVec.size(); ++iJet){
					if (jetVec[iJet]->PT > leadingDetJetPt){
						secondaryDetJetPt = leadingDetJetPt;
						leadingDetJetPt = jetVec[iJet]->PT;
					}
					else if (jetVec[iJet]->PT > secondaryDetJetPt) secondaryDetJetPt = jetVec[iJet]->PT;
				} // closes loop through gen particle vector
				h_["detectorLeadingJet"]->Fill(leadingDetJetPt);
				h_["detectorSecondaryJet"]->Fill(secondaryDetJetPt);

				double leadingBBbarInvmass = invMass_v1(particleVec[bIndices[0]]->E, particleVec[bbarIndices[0]]->E, particleVec[bIndices[0]]->Px, particleVec[bbarIndices[0]]->Px, particleVec[bIndices[0]]->Py, particleVec[bbarIndices[0]]->Py, particleVec[bIndices[0]]->Pz, particleVec[bbarIndices[0]]->Pz, particleVec[bIndices[0]]->Mass, particleVec[bbarIndices[0]]->Mass);
				h_["leadingBBbarInvmass"]->Fill(leadingBBbarInvmass);
				double leadingBBbarSeperation = delR( delPhi(particleVec[bIndices[0]]->Phi, particleVec[bbarIndices[0]]->Phi), delEta(particleVec[bIndices[0]]->Eta, particleVec[bbarIndices[0]]->Eta) );
				h_["leadingBBbarSeperation"]->Fill(leadingBBbarSeperation);
				double secondaryBBbarInvmass = invMass_v1(particleVec[bIndices[1]]->E, particleVec[bbarIndices[1]]->E, particleVec[bIndices[1]]->Px, particleVec[bbarIndices[1]]->Px, particleVec[bIndices[1]]->Py, particleVec[bbarIndices[1]]->Py, particleVec[bIndices[1]]->Pz, particleVec[bbarIndices[1]]->Pz, particleVec[bIndices[1]]->Mass, particleVec[bbarIndices[1]]->Mass);
				h_["secondaryBBbarInvmass"]->Fill(secondaryBBbarInvmass);
				double secondaryBBbarSeperation = delR( delPhi(particleVec[bIndices[1]]->Phi, particleVec[bbarIndices[1]]->Phi), delEta(particleVec[bIndices[1]]->Eta, particleVec[bbarIndices[1]]->Eta) );
				h_["secondaryBBbarSeperation"]->Fill(secondaryBBbarSeperation);
				h2_["leadingBBbarSeperation_massHiggsOverPt"]->Fill(particleVec[higgsIndices[0]]->Mass / (particleVec[higgsIndices[0]]->PT) , leadingBBbarSeperation);
				h2_["secondaryBBbarSeperation_massHiggsOverPt"]->Fill(particleVec[higgsIndices[1]]->Mass / (particleVec[higgsIndices[1]]->PT), secondaryBBbarSeperation);

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
	std::cout << titleName << std::endl;
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

return 0;
} // closes the 'main' function






void CreateHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_)
{
	// set the binning for histograms
    std::vector<double> ptBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 100.0; binLowerEdge+= 100.0) ptBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  100.0; binLowerEdge< 200.0; binLowerEdge+= 50.0) ptBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  200.0; binLowerEdge< 350.0; binLowerEdge+= 25.0) ptBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  350.0; binLowerEdge< 500.0; binLowerEdge+= 50.0) ptBinning.push_back(binLowerEdge);    	
    for(double binLowerEdge=  500.0; binLowerEdge< 600.0; binLowerEdge+= 100.0) ptBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  600.0; binLowerEdge< 800.1; binLowerEdge+= 200.0) ptBinning.push_back(binLowerEdge);

	// create the histograms
    h_["numberOfGluinos"] = new TH1F("numberOfGluinos", ";Number of Gluinos;a.u.", 4, 0, 4);
	
	h_["leadingSquarkPt"] = new TH1F("leadingSquarkPt", ";squark p_{T} (GeV);a.u.", 50, 0, 1000);
	h_["leadingSquarkEta"] = new TH1F("leadingSquarkEta", ";#eta squark;a.u.", 50, -5, 5);
	h_["secondarySquarkPt"] = new TH1F("secondarySquarkPt", ";squark p_{T} (GeV);a.u.", 50, 0, 1000);
	h_["secondarySquarkEta"] = new TH1F("secondarySquarkEta", ";#eta squark;a.u.", 50, -5, 5);
	h2_["leadingSquarkPt_SecondarySquarkPt"] = new TH2F("leadingSquarkPt_SecondarySquarkPt", ";secondary squark p_{T} (GeV);leading squark p_{T} (GeV)", 100, 0, 1000, 100, 0, 1000);
	h2_["leadingSquarkEta_SecondarySquarkEta"] = new TH2F("leadingSquarkEta_SecondarySquarkEta", ";#eta secondary quark;#eta leading quark", 50, -5, 5, 100, 05, 5);
	h2_["leadingSquarkPhi_SecondarySquarkPhi"] = new TH2F("leadingSquarkPhi_SecondarySquarkPhi", ";secondary squark Phi;leading squark Phi", 100, -M_PI, M_PI, 100, -M_PI, M_PI);
	h_["leadingSquarkPt_zeroGluinos"] = new TH1F("leadingSquarkPt_zeroGluinos", ";squark p_{T} (GeV);a.u.", 25, 0, 1500);
	h_["leadingSquarkPt_oneGluinos"] = new TH1F("leadingSquarkPt_oneGluinos", ";squark p_{T} (GeV);a.u.", 25, 0, 1500);
	h_["leadingSquarkPt_twoGluinos"] = new TH1F("leadingSquarkPt_twoGluinos", ";squark p_{T} (GeV);a.u.", 25, 0, 1500);
	h_["leadingSquarkEta_zeroGluinos"] = new TH1F("leadingSquarkEta_zeroGluinos", ";#eta squark;a.u.", 25, -5, 5);
	h_["leadingSquarkEta_oneGluinos"] = new TH1F("leadingSquarkEta_oneGluinos", ";#eta squark;a.u.", 25, -5, 5);
	h_["leadingSquarkEta_twoGluinos"] = new TH1F("leadingSquarkEta_twoGluinos", ";#eta squark;a.u.", 25, -5, 5);
	h_["secondarySquarkPt_zeroGluinos"] = new TH1F("secondarySquarkPt_zeroGluinos", ";squark p_{T} (GeV);a.u.", 25, 0, 1500);
	h_["secondarySquarkPt_oneGluinos"] = new TH1F("secondarySquarkPt_oneGluinos", ";squark p_{T} (GeV);a.u.", 25, 0, 1500);
	h_["secondarySquarkPt_twoGluinos"] = new TH1F("secondarySquarkPt_twoGluinos", ";squark p_{T} (GeV);a.u.", 25, 0, 1500);
	h_["secondarySquarkEta_zeroGluinos"] = new TH1F("secondarySquarkEta_zeroGluinos", ";#eta squark;a.u.", 25, -5, 5);
	h_["secondarySquarkEta_oneGluinos"] = new TH1F("secondarySquarkEta_oneGluinos", ";#eta squark;a.u.", 25, -5, 5);
	h_["secondarySquarkEta_twoGluinos"] = new TH1F("secondarySquarkEta_twoGluinos", ";#eta squark;a.u.", 25, -5, 5);

	h_["leadingQjetPt"] = new TH1F("leadingQjetPt", ";quark p_{T} (GeV);a.u.", 50, 0, 1500);
	h_["leadingQjetEta"] = new TH1F("leadingQjetEta", ";#eta quark;a.u.", 50, -5, 5);
	h_["secondaryQjetPt"] = new TH1F("secondaryQjetPt", ";quark p_{T} (GeV);a.u.", 50, 0, 1500);
	h_["secondaryQjetEta"] = new TH1F("secondaryQjetEta", ";#eta quark;a.u.", 50, -5, 5);
	h2_["leadingQjetPt_secondaryQjetPt"] = new TH2F("leadingQjetPt_secondaryQjetPt", ";secondary quark p_{T} (GeV);leading quark p_{T} (GeV)", 100, 0, 1000, 100, 0, 1000);
	h2_["leadingQjetEta_secondaryQjetEta"] = new TH2F("leadingQjetEta_secondaryQjetEta", ";#eta secondary quark;#eta leading quark", 100, -5, 5, 100, -5, 5);
	h2_["leadingQjetPhi_secondaryQjetPhi"] = new TH2F("leadingQjetPhi_secondaryQjetPhi", ";secondary quark Phi;leading quark Phi", 100, -M_PI, M_PI, 100, -M_PI, M_PI);

	h_["leadingNlspPt"] = new TH1F("leadingNlspPt", ";NLSP p_{T} (GeV);a.u.", 50, 0, 2000);
	h_["leadingNlspEta"] = new TH1F("leadingNlspEta", ";#eta NLSP;a.u.", 50, -5, 5);
	h_["secondaryNlspPt"] = new TH1F("secondaryNlspPt", ";NLSP p_{T} (GeV);a.u.", 50, 0, 2000);
	h_["secondaryNlspEta"] = new TH1F("secondaryNlspEta", ";#eta NLSP;a.u.", 50, -5, 5);
	h2_["leadingNlspPt_secondaryNlspPt"] = new TH2F("leadingNlspPt_secondaryNlspPt", ";secondary NLSP p_{T} (GeV);leading NLSP p_{T} (GeV)", 100, 0, 1000, 100, 0, 1000);
	h2_["leadingNlspEta_secondaryNlspEta"] = new TH2F("leadingNlspEta_secondaryNlspEta", ";#eta secondary NLSP;#eta leading NLSP", 50, -5, 5, 100, 05, 5);
	h2_["leadingNlspPhi_secondaryNlspPhi"] = new TH2F("leadingNlspPhi_secondaryNlspPhi", ";secondary NLSP Phi;leading NLSP Phi", 100, -M_PI, M_PI, 100, -M_PI, M_PI);

	h_["leadingHiggsPt"] = new TH1F("leadingHiggsPt", ";higgs p_{T} (GeV);a.u.", 50, 0, 2000);
	h_["leadingHiggsEta"] = new TH1F("leadingHiggsEta", "; #eta higgs;a.u.", 50, -5, 5);
	h_["secondaryHiggsPt"] = new TH1F("secondaryHiggsPt", ";higgs p_{T} (GeV);a.u.", 50, 0, 2000);
	h_["secondaryHiggsEta"] = new TH1F("secondaryHiggsEta", ";#eta higgs;a.u.", 50, -5, 5);
	h2_["leadingQjetPt_leadingHiggsPt"] = new TH2F("leadingQjetPt_leadingHiggsPt", ";higgs p_{T} (GeV);quark p_{T} (GeV)", 100, 0, 1000, 100, 0, 1000);
	h2_["leadingQjetEta_leadingHiggsEta"] = new TH2F("leadingQjetEta_leadingHiggsEta", ";#eta higgs;#eta quark", 100, -5, 5, 100, -5, 5);
	h2_["leadingQjetPhi_leadingHiggsPhi"] = new TH2F("leadingQjetPhi_leadingHiggsPhi", ";higgs Phi;quark Phi", 100, -M_PI, M_PI, 100, -M_PI, M_PI);
	h2_["secondaryQjetPt_secondaryHiggsPt"] = new TH2F("secondaryQjetPt_secondaryHiggsPt", ";higgs p_{T};quark p_{T}", 100, 0, 1000, 100, 0, 1000);
	h2_["secondaryQjetEta_secondaryHiggsEta"] = new TH2F("secondaryQjetEta_secondaryHiggsEta", ";#eta higgs;#eta quark", 100, -5, 5, 100, -5, 5);
	h2_["secondaryQjetPhi_secondaryHiggsPhi"] = new TH2F("secondaryQjetPhi_secondaryHiggsPhi", ";higgs Phi;quark Phi", 100, -M_PI, M_PI, 100, -M_PI, M_PI);

	h_["leadingLspPt"] = new TH1F("leadingLspPt", ";LSP p_{T} (GeV);a.u.", 50, 0, 300);
	h_["leadingLspEta"] = new TH1F("leadingLspEta", ";#eta LSP;a.u.", 50, -5, 5);
	h_["secondaryLspPt"] = new TH1F("secondaryLspPt", ";LSP p_{T} (GeV);a.u.", 50, 0, 300);
	h_["secondaryLspEta"] = new TH1F("secondaryLspEta", ";#eta LSP;a.u.", 50, -5, 5);
	h_["lspMET"] = new TH1F("lspMET", ";LSP E_{T}^{miss};a.u.", 50, 0, 300);
	h_["detectorMET"] = new TH1F("detectorMET", ";detector E_{T}^{miss} (GeV);a.u.", 10, 0, 600);
	h_["detectorHT"] = new TH1F("detectorHT", ";detector HT (GeV);a.u.", 50, 0, 5000);
	h_["detectorLeadingJet"] = new TH1F("detectorLeadingJet", ";detector Leading Jet p_{T} (GeV);a.u.", 50, 0, 2000);
	h_["detectorSecondaryJet"] = new TH1F("detectorSecondaryJet", ";detector Secondary Jet p_{T} (GeV);a.u.", 50, 0, 2000);

	h_["leadingBBbarSeperation"] = new TH1F("leadingBBbarSeperation", ";dR_bb;a.u.", 50, 0, 2.5);
	h_["leadingBBbarInvmass"] = new TH1F("leadingBBbarInvmass", ";mass bb (GeV);a.u.", 100, 60, 140);
	h_["secondaryBBbarSeperation"] = new TH1F("secondaryBBbarSeperation", ";dR_bb;a.u.", 50, 0, 2.5);
	h_["secondaryBBbarInvmass"] = new TH1F("secondaryBBbarInvmass", ";mass bb (GeV);a.u.", 100, 60, 140);
	h2_["leadingBBbarSeperation_massHiggsOverPt"] = new TH2F("leadingBBbarSeperation_massHiggsOverPt", ";(higgs) mass / p_{T};dR_bb", 100, 0, 1.5, 100, 0, 2.5);
	h2_["secondaryBBbarSeperation_massHiggsOverPt"] = new TH2F("secondaryBBbarSeperation_massHiggsOverPt",";(higgs) mass / p_{T};dR_bb", 100, 0, 1.5, 100, 0, 2.5);

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
// nb: currently these are set in madgraph to only be u or d quarks
// leading arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( abs(particleVec[iPar]->PID) == 1 || abs(particleVec[iPar]->PID) == 2 ){
		if (particleVec[iPar]->M1 == squarkIndices[0]){
			qjetIndices.push_back(iPar);  
			break;
		}
	}
} // closes loop through gen particle vector
// secondary arm quark
for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
	if ( abs(particleVec[iPar]->PID) == 1 || abs(particleVec[iPar]->PID) == 2 ){
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