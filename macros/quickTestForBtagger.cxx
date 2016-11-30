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
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"

// Headers from this package
#include "Analysis/Analysis_boostedNmssmHiggs/interface/Kinematics.h"


/*
run with 
root -q -b -l $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/macros/quickTestForBtagger.cxx
*/


std::vector<reco::GenParticle> higgsBbGenParticles(edm::Handle<std::vector<reco::GenParticle>>);
bool isThereAFatJetMatch(edm::Handle<std::vector<pat::Jet>>, reco::GenParticle, double, pat::Jet&);
std::string getOutputDirFromOutputFile(std::string);



void quickTestForBtagger() 
{
	gSystem->Load("libFWCoreFWLite.so");
	FWLiteEnabler::enable();
	
	////////////////////
	////////////////////	/////////////////////
	// Set parameters //	/////////////////////

	int maxEvents_ = 10;
	unsigned int outputEvery_ = 1;
	std::vector<std::string> inputFiles_ = {"/users/jt15104/testCMSSW/CMSSW_8_0_3_patch1/src/bTagPatTuple.root"};
	double dRMaxMatch = 0.8; // max dR between higgs boson and fatJet to claim a match

	double doubleBtagWP = 0.3; // loose
	// double doubleBtagWP = 0.6; // medium
	// double doubleBtagWP = 0.8; // tight
	// double doubleBtagWP = 0.9; // very tight

	/////////////////////	/////////////////////
	/////////////////////	/////////////////////
	/////////////////////
		
	// Loop through the input files
	int ievt=0;
	for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
	    
	    // Open input file
	    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
	    if( inFile ){

	    	std::cout << "debug1" << std::endl;		
			fwlite::Event ev(inFile);

			// Loop through the events for this file
			for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
				edm::EventBase const & event = ev;
				std::cout << "debug2" << std::endl;

				// break loop if maximal number of events is reached 
				if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
				
				// event counter
				if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false){
					std::cout << "   File " << iFile+1 << " of " << inputFiles_.size() << ":";
					std::cout << "  processing event: " << ievt << std::endl;
				}

				// Handle to the fatJet collection
				edm::Handle<std::vector<pat::Jet>> fatJets;
				event.getByLabel(std::string("selectedPatJetsAK8PFCHS"), fatJets);

				// Handle to the genParticle collection
				edm::Handle<std::vector<reco::GenParticle>> genParticles;
				event.getByLabel(std::string("genParticles"), genParticles);		 


				// A N A L Y S I S
				// ------------------------------------------------------------------------------------------------------------//
				// Find the higgs to bb particles, store the objects in a vector

				std::vector<reco::GenParticle> higgsBbParticles = higgsBbGenParticles(genParticles); 
				// Loop through the higgsBb particles
				for (size_t iH = 0; iH < higgsBbParticles.size(); ++iH){

					std::cout << "debug: we have a non SM higgs" << std::endl;
					const reco::GenParticle higgsBbParticle = higgsBbParticles[iH];

					// See if there is a fatJet that matches to the higgsBb (closest in dR, must have dR<dRMaxMatch)
					pat::Jet fatJetMatch; // if there is a matching fatJet, this object will contain it
					bool isMatch =isThereAFatJetMatch(fatJets, higgsBbParticle, dRMaxMatch, fatJetMatch);
					std::cout << "debug: we have a matched fat jet" << std::endl;

				} // closes loop through higgsBb Particles
				// ------------------------------------------------------------------------------------------------------------//


			} // closes loop through events for this file
			inFile->Close();
	    } // closes 'if' the file exists

	    // break loop if maximal number of events is reached:
	    // this has to be done twice to stop the file loop as well
	    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;

	} // closes loop through files
	return 0;
} // closes the 'main' function






std::vector<reco::GenParticle> higgsBbGenParticles(edm::Handle<std::vector<reco::GenParticle>> genParticles)
{
	// this is a vector containing the h->bb higgs'
	std::vector<reco::GenParticle> hBbGenParticles;

	for (size_t iGen = 0; iGen < genParticles->size(); ++iGen){

		const reco::GenParticle & genParticle = (*genParticles)[iGen];

			if (genParticle.pdgId()==35){ // particle is a higgs						
				if (genParticle.numberOfDaughters()==2 && abs(genParticle.daughter(0)->pdgId())==5 && abs(genParticle.daughter(1)->pdgId())==5){ // higgs decays to two b-quarks

					hBbGenParticles.push_back(genParticle);

				} // closes 'if' higgs decays to two b-quarks			
			} // closes 'if' genParticle is a higgs
		} // closes loop through genParticles
	return hBbGenParticles;
}	






bool isThereAFatJetMatch(edm::Handle<std::vector<pat::Jet>> fatJetsD, reco::GenParticle higssBbGenParticleD, double dRMaxMatchD, pat::Jet & fatJetMatch)
{
	size_t closestFatJetIndex = 999;
	double dRMin = 9999.99;
							
	for (size_t iFJ = 0; iFJ < fatJetsD->size(); ++iFJ){
		const pat::Jet & fatJet = (*fatJetsD)[iFJ];

		double dR = delR( delPhi( fatJet.phi(),higssBbGenParticleD.phi() ), delEta( fatJet.eta(),higssBbGenParticleD.eta() ) );
		if (dR < dRMin){
			dRMin = dR;
			closestFatJetIndex = iFJ; 
		}
	} // closes loop through fatJets entries for the event 

	// if there is a fatJet that matches, update fatJetMatch with this object and return 'true' 
	if (dRMin < dRMaxMatchD){
		fatJetMatch = (*fatJetsD)[closestFatJetIndex]; 
		return true; 
	} 

	else return false;
} // closes the function 'isThereAFatJetMatch'