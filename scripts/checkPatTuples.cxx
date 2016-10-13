#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"


#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"



void checkPatTuples(){

gSystem->Load( "libFWCoreFWLite" );
AutoLibraryLoader::enable();

TFile * inFile = TFile::Open("/hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph/CRAB3_testingMCproduction/161006_125105/0000/bTagPatTuple_1.root");
fwlite::Event ev(inFile);

// TH1F * hFatJet_pt = new TH1F("fatJet_pt", "pt (GeV)", 50, 0, 500);
// TH1F * hFatJet_bDisc = new TH1F("fatJet_bDisc", "discrim", 50, -2, 2);
// TH1F * hFatJet_mass = new TH1F("fatJet_mass", "mass (GeV)", 100, 0, 500);

unsigned int ievt = 0;
for (ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
	if (ievt%1==0) std::cout << "Processing Event: " << ievt << std::endl;

edm::EventBase const & event = ev;

		// // Handle to the genParticle collection
		// edm::Handle<std::vector<reco::GenParticle>> genParticles;
		// event.getByLabel(std::string("genParticles"), genParticles);		 

		// // loop genParticle collection
		// for (std::vector<reco::GenParticle>::const_iterator gen=genParticles->begin(); gen!=genParticles->end(); ++gen){
		// 	if (gen->pdgId()==25){ std::cout << "higgs" << std::endl;



  //    int n = gen->numberOfDaughters();
  //    for(size_t j = 0; j < n; ++ j) {
  //    	std::cout << gen->daughter(j)->pdgId() << std::endl;
		// }
	}
}






// edm::Handle<std::vector<pat::Jet>> fatJets;
// event.getByLabel(std::string("selectedPatJetsAK8PFCHS"), fatJets);

// std::cout << fatJets->size() << endl;
// std::vector<pat::Jet>::const_iterator fJ1=fatJets->begin();

		// fJ1->

// const reco::GenJet * gj = fJ1->genJet();

		// edm::Handle<std::vector<reco::GenJet>> GEN;
		// event.getByLabel(std::string("selectedPatJets/genJets"), GEN);

		// std::cout << GEN->size() << std::endl;




	// 	for (std::vector<pat::Jet>::const_iterator fJ1=fatJets->begin(); fJ1!=fatJets->end(); ++fJ1){

	// 		hFatJet_pt->Fill(fJ1->pt());
	// 		hFatJet_bDisc->Fill(fJ1->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));
	// 		hFatJet_mass->Fill(fJ1->mass());

	// 		double abc = fJ1->

	// 	} //closes loop through fatJets entries for the event 
	// } //closes loop through events

	// TFile g( "output_checkPatTuples.root" , "new");
	// hFatJet_pt->Write();
	// hFatJet_bDisc->Write();
	// hFatJet_mass->Write();

// } //closes macro function