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

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"


void checkPatTuples(){

	gSystem->Load( "libFWCoreFWLite" );
	AutoLibraryLoader::enable();

	TFile * inFile = TFile::Open("../python/bTagPatTuple.root");
	fwlite::Event ev(inFile);

	TH1F * hFatJet_pt = new TH1F("fatJet_pt", "pt (GeV)", 50, 0, 500);
	TH1F * hFatJet_bDisc = new TH1F("fatJet_bDisc", "discrim", 50, -2, 2);
	TH1F * hFatJet_mass = new TH1F("fatJet_mass", "mass (GeV)", 100, 0, 500);

	unsigned int ievt = 0;
	for (ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
		if (ievt%10==0) std::cout << "Processing Event: " << ievt << std::endl;
		edm::EventBase const & event = ev;

		edm::Handle<std::vector<pat::Jet>> fatJets;
		event.getByLabel(std::string("selectedPatJetsAK8PFCHS"), fatJets);

		for (std::vector<pat::Jet>::const_iterator fJ1=fatJets->begin(); fJ1!=fatJets->end(); ++fJ1){

			hFatJet_pt->Fill(fJ1->pt());
			hFatJet_bDisc->Fill(fJ1->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));
			hFatJet_mass->Fill(fJ1->mass());
		} //closes loop through fatJets entries for the event 
	} //closes loop through events

	TFile g( "output_checkPatTuples.root" , "new");
	hFatJet_pt->Write();
	hFatJet_bDisc->Write();
	hFatJet_mass->Write();

} //closes macro function