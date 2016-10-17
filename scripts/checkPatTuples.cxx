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

unsigned int ievt = 0;
for (ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
	if (ievt%1==0) std::cout << "Processing Event: " << ievt << std::endl;
		edm::EventBase const & event = ev;
	}
}