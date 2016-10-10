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
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

// Headers from this package

// command to make work!
// ~/CMSSW_8_0_20/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies/DoubleBTaggerEfficiencyStudies inputFiles=/users/jt15104/CMSSW_8_0_20/src/Analysis/Analysis_boostedNmssmHiggs/python/bTagPatTuple_testingv1.root 




void CreateHistograms(const std::vector<std::string> doubleBtagWPnameD,
	                  // std::map<std::string, TH1F*> & h_, std::map<std::string, TH2F*> & h2_) // both of these lines seem to work...
	                  std::map<std::string, TH1F*> h_, std::map<std::string, TH2F*> h2_) // both of these lines seem to work...
{

	const std::vector<double> etaBinning = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};

    std::vector<double> ptBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 800.1; binLowerEdge+= 20.0) ptBinning.push_back(binLowerEdge);

    std::vector<double> ptScatXBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 800.1; binLowerEdge+= 20.0) ptScatXBinning.push_back(binLowerEdge);

    std::vector<double> ptScatYBinning = ptScatXBinning;
    // for(double binLowerEdge=  0.0; binLowerEdge< 800.1; binLowerEdge+= 20.0) ptScatYBinning.push_back(binLowerEdge);	

    std::vector<double> massBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 200.1; binLowerEdge+= 10.0) massBinning.push_back(binLowerEdge);


	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPnameD.size(); ++iWP){

		h2_[Form("ptScatter_%s", doubleBtagWPnameD[iWP].c_str())] = new TH2F(
			Form("ptScatter_%s", doubleBtagWPnameD[iWP].c_str()), ";higgs p_{T} (GeV); fatJet p_{T} (GeV)",
			ptScatXBinning.size()-1, &(ptScatXBinning)[0], ptScatYBinning.size()-1, &(ptScatYBinning)[0]);

		h_[Form("fatJetMass_%s", doubleBtagWPnameD[iWP].c_str())] = new TH1F(
		   Form("fatJetMass_%s", doubleBtagWPnameD[iWP].c_str()), ";m_{Jet} (GeV); a.u.", massBinning.size()-1, &(massBinning)[0]);

   		for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinning.size()-1; ++iEtaBin){

   			h_[Form("eff_%sDoubleBTagWP_eta%f-%f", doubleBtagWPnameD[iWP].c_str(), etaBinning[iEtaBin], etaBinning[iEtaBin+1] )] = new TH1F(
   			   Form("eff_%sDoubleBTagWP_eta%f-%f", doubleBtagWPnameD[iWP].c_str(), etaBinning[iEtaBin], etaBinning[iEtaBin+1] ),
   			   ";p_{T} (GeV); efficiency", ptBinning.size()-1, &(ptBinning)[0]); 

   		} // closes loop through etaBins
	} // closes loop through Btag WPs 

} //closes CreateHistograms() function


void WriteHistograms(const std::string outputFileD, std::map<std::string, TH1F*> h_, std::map<std::string, TH2F*> h2_)
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





int main(int argc, char* argv[]) 
{
  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views

  // ----------------------------------------------------------------------
  // First Part: 
  //
  //  * enable the AutoLibraryLoader 
  //  * book the histograms of interest 
  //  * open the input file
  // ----------------------------------------------------------------------

	// load framework libraries...fix this up to current version...
gSystem->Load( "libFWCoreFWLite" );
AutoLibraryLoader::enable();
	
	// Set parameters
	const std::vector<double> doubleBtagWP = {0.3, 0.6, 0.8, 0.9};
	const std::vector<std::string> doubleBtagWPname = {"loose", "medium", "tight", "extraTight"};
	const double dRMax = 0.8; // dR between higgs boson and bTaggedFatJet

	// Create histograms, they are accessed by eg: h_["fatJetMass_loose"]->Fill(125.0);
	std::map<std::string, TH1F*> h_;
	std::map<std::string, TH2F*> h2_;
	CreateHistograms(doubleBtagWPname, h_, h2_);


	// initialize command line parser
	optutl::CommandLineParser parser ("Analyze DoubleBTagger Efficiencies");

	// set defaults
	parser.integerValue ("maxEvents"  ) = -1;
	parser.integerValue ("outputEvery") =   10;
	parser.stringValue  ("outputFile" ) = "output_DoubleBTaggerEfficiencyStudies.root";
	// parser.stringValue  ("inputFiles" ) = {"../python/bTagPatTuple_testingv1.root"};
// I kind of want an output directory...and the name of this controls what is going on

	// parse arguments
	parser.parseArguments (argc, argv);
	int maxEvents_ = parser.integerValue("maxEvents");
	unsigned int outputEvery_ = parser.integerValue("outputEvery");
	std::string outputFile_ = parser.stringValue("outputFile");
	std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");



		// Loop through the input files
		int ievt=0;  
		for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
	    
	    // Open input files
	    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
	    if( inFile ){
		
			fwlite::Event ev(inFile);
			// Loop through the events for this file
			for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
				edm::EventBase const & event = ev;

				// break loop if maximal number of events is reached 
				if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
				// simple event counter
				if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false) 
				std::cout << "  processing event: " << ievt << std::endl;


				// Handle to the fatJet collection
				edm::Handle<std::vector<pat::Jet>> fatJets;
				event.getByLabel(std::string("selectedPatJetsAK8PFCHS"), fatJets);

				// Handle to the genParticle collection
				edm::Handle<std::vector<reco::GenParticle>> genParticles;
				event.getByLabel(std::string("genParticles"), genParticles);		 


				// Loop through genParticle collection to find the h->bb
				// std::vector<size_t> hBbGenIndices;
				for (size_t iGen = 0; iGen < genParticles->size(); ++iGen){
				const reco::GenParticle & genParticle = (*genParticles)[iGen];

					if (genParticle.pdgId()==25){ // particle is a higgs
						
						if (genParticle.numberOfDaughters()==2 && abs(genParticle.daughter(0)->pdgId())==5 && abs(genParticle.daughter(1)->pdgId())==5){ // higgs decays to two b-quarks
							
							// hBbGenIndices.push_back(iGen); // this vector indexes the higgs bosons of interest


							// Loop through fatJet collection
							// find closest matching jet to the higgs, check it is a match, see if it passes btag WP cuts
							

							size_t closestFatJetIndex;
							for (size_t iFJ = 0; iFJ < fatJets->size(); ++iFJ){
									const pat::Jet & fatJet = (*fatJets)[iFJ];




							} // closes loop through fatJets entries for the event 


							//now do histogram filling...


						} // closes 'if' higgs decays to two b-quarks			
					} // closes 'if' genParticle is a higgs
				} // closes loop through genParticles










			} // closes loop through events for this file
			inFile->Close();
	    } // closes 'if' the file exists

	    // break loop if maximal number of events is reached:
	    // this has to be done twice to stop the file loop as well
	    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;

	} // closes loop through files

WriteHistograms(outputFile_.c_str(), h_, h2_);

return 0;
} // closes the 'main' function