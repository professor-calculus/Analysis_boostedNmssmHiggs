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




// could put all these vectors in the function
void CreateHistograms(const std::vector<std::string> doubleBtagWPname, const std::vector<double> etaBinning, const std::vector<double> ptBinning,
	                  const std::vector<double> ptScatXBinning, const std::vector<double> ptScatYBinning, const std::vector<double> massBinning,
	                  // std::map<std::string, TH1F*> & h_, std::map<std::string, TH2F*> & h2_)
	                  std::map<std::string, TH1F*> h_, std::map<std::string, TH2F*> h2_) // both of these lines seem to work...

{

	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){

		h2_[Form("ptScatter_%s", doubleBtagWPname[iWP].c_str())] = new TH2F(
			Form("ptScatter_%s", doubleBtagWPname[iWP].c_str()), ";higgs p_{T} (GeV); fatJet p_{T} (GeV)",
			ptScatXBinning.size()-1, &(ptScatXBinning)[0], ptScatYBinning.size()-1, &(ptScatYBinning)[0]);

		h_[Form("fatJetMass_%s", doubleBtagWPname[iWP].c_str())] = new TH1F(
		   Form("fatJetMass_%s", doubleBtagWPname[iWP].c_str()), ";m_{Jet} (GeV); a.u.", massBinning.size()-1, &(massBinning)[0]);

   		for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinning.size()-1; ++iEtaBin){

   			h_[Form("eff_%sDoubleBTagWP_eta%f-%f", doubleBtagWPname[iWP].c_str(), etaBinning[iEtaBin], etaBinning[iEtaBin+1] )] = new TH1F(
   			   Form("eff_%sDoubleBTagWP_eta%f-%f", doubleBtagWPname[iWP].c_str(), etaBinning[iEtaBin], etaBinning[iEtaBin+1] ),
   			   ";p_{T} (GeV); efficiency", ptBinning.size()-1, &(ptBinning)[0]); 

   		} // closes loop through etaBins
   } // closes loop through Btag WPs 



} //closes CreateHistograms() function


// void WriteHistograms(const std::string & filenameD)
// {
//    TFile * outFile = new TFile(filenameD.c_str(),"RECREATE");
//    for ( auto & h : h_ ){
// 		h.second->Write();
// 	}
// 	for ( auto & h : h2_ ){
// 		h.second->Write();
// 	}
//    outFile -> Close();
//    delete outFile;
// }





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

	// load framework libraries
gSystem->Load( "libFWCoreFWLite" );
AutoLibraryLoader::enable();
	// set parameters
	const std::vector<double> doubleBtagWP = {0.3, 0.6, 0.8, 0.9};
	const std::vector<std::string> doubleBtagWPname = {"loose", "medium", "tight", "extraTight"};

	const std::vector<double> etaBinning = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};

    std::vector<double> ptBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 800.0; binLowerEdge+= 20.0) ptBinning.push_back(binLowerEdge);
    // for(double binLowerEdge= 40.0; binLowerEdge< 70.0; binLowerEdge+= 2.0) ptBinning.push_back(binLowerEdge);
    // for(double binLowerEdge= 70.0; binLowerEdge<100.0; binLowerEdge+= 2.0) ptBinning.push_back(binLowerEdge);
    // for(double binLowerEdge=100.0; binLowerEdge<160.0; binLowerEdge+= 5.0) ptBinning.push_back(binLowerEdge);
    // for(double binLowerEdge=160.0; binLowerEdge<200.1; binLowerEdge+=10.0) ptBinning.push_back(binLowerEdge);

    std::vector<double> ptScatXBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 800.0; binLowerEdge+= 20.0) ptScatXBinning.push_back(binLowerEdge);
    // for(double binLowerEdge= 40.0; binLowerEdge< 70.0; binLowerEdge+= 2.0) ptScatXBinning.push_back(binLowerEdge);
    // for(double binLowerEdge= 70.0; binLowerEdge<100.0; binLowerEdge+= 2.0) ptScatXBinning.push_back(binLowerEdge);
    // for(double binLowerEdge=100.0; binLowerEdge<160.0; binLowerEdge+= 5.0) ptScatXBinning.push_back(binLowerEdge);
    // for(double binLowerEdge=160.0; binLowerEdge<200.1; binLowerEdge+=10.0) ptScatXBinning.push_back(binLowerEdge);

    std::vector<double> ptScatYBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 800.0; binLowerEdge+= 20.0) ptScatYBinning.push_back(binLowerEdge);
    // for(double binLowerEdge= 40.0; binLowerEdge< 70.0; binLowerEdge+= 2.0) ptScatYBinning.push_back(binLowerEdge);
    // for(double binLowerEdge= 70.0; binLowerEdge<100.0; binLowerEdge+= 2.0) ptScatYBinning.push_back(binLowerEdge);
    // for(double binLowerEdge=100.0; binLowerEdge<160.0; binLowerEdge+= 5.0) ptScatYBinning.push_back(binLowerEdge);
    // for(double binLowerEdge=160.0; binLowerEdge<200.1; binLowerEdge+=10.0) ptScatYBinning.push_back(binLowerEdge);	

    std::vector<double> massBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 200.0; binLowerEdge+= 10.0) massBinning.push_back(binLowerEdge);
    // for(double binLowerEdge= 40.0; binLowerEdge< 70.0; binLowerEdge+= 2.0) massBinning.push_back(binLowerEdge);
    // for(double binLowerEdge= 70.0; binLowerEdge<100.0; binLowerEdge+= 2.0) massBinning.push_back(binLowerEdge);
    // for(double binLowerEdge=100.0; binLowerEdge<160.0; binLowerEdge+= 5.0) massBinning.push_back(binLowerEdge);
    // for(double binLowerEdge=160.0; binLowerEdge<200.1; binLowerEdge+=10.0) massBinning.push_back(binLowerEdge);

	// Histograms
	std::map<std::string, TH1F*> h_;
	std::map<std::string, TH2F*> h2_;
	CreateHistograms(doubleBtagWPname, etaBinning, ptBinning, ptScatXBinning, ptScatYBinning, massBinning, h_, h2_);

	h_["fatJetMass_loose"]->Fill(125.0);

	// for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){

	// 	h2_[Form("ptScatter_%s", doubleBtagWPname[iWP].c_str())] = new TH2F(
	// 		Form("ptScatter_%s", doubleBtagWPname[iWP].c_str()), ";higgs p_{T} (GeV); fatJet p_{T} (GeV)",
	// 		ptScatXBinning.size()-1, &(ptScatXBinning)[0], ptScatYBinning.size()-1, &(ptScatYBinning)[0]);

	// 	h_[Form("fatJetMass_%s", doubleBtagWPname[iWP].c_str())] = new TH1F(
	// 	   Form("fatJetMass_%s", doubleBtagWPname[iWP].c_str()), ";m_{Jet} (GeV); a.u.", massBinning.size()-1, &(massBinning)[0]);

 //   		for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinning.size()-1; ++iEtaBin){

 //   			h_[Form("eff_%sDoubleBTagWP_eta%f-%f", doubleBtagWPname[iWP].c_str(), etaBinning[iEtaBin], etaBinning[iEtaBin+1] )] = new TH1F(
 //   			   Form("eff_%sDoubleBTagWP_eta%f-%f", doubleBtagWPname[iWP].c_str(), etaBinning[iEtaBin], etaBinning[iEtaBin+1] ),
 //   			   ";p_{T} (GeV); efficiency", ptBinning.size()-1, &(ptBinning)[0]); 

 //   		} // closes loop through etaBins
 //   } // closes loop through Btag WPs 


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



	// loop the events
	int ievt=0;  
	for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
    // open input file (can be located on castor)
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if( inFile ){
      // ----------------------------------------------------------------------
      // Second Part: 
      //
      //  * loop the events in the input file 
      //  * receive the collections of interest via fwlite::Handle
      //  * fill the histograms
      //  * after the loop close the input file
      // ----------------------------------------------------------------------      
		fwlite::Event ev(inFile);
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

		// loop genParticle collection
		for (std::vector<reco::GenParticle>::const_iterator gen=genParticles->begin(); gen!=genParticles->end(); ++gen){
			if (gen->pdgId()==25){
				std::cout << "higgs" << std::endl;

				for(size_t j = 0; j < gen->numberOfDaughters(); ++ j){
					std::cout << gen->daughter(j)->pdgId() << std::endl;
				}
			}
		}


	// // loop jet collection and fill histograms
	// 	for (std::vector<pat::Jet>::const_iterator fJ1=fatJets->begin(); fJ1!=fatJets->end(); ++fJ1){
	// 		hFatJet_pt->Fill(fJ1->pt());
	// 		hFatJet_bDisc->Fill(fJ1->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));
	// 		hFatJet_mass->Fill(fJ1->mass());
	// 	} //closes loop through fatJets entries for the event 








      }  
      // close input file
      inFile->Close();
    }
    // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }



// WriteHistograms(outputFile_.c_str());




  return 0;
}

