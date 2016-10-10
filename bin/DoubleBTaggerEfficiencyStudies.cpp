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



// function: calculate dPhi for two phi inputs
// NB:this is phi1-phi2
double calc_dPhi(double phi1, double phi2)
{
	double dPhi = phi1 - phi2;
	if (dPhi>M_PI) dPhi = dPhi - 2*M_PI;
	if (dPhi<-M_PI) dPhi = dPhi + 2*M_PI;
	return dPhi;
}
  
// function: calculate dEta for two eta inputs
// NB:this is eta1-eta2
double calc_dEta(double eta1, double eta2)
{
	double dEta = eta1 - eta2;
	return dEta;
}
    
// function: calculates dR given the dPhi and dEta
double calc_dR(double dPhi, double dEta)
{
  double dR = sqrt(dPhi*dPhi + dEta*dEta);
  return dR;
}





// need to take eta binning out of here
void CreateHistograms(const std::vector<std::string> doubleBtagWPnameD, const std::vector<double> etaBinningD,
	                  std::map<std::string, TH1F*> & h_, std::map<std::string, TH2F*> & h2_)
{

    std::vector<double> ptBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 800.1; binLowerEdge+= 20.0) ptBinning.push_back(binLowerEdge);

    std::vector<double> ptScatXBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 800.1; binLowerEdge+= 20.0) ptScatXBinning.push_back(binLowerEdge);

    std::vector<double> ptScatYBinning = ptScatXBinning;
    // for(double binLowerEdge=  0.0; binLowerEdge< 800.1; binLowerEdge+= 20.0) ptScatYBinning.push_back(binLowerEdge);	

    std::vector<double> massBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 200.1; binLowerEdge+= 10.0) massBinning.push_back(binLowerEdge);


	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPnameD.size(); ++iWP){

		h2_[Form("ptScatter_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())] = new TH2F(
			Form("ptScatter_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str()), ";higgs p_{T} (GeV); fatJet p_{T} (GeV)",
			ptScatXBinning.size()-1, &(ptScatXBinning)[0], ptScatYBinning.size()-1, &(ptScatYBinning)[0]);

		h_[Form("fatJetMass_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())] = new TH1F(
		   Form("fatJetMass_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str()), ";m_{Jet} (GeV); a.u.", massBinning.size()-1, &(massBinning)[0]);

   		for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinningD.size()-1; ++iEtaBin){

   			h_[Form("effNumerator_%sDoubleBTagWP_eta%f-%f", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )] = new TH1F(
   			   Form("effNumerator_%sDoubleBTagWP_eta%f-%f", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] ),
   			   ";p_{T} (GeV); efficiency", ptBinning.size()-1, &(ptBinning)[0]); 

   			h_[Form("effDenominator_%sDoubleBTagWP_eta%f-%f", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )] = new TH1F(
   			   Form("effDenominator_%sDoubleBTagWP_eta%f-%f", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] ),
   			   ";p_{T} (GeV); efficiency", ptBinning.size()-1, &(ptBinning)[0]); 



   		} // closes loop through etaBins
	} // closes loop through Btag WPs 

} //closes CreateHistograms() function


void WriteHistograms(const std::string outputFileD, std::map<std::string, TH1F*> & h_, std::map<std::string, TH2F*> & h2_)
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
	const double dRMaxMatch = 0.8; // dR between higgs boson and bTaggedFatJet
	const std::vector<double> etaBinning = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0}; // need to put this in the create function prooperly

	// Create histograms, they are accessed by eg: h_["fatJetMass_loose"]->Fill(125.0);
	std::map<std::string, TH1F*> h_;
	std::map<std::string, TH2F*> h2_;
	CreateHistograms(doubleBtagWPname, etaBinning, h_, h2_);


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




// could have a function that returns the indices of the higss->bbs...
				// Loop through genParticle collection to find the h->bb
				// std::vector<size_t> hBbGenIndices;
				for (size_t iGen = 0; iGen < genParticles->size(); ++iGen){
				const reco::GenParticle & genParticle = (*genParticles)[iGen];

					if (genParticle.pdgId()==25){ // particle is a higgs
						
						if (genParticle.numberOfDaughters()==2 && abs(genParticle.daughter(0)->pdgId())==5 && abs(genParticle.daughter(1)->pdgId())==5){ // higgs decays to two b-quarks
							
							// hBbGenIndices.push_back(iGen); // this vector indexes the higgs bosons of interest


							// Loop through fatJet collection
							// find closest matching jet to the higgs, check it is a match, see if it passes btag WP cuts
							


// could have a function that takes the jatJets object and finds the closest match, returning the closest full fatJet object...:)
// as long as it is less than dRMAXMATCH, if it doesn't you could just return a 0


							size_t closestFatJetIndex = 999;
							double dRMin = 999.99;
							for (size_t iFJ = 0; iFJ < fatJets->size(); ++iFJ){
								const pat::Jet & fatJet = (*fatJets)[iFJ];

								double dR = calc_dR( calc_dPhi( fatJet.phi(),genParticle.phi() ), calc_dEta( fatJet.eta(),genParticle.eta() ) );
								if (dR < dRMin){

									dRMin = dR;
									closestFatJetIndex = iFJ; 
								}

							} // closes loop through fatJets entries for the event 

							//now do histogram filling...
							if (dRMin < dRMaxMatch){

const pat::Jet & fatJetMatch = (*fatJets)[closestFatJetIndex];
// // this could also be a function...call it fill histograms
	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPname.size(); ++iWP){

		if (fatJetMatch.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > doubleBtagWP[iWP]){


			std::cout << genParticle.pt() << "   " << fatJetMatch.mass() << std::endl;


		h2_[Form("ptScatter_%sDoubleBTagWP", doubleBtagWPname[iWP].c_str())]->Fill(genParticle.pt(), fatJetMatch.pt()); // check two-d histo filling!!

		h_[Form("fatJetMass_%sDoubleBTagWP", doubleBtagWPname[iWP].c_str())]->Fill(fatJetMatch.mass());

   		for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinning.size()-1; ++iEtaBin){

   			if( abs(fatJetMatch.eta()) >= etaBinning[iEtaBin] && abs(fatJetMatch.eta()) < etaBinning[iEtaBin+1]){

   			h_[Form("eff_%sDoubleBTagWP_eta%f-%f", doubleBtagWPname[iWP].c_str(), etaBinning[iEtaBin], etaBinning[iEtaBin+1] )]->Fill(fatJetMatch.pt());

} // if eta in the right bin like
}// closes loop through etaBins
   		} //if btag discriminator 
	} // closes loop through Btag WPs




// // need to do efficiencies properly, with a denom and numerator...deal with them later...


							}//if dr match

					


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