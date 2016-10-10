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
#include "Analysis/Analysis_boostedNmssmHiggs/interface/Kinematics.h"

// command to make work!
// ~/CMSSW_8_0_20/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies/DoubleBTaggerEfficiencyStudies inputFiles=/users/jt15104/CMSSW_8_0_20/src/Analysis/Analysis_boostedNmssmHiggs/python/bTagPatTuple_testingv1.root 

void CreateHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, const std::vector<std::string>, const std::vector<double>);
void FillHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, bool, pat::Jet, reco::GenParticle, std::vector<std::string>, std::vector<double>, std::vector<double>);
void WriteHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, const std::string);
std::vector<reco::GenParticle> higgsBbGenParticles(edm::Handle<std::vector<reco::GenParticle>>);
bool isThereAFatJetMatch(edm::Handle<std::vector<pat::Jet>>, reco::GenParticle, double, pat::Jet&);




int main(int argc, char* argv[]) 
{
	gSystem->Load("libFWCoreFWLite");//update THIS
	AutoLibraryLoader::enable();
	
	// Set parameters
	const std::vector<double> doubleBtagWP = {0.3, 0.6, 0.8, 0.9};
	const std::vector<std::string> doubleBtagWPname = {"loose", "medium", "tight", "veryTight"};
	const double dRMaxMatch = 0.8; // max dR between higgs boson and fatJet to claim a match
	const std::vector<double> etaBinning = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};

	// Create histograms, they are accessed by eg: h_["fatJetMass_loose"]->Fill(125.0);
	std::map<std::string, TH1F*> h_;
	std::map<std::string, TH2F*> h2_;
	CreateHistograms(h_, h2_, doubleBtagWPname, etaBinning);

	// Initialize command line parser
	optutl::CommandLineParser parser ("Analyze DoubleBTagger Efficiencies");

	// Set defaults
	parser.integerValue ("maxEvents"  ) = -1;
	parser.integerValue ("outputEvery") =   10;
	parser.stringValue  ("outputFile" ) = "output_DoubleBTaggerEfficiencyStudies.root";
	// parser.stringValue  ("inputFiles" ) = {"../python/bTagPatTuple_testingv1.root"};
// I kind of want an output directory...and the name of this controls what is going on...WORK

	// Parse arguments
	parser.parseArguments (argc, argv);
	int maxEvents_ = parser.integerValue("maxEvents");
	unsigned int outputEvery_ = parser.integerValue("outputEvery");
	std::string outputFile_ = parser.stringValue("outputFile");
	std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");



	// Loop through the input files
	int ievt=0;  
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
				if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false) 
				std::cout << "  processing event: " << ievt << std::endl;


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
					const reco::GenParticle higgsBbParticle = higgsBbParticles[iH];

					// See if there is a fatJet that matches to the higgsBb (closest in dR, must have dR<dRMaxMatch)
					pat::Jet fatJetMatch; // if there is a matching fatJet, this object will contain it
					bool isMatch =isThereAFatJetMatch(fatJets, higgsBbParticle, dRMaxMatch, fatJetMatch);


					FillHistograms(h_, h2_, isMatch, fatJetMatch, higgsBbParticle, doubleBtagWPname, doubleBtagWP, etaBinning);

				} // closes loop through higgsBb Particles
				// ------------------------------------------------------------------------------------------------------------//


			} // closes loop through events for this file
			inFile->Close();
	    } // closes 'if' the file exists

	    // break loop if maximal number of events is reached:
	    // this has to be done twice to stop the file loop as well
	    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;

	} // closes loop through files

WriteHistograms(h_, h2_, outputFile_.c_str());

return 0;
} // closes the 'main' function






void CreateHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_, const std::vector<std::string> doubleBtagWPnameD, const std::vector<double> etaBinningD)
{
	// set the binning for histograms
    std::vector<double> ptBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 800.1; binLowerEdge+= 20.0) ptBinning.push_back(binLowerEdge);

    std::vector<double> ptScatXBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 800.1; binLowerEdge+= 20.0) ptScatXBinning.push_back(binLowerEdge);

    std::vector<double> ptScatYBinning = ptScatXBinning;
    // for(double binLowerEdge=  0.0; binLowerEdge< 800.1; binLowerEdge+= 20.0) ptScatYBinning.push_back(binLowerEdge);	

    std::vector<double> massBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 200.1; binLowerEdge+= 10.0) massBinning.push_back(binLowerEdge);


    // create the histograms
	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPnameD.size(); ++iWP){

		h2_[Form("ptScatter_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())] = new TH2F(
			Form("ptScatter_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str()), ";higgs p_{T} (GeV); fatJet p_{T} (GeV)",
			ptScatXBinning.size()-1, &(ptScatXBinning)[0], ptScatYBinning.size()-1, &(ptScatYBinning)[0]);

		h_[Form("fatJetMass_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())] = new TH1F(
		   Form("fatJetMass_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str()), ";m_{Jet} (GeV); a.u.", massBinning.size()-1, &(massBinning)[0]);

   		for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinningD.size()-1; ++iEtaBin){

   			h_[Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )] = new TH1F(
   			   Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] ),
   			   ";p_{T} (GeV); efficiency", ptBinning.size()-1, &(ptBinning)[0]); 

   			if (iWP==0) // so we only create the denominators once
   			h_[Form("effDenominator_eta%.2f-%.2f", etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )] = new TH1F(
   			   Form("effDenominator_eta%.2f-%.2f", etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] ),
   			   ";p_{T} (GeV); efficiency", ptBinning.size()-1, &(ptBinning)[0]); 


   		} // closes loop through etaBins
	} // closes loop through Btag WPs 

} //closes the function 'CreateHistograms'






void FillHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_, bool isMatch, pat::Jet fatJetMatchD, reco::GenParticle higssBbGenParticleD, std::vector<std::string> doubleBtagWPnameD, std::vector<double> doubleBtagWPD, std::vector<double> etaBinningD)
{
	// fill the efficiency denominators
	for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinningD.size()-1; ++iEtaBin){

		if( abs(fatJetMatchD.eta()) >= etaBinningD[iEtaBin] && abs(higssBbGenParticleD.eta()) < etaBinningD[iEtaBin+1])
		h_[Form("effDenominator_eta%.2f-%.2f", etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )]->Fill(higssBbGenParticleD.pt());

	}

	// fill other histograms if there is a match
	if (isMatch){

		for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPnameD.size(); ++iWP){
			if (fatJetMatchD.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > doubleBtagWPD[iWP]){

				h2_[Form("ptScatter_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())]->Fill(higssBbGenParticleD.pt(), fatJetMatchD.pt());
				h_[Form("fatJetMass_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())]->Fill(fatJetMatchD.mass());

		   		for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinningD.size()-1; ++iEtaBin){
		   		
		   			if( abs(fatJetMatchD.eta()) >= etaBinningD[iEtaBin] && abs(higssBbGenParticleD.eta()) < etaBinningD[iEtaBin+1])
					h_[Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )]->Fill(higssBbGenParticleD.pt());
				
				} // closes loop through etaBins

	   		} // closes 'if' Btag discriminator greater than WP 
		} // closes loop through Btag WPs
	} // closes 'if' there is a matching fatJet

} // closes the function 'FillHistograms'






void WriteHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_, const std::string outputFileD)
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






std::vector<reco::GenParticle> higgsBbGenParticles(edm::Handle<std::vector<reco::GenParticle>> genParticles)
{
	std::vector<reco::GenParticle> hBbGenParticles;
	for (size_t iGen = 0; iGen < genParticles->size(); ++iGen){
		const reco::GenParticle & genParticle = (*genParticles)[iGen];

			if (genParticle.pdgId()==25){ // particle is a higgs						
				if (genParticle.numberOfDaughters()==2 && abs(genParticle.daughter(0)->pdgId())==5 && abs(genParticle.daughter(1)->pdgId())==5){ // higgs decays to two b-quarks
							
					hBbGenParticles.push_back(genParticle); // this vector indexes the higgs bosons of interest	

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