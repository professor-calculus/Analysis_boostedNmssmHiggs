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
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

// Headers from this package
#include "Analysis/Analysis_boostedNmssmHiggs/interface/Kinematics.h"
#include "Analysis/Analysis_boostedNmssmHiggs/interface/PlottingDoubleBTaggerEfficiencyStudies.h"

// preliminary running, compile with scram b and then
// $ DoubleBTaggerEfficiencyStudies inputfiles=XYZ outputfile=ABC orderedsecondaryfiles=0
// ($CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies/DoubleBTaggerEfficiencyStudies)
// nb. if you are running it on DICE, include the word runOnDice at the end of the arguments of the executable

/* Notes on runOnDice mode
watchout...with this toggle the executable can now overwrite outputs!!!
hacks the original script abit so it can do some things slightly differently
it does so using the quantity: bool runOnDice
*/


void CreateHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::vector<std::string>, std::vector<double>);
void FillHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, bool, pat::Jet, reco::GenParticle, double, std::vector<std::string>, std::vector<double>, std::vector<double>);
void WriteHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::string);
void WriteHistogramsDICE(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::string);
std::vector<reco::GenParticle> higgsBbGenParticles(edm::Handle<std::vector<reco::GenParticle>>, std::vector<double> &);
bool isThereAFatJetMatch(edm::Handle<std::vector<pat::Jet>>, reco::GenParticle, double, pat::Jet&);
std::string getOutputDirFromOutputFile(std::string);


int main(int argc, char* argv[]) 
{
	gSystem->Load("libFWCoreFWLite.so");
	FWLiteEnabler::enable();

	////////////////////
	////////////////////
	// Set parameters //
	std::vector<double> doubleBtagWP = {0.3, 0.6, 0.8, 0.9}; // these WP vectors must correspond to one-another
	std::vector<std::string> doubleBtagWPname = {"loose", "medium", "tight", "veryTight"};
	double dRMaxMatch = 0.8; // max dR between higgs boson and fatJet to claim a match
	// std::vector<double> etaBinning = {0.00, 0.80, 1.60, 2.40};
	std::vector<double> etaBinning = {0.00, 2.40};
	/////////////////////
	/////////////////////

	// see if it is in runOnDice mode
	bool runOnDice = false;
    for (int i = 0; i < argc; ++i) {
        std::string argvString(argv[i]);
        if (argvString == "runOnDice") runOnDice = true;
    }

	// Create histograms, they are accessed by eg: h_["fatJetMass_loose"]->Fill(125.0);
	std::map<std::string, TH1F*> h_;
	std::map<std::string, TH2F*> h2_;
	CreateHistograms(h_, h2_, doubleBtagWPname, etaBinning);

	// Initialize command line parser
	optutl::CommandLineParser parser ("Analyze DoubleBTagger Efficiencies");

	//////////////////
	// Set defaults //
	parser.integerValue ("maxevents"      ) = -1; // -1 for all events
	parser.integerValue ("outputevery"    ) =   100;
	// parser.stringVector ("inputfiles"    ) = {"/hdfs/user/jt15104/Analysis_boostedNmssmHiggs/patTuples/CMSSW_8_0_20/signalSamples/nmssmSignalCascadeV05_13TeV_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/nmssmSignalCascadeV05_13TeV_patTupleAddBTag_ed12_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/bTagPatTuple_888.root"};
	// parser.stringValue  ("outputfile"     ) = "output_DoubleBTaggerEfficiencyStudies/histos.root";
	parser.boolValue    ("orderedsecondaryfiles") = false;
	//////////////////

	// Parse arguments_
	if (runOnDice) parser.parseArguments (argc-1, argv);
	else parser.parseArguments (argc, argv);
	int maxEvents_ = parser.integerValue("maxevents");
	unsigned int outputEvery_ = parser.integerValue("outputevery");
	std::vector<std::string> inputFiles_ = parser.stringVector("inputfiles");
	std::string outputFile_ = parser.stringValue("outputfile");
	bool justDoPlotting_ = parser.boolValue("orderedsecondaryfiles");

	if (runOnDice==false){ // don't do this if on DICE as it cannot do this command locally

		std::string outputDirectory_ = getOutputDirFromOutputFile(outputFile_);
		bool makeDir = !(std::system(Form("mkdir %s",outputDirectory_.c_str())));
		if (justDoPlotting_ == false && makeDir == false){
			std::cout << "The chosen output directory already exists, or the parent directory does not exist:" << std::endl;
			std::cout << "Do not wish to overwrite ROOT file: Exiting..." << std::endl;
			return 1;
		}
		
		// copy the code used to make the histogram ROOT file into the same directory (this could be out of sync if you edit after compilation)
		// also copy the parser values that were used to a .txt file (ie the input data used)
		if (justDoPlotting_ == false){
			
			std::system(Form("cp $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies.cpp %s",outputDirectory_.c_str()));
			
			std::ofstream parserRecord;
			parserRecord.open(Form("%soptionsUsed.txt",outputDirectory_.c_str()));
			parserRecord << "MaxEvents: " << maxEvents_ << "\n";
			parserRecord << "InputFiles: " << "\n";
			for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
				parserRecord << inputFiles_[iFile] << "\n";
			}
			parserRecord.close();
		}
	} // closes 'if' runOnDice


	// Loop through the input files
	int ievt=0;
	// but we are just doing the plotting goto the end
	if (justDoPlotting_) goto plottingLabel;

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
				std::vector<double> dRbbVec; // gets updated by 'higgsBbGenParticles' function
				std::vector<reco::GenParticle> higgsBbParticles = higgsBbGenParticles(genParticles, dRbbVec); 
				// Loop through the higgsBb particles
				for (size_t iH = 0; iH < higgsBbParticles.size(); ++iH){

					const reco::GenParticle higgsBbParticle = higgsBbParticles[iH];

					// See if there is a fatJet that matches to the higgsBb (closest in dR, must have dR<dRMaxMatch)
					pat::Jet fatJetMatch; // if there is a matching fatJet, this object will contain it
					bool isMatch =isThereAFatJetMatch(fatJets, higgsBbParticle, dRMaxMatch, fatJetMatch);
					FillHistograms(h_, h2_, isMatch, fatJetMatch, higgsBbParticle, dRbbVec[iH] ,doubleBtagWPname, doubleBtagWP, etaBinning);

				} // closes loop through higgsBb Particles
				// ------------------------------------------------------------------------------------------------------------//

			} // closes loop through events for this file
			inFile->Close();
	    } // closes 'if' the file exists

	    // break loop if maximal number of events is reached:
	    // this has to be done twice to stop the file loop as well
	    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;

	} // closes loop through files


	// don't do the .pdf plotting on DICE (we will want to hadd outputs first)
	if (runOnDice==true){ 
		WriteHistogramsDICE(h_, h2_, outputFile_.c_str());
		return 0;
	}  

	WriteHistograms(h_, h2_, outputFile_.c_str());
	plottingLabel:
	if (std::system(Form("test -e %s",outputFile_.c_str())) == 0) // check the file exists
		PlottingDoubleBTaggerEfficiencyStudies createPlots(outputFile_.c_str(), doubleBtagWPname, etaBinning);
	else{
		std::cout << "the following output root file does not exist to make plots from:" << std::endl;
		std::cout << outputFile_ << std::endl;
	}

return 0;
} // closes the 'main' function






void CreateHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_, std::vector<std::string> doubleBtagWPnameD, const std::vector<double> etaBinningD)
{
	// set the binning for histograms
    std::vector<double> ptBinning;
    // for(double binLowerEdge=  0.0; binLowerEdge< 100.0; binLowerEdge+= 50.0) ptBinning.push_back(binLowerEdge);
    // for(double binLowerEdge=  100.0; binLowerEdge< 150.0; binLowerEdge+= 50.0) ptBinning.push_back(binLowerEdge);
    // for(double binLowerEdge=  150.0; binLowerEdge< 350.0; binLowerEdge+= 25.0) ptBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  0.0; binLowerEdge< 1000.0; binLowerEdge+= 50.0) ptBinning.push_back(binLowerEdge);    	
    for(double binLowerEdge=  1000.0; binLowerEdge< 1800.1; binLowerEdge+= 100.0) ptBinning.push_back(binLowerEdge);
    // for(double binLowerEdge=  600.0; binLowerEdge< 800.1; binLowerEdge+= 200.0) ptBinning.push_back(binLowerEdge);

    std::vector<double> ptScatXBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 1800.1; binLowerEdge+= 5.0) ptScatXBinning.push_back(binLowerEdge);

    std::vector<double> ptScatYBinning = ptScatXBinning;
    // for(double binLowerEdge=  0.0; binLowerEdge< 800.1; binLowerEdge+= 20.0) ptScatYBinning.push_back(binLowerEdge);	

    std::vector<double> massBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 200.1; binLowerEdge+= 5.0) massBinning.push_back(binLowerEdge);

    std::vector<double> etaDistBinning;
    for(double binLowerEdge=  -4.00; binLowerEdge< 4.0001; binLowerEdge+= 0.20) etaDistBinning.push_back(binLowerEdge);    

    std::vector<double> matchDeltaRDistBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 0.8001; binLowerEdge+= 0.02) matchDeltaRDistBinning.push_back(binLowerEdge); 

    std::vector<double> bbDeltaRDistBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 2.501; binLowerEdge+= 0.05) bbDeltaRDistBinning.push_back(binLowerEdge); 

    std::vector<double> dREffBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 0.40; binLowerEdge+= 0.10) dREffBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  0.40; binLowerEdge< 0.70; binLowerEdge+= 0.10) dREffBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  0.70; binLowerEdge< 1.00; binLowerEdge+= 0.10) dREffBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  1.00; binLowerEdge< 1.10; binLowerEdge+= 0.10) dREffBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  1.10; binLowerEdge< 1.70; binLowerEdge+= 0.20) dREffBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  1.70; binLowerEdge< 2.501; binLowerEdge+= 0.40) dREffBinning.push_back(binLowerEdge);

	// create the debugging histograms
	h_["DEBUG_higgsBbDRpreMatching"] = new TH1F("DEBUG_higgsBbDRpreMatching", ";dR_bb;a.u.", 50, 0, 2.50);

    // create the histograms (includes histograms with a fat jet mass cut)
	for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPnameD.size(); ++iWP){

		h2_[Form("ptScatter_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())] = new TH2F(
			Form("ptScatter_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str()), ";MC Higgs p_{T} (GeV); doubleBTagJet p_{T} (GeV)",
			ptScatXBinning.size()-1, &(ptScatXBinning)[0], ptScatYBinning.size()-1, &(ptScatYBinning)[0]);

		h_[Form("fatJetMass_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())] = new TH1F(
		   Form("fatJetMass_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str()), ";mass_doubleBTagJet (GeV); a.u.", massBinning.size()-1, &(massBinning)[0]);

		h_[Form("fatJetEta_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())] = new TH1F(
		   Form("fatJetEta_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str()), ";doubleBTagJet #eta; a.u.", etaDistBinning.size()-1, &(etaDistBinning)[0]);

		h_[Form("matchDeltaR_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())] = new TH1F(
		   Form("matchDeltaR_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str()), ";dR_match; a.u.", matchDeltaRDistBinning.size()-1, &(matchDeltaRDistBinning)[0]);

		h_[Form("bbDeltaR_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())] = new TH1F(
		   Form("bbDeltaR_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str()), ";MC dR_bb; a.u.", bbDeltaRDistBinning.size()-1, &(bbDeltaRDistBinning)[0]);

   		for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinningD.size()-1; ++iEtaBin){

   			h_[Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )] = new TH1F(
   			   Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] ),
   			   ";MC Higgs p_{T} (GeV); matchCount", ptBinning.size()-1, &(ptBinning)[0]);


   			if (iWP==0) // so we only create the denominators once
   			h_[Form("effDenominator_eta%.2f-%.2f", etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )] = new TH1F(
   			   Form("effDenominator_eta%.2f-%.2f", etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] ),
   			   ";MC Higgs p_{T} (GeV); totalCount", ptBinning.size()-1, &(ptBinning)[0]); 

   			// same as above but as a function of dr rather than pt 
   			h_[Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f_fcnDR", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )] = new TH1F(
   			   Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f_fcnDR", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] ),
   			   ";MC dR_bb; matchCount", dREffBinning.size()-1, &(dREffBinning)[0]);

   			if (iWP==0) // so we only create the denominators once
   			h_[Form("effDenominator_eta%.2f-%.2f_fcnDR", etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )] = new TH1F(
   			   Form("effDenominator_eta%.2f-%.2f_fcnDR", etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] ),
   			   ";MC dR_bb; totalCount", dREffBinning.size()-1, &(dREffBinning)[0]); 

   		} // closes loop through etaBins
	} // closes loop through Btag WPs 
} //closes the function 'CreateHistograms'







void FillHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_, bool isMatch, pat::Jet fatJetMatchD, reco::GenParticle higssBbGenParticleD, double dRbb, std::vector<std::string> doubleBtagWPnameD, std::vector<double> doubleBtagWPD, std::vector<double> etaBinningD)
{


	// fill the efficiency denominators
	for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinningD.size()-1; ++iEtaBin){

		if( fabs(higssBbGenParticleD.eta()) >= etaBinningD[iEtaBin] && fabs(higssBbGenParticleD.eta()) < etaBinningD[iEtaBin+1]){
			h_[Form("effDenominator_eta%.2f-%.2f", etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )]->Fill(higssBbGenParticleD.pt());
			h_[Form("effDenominator_eta%.2f-%.2f_fcnDR", etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )]->Fill(dRbb);
		} // closes 'if' eta within the set bin
	}

	h_["DEBUG_higgsBbDRpreMatching"]->Fill(dRbb);

	// fill other histograms if there is a match
	if (isMatch){

		for (std::vector<std::string>::size_type iWP=0; iWP<doubleBtagWPnameD.size(); ++iWP){
			if (fatJetMatchD.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > doubleBtagWPD[iWP]){

				h2_[Form("ptScatter_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())]->Fill(higssBbGenParticleD.pt(), fatJetMatchD.pt());
				h_[Form("fatJetMass_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())]->Fill(fatJetMatchD.mass());
				h_[Form("fatJetEta_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())]->Fill(fatJetMatchD.eta());
				h_[Form("bbDeltaR_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())]->Fill(dRbb);
				double dRmatch = delR( delPhi( fatJetMatchD.phi(),higssBbGenParticleD.phi() ), delEta( fatJetMatchD.eta(),higssBbGenParticleD.eta() ) );
				h_[Form("matchDeltaR_%sDoubleBTagWP", doubleBtagWPnameD[iWP].c_str())]->Fill(dRmatch);				

		   		for (std::vector<double>::size_type iEtaBin=0; iEtaBin<etaBinningD.size()-1; ++iEtaBin){

		   			if( fabs(higssBbGenParticleD.eta()) >= etaBinningD[iEtaBin] && fabs(higssBbGenParticleD.eta()) < etaBinningD[iEtaBin+1]){ //this is where the eta binning is used
						h_[Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )]->Fill(higssBbGenParticleD.pt());
						h_[Form("effNumerator_%sDoubleBTagWP_eta%.2f-%.2f_fcnDR", doubleBtagWPnameD[iWP].c_str(), etaBinningD[iEtaBin], etaBinningD[iEtaBin+1] )]->Fill(dRbb);
					} // closes 'if' the eta within the set bin
				} // closes loop through etaBins

	   		} // closes 'if' Btag discriminator greater than WP 
		} // closes loop through Btag WPs
	} // closes 'if' there is a matching fatJet

} // closes the function 'FillHistograms'







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







void WriteHistogramsDICE(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_, std::string outputFileD)
{
	TFile outFile(outputFileD.c_str(),"RECREATE");
	for ( auto & h : h_ ){
		h.second->Write();
	}
	for ( auto & h : h2_ ){
		h.second->Write();
	}
	outFile.Close();
	std::system("cp	*.root	../../.	"); // get the output back on soolin
}








std::vector<reco::GenParticle> higgsBbGenParticles(edm::Handle<std::vector<reco::GenParticle>> genParticles, std::vector<double> & dRbbVec)
{
	// this is a vector containing the h->bb higgs'
	std::vector<reco::GenParticle> hBbGenParticles;

	for (size_t iGen = 0; iGen < genParticles->size(); ++iGen){

		const reco::GenParticle & genParticle = (*genParticles)[iGen];

			if (genParticle.pdgId()==35){ // particle is a higgs						
				if (genParticle.numberOfDaughters()==2 && abs(genParticle.daughter(0)->pdgId())==5 && abs(genParticle.daughter(1)->pdgId())==5){ // higgs decays to two b-quarks

					hBbGenParticles.push_back(genParticle);
					dRbbVec.push_back( delR( delPhi( genParticle.daughter(0)->phi(),genParticle.daughter(1)->phi() ), delEta( genParticle.daughter(0)->eta(),genParticle.daughter(1)->eta() ) ) );

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







std::string getOutputDirFromOutputFile(std::string outputFile)
{
	std::string forwardSlash = "/";
	std::string outputDirectory = outputFile;
	// strip the directory from the outputfile name
	for (size_t c = outputFile.size()-1; c > 0; --c){
		if (outputFile[c] == forwardSlash[0]){
			outputDirectory = outputFile.substr(0, c+1);
			break;
		}
	}
	return outputDirectory;
}
