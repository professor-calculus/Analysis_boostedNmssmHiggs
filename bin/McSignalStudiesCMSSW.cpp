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
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

// Headers from this package
#include "Analysis/Analysis_boostedNmssmHiggs/interface/Kinematics.h"
#include "Analysis/Analysis_boostedNmssmHiggs/interface/PlottingMcSignalStudiesCMSSW.h"

// TODO: use ak4pfCHS in future when it is available!!!

// PAT TUPLE FORMAT ONLY !!!

// preliminary running, compile with scram b and then
// $ McSignalStudiesCMSSW inputfiles=XYZ outputfile=ABC orderedsecondaryfiles=0
// ($CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/McSignalStudiesCMSSW/McSignalStudiesCMSSW)
// nb. if you are running it on DICE, include the word runOnDice at the end of the arguments of the executable

/* Notes on runOnDice mode
watchout...with this toggle the executable can now overwrite outputs!!!
hacks the original script abit so it can do some things slightly differently
it does so using the quantity: bool runOnDice
*/

// NOTE: the first element is leading arm (decay from highest pt squark)
// NOTE: the second element is for the secondary arm (decay from secondary pt squark)

// NB: *known bug* if running on just plotting mode you still need to give it an inputfiles argument so that the script can get a title without seg faulting!!!

void CreateHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&);
void WriteHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::string);
void WriteHistogramsDICE(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::string);
std::string getOutputDirFromOutputFile(std::string);
bool getCascadeParticles(edm::Handle<std::vector<reco::GenParticle>>,int,std::string,unsigned int&,std::vector<reco::GenParticle>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&,std::vector<const reco::Candidate*>&);


int main(int argc, char* argv[]) 
{
	gSystem->Load("libFWCoreFWLite.so");
	FWLiteEnabler::enable();
	
	////////////////////
	////////////////////
	// SET PARAMETERS //
	std::vector<double> doubleBtagWP = {0.3, 0.6, 0.9}; // these WP's correspond to loose, medium, tight, veryTight
	double dRMaxMatch = 0.5; // max dR between higgs boson and fatJet to claim a match
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
	CreateHistograms(h_, h2_);

	// Initialize command line parser
	optutl::CommandLineParser parser ("McSignalStudiesCMSSW ");
	//////////////////
	// Set defaults //
	parser.integerValue ("maxevents"      ) = -1; // -1 for all events
	parser.integerValue ("outputevery"    ) = 100;
	// parser.stringVector ("inputfiles"   ) = {"/hdfs/user/jt15104/Analysis_boostedNmssmHiggs/patTuples/CMSSW_8_0_20/signalSamples/nmssmSignalCascadeV05_13TeV_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/nmssmSignalCascadeV05_13TeV_patTupleAddBTag_ed12_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/bTagPatTuple_888.root"};
	// parser.stringValue  ("outputfile"     ) = "output_McSignalStudiesCMSSW/histos.root";
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
			
			std::system(Form("cp $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/bin/McSignalStudiesCMSSW.cpp %s",outputDirectory_.c_str()));
			
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

	///////////////////////////////
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

				// Handle to the ak4Jet collection
				edm::Handle<std::vector<pat::Jet>> ak4Jets;
				event.getByLabel(std::string("selectedPatJets"), ak4Jets);
				// event.getByLabel(std::string("selectedPatJetsAK4PFCHS"), ak4Jets); // TODO: use this in future

				// Handle to the fatJet collection
				edm::Handle<std::vector<pat::Jet>> fatJets;
				event.getByLabel(std::string("selectedPatJetsAK8PFCHS"), fatJets);

				// Handle to the MET collection
				edm::Handle<std::vector<pat::MET>> METvec;
				event.getByLabel(std::string("patMETs"), METvec);

				// Handle to the Gen Particle collection
 				edm::Handle<std::vector<reco::GenParticle>> genParticles;
				event.getByLabel(std::string("genParticles"), genParticles);	


				// A N A L Y S I S
				// ------------------------------------------------------------------------------------------------------------//
				

				// first get all the Gen Particles in the cascade
				//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				// count how many gluinos are involved, we know we expect 2 squarks
				unsigned int gluinoCount = 0;
				// the first element is leading arm version of particle(decay from highest pt squark)
				// the second element is for the secondary arm version of particle (decay from secondary pt squark)
				std::vector<reco::GenParticle> squarkVec; // note that this object is slightly different to the others
				std::vector<const reco::Candidate*> qjetVec;
				std::vector<const reco::Candidate*> nlspVec;
				std::vector<const reco::Candidate*> lspVec;
				std::vector<const reco::Candidate*> higgsVec;
				std::vector<const reco::Candidate*> bVec;
				std::vector<const reco::Candidate*> bbarVec;

				bool gotCascadeParticles = getCascadeParticles(genParticles,ievt,inputFiles_[iFile],gluinoCount,squarkVec,qjetVec,nlspVec,lspVec,higgsVec,bVec,bbarVec);
				if (gotCascadeParticles == false) continue;

				//  e.g.
				// std::cout << "leading squark pt: " << squarkVec[0].pt() << "   secondary squark pt: " << squarkVec[1].pt() << std::endl;
				// std::cout << "leading lsp pt: " << lspVec[0]->pt() << "   second lsp pt: " << lspVec[1]->pt();
				//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

				// Gen Particle Plots
				h_["numberOfGluinos"]->Fill(gluinoCount);
				
				h_["leadingSquarkPt"]->Fill(squarkVec[0].pt());
				if (gluinoCount == 0) h_["leadingSquarkPt_zeroGluinos"]->Fill(squarkVec[0].pt());
				if (gluinoCount == 1) h_["leadingSquarkPt_oneGluinos"]->Fill(squarkVec[0].pt());
				if (gluinoCount == 2) h_["leadingSquarkPt_twoGluinos"]->Fill(squarkVec[0].pt());

				h_["leadingSquarkEta"]->Fill(squarkVec[0].eta());
				if (gluinoCount == 0) h_["leadingSquarkEta_zeroGluinos"]->Fill(squarkVec[0].eta());
				if (gluinoCount == 1) h_["leadingSquarkEta_oneGluinos"]->Fill(squarkVec[0].eta());
				if (gluinoCount == 2) h_["leadingSquarkEta_twoGluinos"]->Fill(squarkVec[0].eta());

				h_["secondarySquarkPt"]->Fill(squarkVec[1].pt());
				if (gluinoCount == 0) h_["secondarySquarkPt_zeroGluinos"]->Fill(squarkVec[1].pt());
				if (gluinoCount == 1) h_["secondarySquarkPt_oneGluinos"]->Fill(squarkVec[1].pt());
				if (gluinoCount == 2) h_["secondarySquarkPt_twoGluinos"]->Fill(squarkVec[1].pt());

				h_["secondarySquarkEta"]->Fill(squarkVec[1].eta());
				if (gluinoCount == 0) h_["secondarySquarkEta_zeroGluinos"]->Fill(squarkVec[1].eta());
				if (gluinoCount == 1) h_["secondarySquarkEta_oneGluinos"]->Fill(squarkVec[1].eta());
				if (gluinoCount == 2) h_["secondarySquarkEta_twoGluinos"]->Fill(squarkVec[1].eta());

				h2_["leadingSquarkPt_SecondarySquarkPt"]->Fill(squarkVec[1].pt(),squarkVec[0].pt());
				h2_["leadingSquarkEta_SecondarySquarkEta"]->Fill(squarkVec[1].eta(),squarkVec[0].eta());
				h2_["leadingSquarkPhi_SecondarySquarkPhi"]->Fill(squarkVec[1].phi(),squarkVec[0].phi());

				h_["leadingQjetPt"]->Fill(qjetVec[0]->pt());
				h_["leadingQjetEta"]->Fill(qjetVec[0]->eta());
				h_["secondaryQjetPt"]->Fill(qjetVec[1]->pt());
				h_["secondaryQjetEta"]->Fill(qjetVec[1]->eta());
				h2_["leadingQjetPt_secondaryQjetPt"]->Fill(qjetVec[1]->pt(), qjetVec[0]->pt());
				h2_["leadingQjetEta_secondaryQjetEta"]->Fill(qjetVec[1]->eta(), qjetVec[0]->eta());
				h2_["leadingQjetPhi_secondaryQjetPhi"]->Fill(qjetVec[1]->phi(), qjetVec[0]->phi());

				h_["leadingNlspPt"]->Fill(nlspVec[0]->pt());
				h_["leadingNlspEta"]->Fill(nlspVec[0]->eta());
				h_["secondaryNlspPt"]->Fill(nlspVec[1]->pt());
				h_["secondaryNlspEta"]->Fill(nlspVec[1]->eta());

				h2_["leadingQjetPt_leadingHiggsPt"]->Fill(higgsVec[0]->pt(), qjetVec[0]->pt());
				h2_["leadingQjetEta_leadingHiggsEta"]->Fill(higgsVec[0]->eta(), qjetVec[0]->eta());
				h2_["leadingQjetPhi_leadingHiggsPhi"]->Fill(higgsVec[0]->phi(), qjetVec[0]->phi());
				h2_["secondaryQjetPt_secondaryHiggsPt"]->Fill(higgsVec[1]->pt(), qjetVec[1]->pt());
				h2_["secondaryQjetEta_secondaryHiggsEta"]->Fill(higgsVec[1]->eta(), qjetVec[1]->eta());
				h2_["secondaryQjetPhi_secondaryHiggsPhi"]->Fill(higgsVec[1]->phi(), qjetVec[1]->phi());

				h_["leadingHiggsPt"]->Fill(higgsVec[0]->pt());
				h_["leadingHiggsEta"]->Fill(higgsVec[0]->eta());
				h_["secondaryHiggsPt"]->Fill(higgsVec[1]->pt());
				h_["secondaryHiggsEta"]->Fill(higgsVec[1]->eta());
				h2_["leadingHiggsPt_secondaryHiggsPt"]->Fill(higgsVec[1]->pt(), higgsVec[0]->pt());
				h2_["leadingHiggsEta_secondaryHiggsEta"]->Fill(higgsVec[1]->eta(), higgsVec[0]->eta());
				h2_["leadingHiggsPhi_secondaryHiggsPhi"]->Fill(higgsVec[1]->phi(), higgsVec[0]->phi());
				double dphiLeadingHiggsQjet = delPhi(higgsVec[0]->phi(),qjetVec[0]->phi());
				if (dphiLeadingHiggsQjet < 0) dphiLeadingHiggsQjet += 2*M_PI;
				h_["leadingHiggsQjetDphi"]->Fill(dphiLeadingHiggsQjet);
				double dphiSecondaryHiggsQjet = delPhi(higgsVec[1]->phi(),qjetVec[1]->phi());
				if (dphiSecondaryHiggsQjet < 0) dphiSecondaryHiggsQjet += 2*M_PI;
				h_["secondaryHiggsQjetDphi"]->Fill(dphiSecondaryHiggsQjet);
				double dphiLeadingHiggsSecondaryHiggs = delPhi(higgsVec[0]->phi(),higgsVec[1]->phi());
				if (dphiLeadingHiggsSecondaryHiggs < 0) dphiLeadingHiggsSecondaryHiggs += 2*M_PI;
				h_["leadingHiggsSecondaryHiggsDphi"]->Fill(dphiLeadingHiggsSecondaryHiggs);
				double dphiLeadingQjetSecondaryQjet = delPhi(qjetVec[0]->phi(),qjetVec[1]->phi());
				if (dphiLeadingQjetSecondaryQjet < 0) dphiLeadingQjetSecondaryQjet += 2*M_PI;
				h_["leadingQjetSecondaryQjetDphi"]->Fill(dphiLeadingQjetSecondaryQjet);

				h_["leadingLspPt"]->Fill(lspVec[0]->pt());
				h_["leadingLspEta"]->Fill(lspVec[0]->eta());
				h_["secondaryLspPt"]->Fill(lspVec[1]->pt());
				h_["secondaryLspEta"]->Fill(lspVec[1]->eta());
				//function: adds two pt vectors together in the (2d) transverse plane
				//output: first element is magnitude, second element is phi
				std::vector<double> lspMet = addTwoPtVectors(lspVec[0]->pt(), lspVec[0]->phi(), lspVec[1]->pt(), lspVec[1]->phi());
				h_["lspMET"]->Fill(lspMet[0]);

				double leadingBBbarSeperation = delR( delPhi(bVec[0]->phi(), bbarVec[0]->phi()), delEta(bVec[0]->eta(), bbarVec[0]->eta()) );
				h_["leadingBBbarSeperation"]->Fill(leadingBBbarSeperation);
				double secondaryBBbarSeperation = delR( delPhi(bVec[1]->phi(), bbarVec[1]->phi()), delEta(bVec[1]->eta(), bbarVec[1]->eta()) );
				h_["secondaryBBbarSeperation"]->Fill(secondaryBBbarSeperation);
				h2_["leadingBBbarSeperation_secondaryBBbarSeperation"]->Fill(secondaryBBbarSeperation,leadingBBbarSeperation);

				if (leadingBBbarSeperation <= 0.8 && secondaryBBbarSeperation <= 0.8) h2_["logicBbDr"]->Fill(0.1,0.1);
				if (leadingBBbarSeperation <= 0.8 && secondaryBBbarSeperation > 0.8) h2_["logicBbDr"]->Fill(1.1,0.1);
				if (leadingBBbarSeperation > 0.8 && secondaryBBbarSeperation <= 0.8) h2_["logicBbDr"]->Fill(0.1,1.1);
				if (leadingBBbarSeperation > 0.8 && secondaryBBbarSeperation > 0.8) h2_["logicBbDr"]->Fill(1.1,1.1);

				if (higgsVec[0]->pt() <= 170 && higgsVec[1]->pt() <= 170) h2_["logicHiggsPt"]->Fill(0.1,0.1);
				if (higgsVec[0]->pt() <= 170 && higgsVec[1]->pt() > 170) h2_["logicHiggsPt"]->Fill(1.1,0.1);
				if (higgsVec[0]->pt() > 170 && higgsVec[1]->pt() <= 170) h2_["logicHiggsPt"]->Fill(0.1,1.1);
				if (higgsVec[0]->pt() > 170 && higgsVec[1]->pt() > 170) h2_["logicHiggsPt"]->Fill(1.1,1.1);

				if (higgsVec[0]->pt() <= 170 && leadingBBbarSeperation <= 0.8) h2_["logicBbDrHiggsPt"]->Fill(0.1,0.1);
				if (higgsVec[0]->pt() > 170 && leadingBBbarSeperation <= 0.8) h2_["logicBbDrHiggsPt"]->Fill(1.1,0.1);
				if (higgsVec[0]->pt() <= 170 && leadingBBbarSeperation > 0.8) h2_["logicBbDrHiggsPt"]->Fill(0.1,1.1);
				if (higgsVec[0]->pt() > 170 && leadingBBbarSeperation > 0.8) h2_["logicBbDrHiggsPt"]->Fill(1.1,1.1);
				if (higgsVec[1]->pt() <= 170 && secondaryBBbarSeperation <= 0.8) h2_["logicBbDrHiggsPt"]->Fill(0.1,0.1);
				if (higgsVec[1]->pt() > 170 && secondaryBBbarSeperation <= 0.8) h2_["logicBbDrHiggsPt"]->Fill(1.1,0.1);
				if (higgsVec[1]->pt() <= 170 && secondaryBBbarSeperation > 0.8) h2_["logicBbDrHiggsPt"]->Fill(0.1,1.1);
				if (higgsVec[1]->pt() > 170 && secondaryBBbarSeperation > 0.8) h2_["logicBbDrHiggsPt"]->Fill(1.1,1.1);

				double numHiggsHighPt = 0.0;
				double numHighBbSep = 0.0;
				if (higgsVec[0]->pt() > 170) numHiggsHighPt++;
				if (higgsVec[1]->pt() > 170) numHiggsHighPt++;
				if (leadingBBbarSeperation > 0.8) numHighBbSep++;
				if (secondaryBBbarSeperation > 0.8) numHighBbSep++;
				numHiggsHighPt += 0.01;
				numHighBbSep += 0.01;
				h2_["logicNumberHiggsPtNumberBbDr"]->Fill(numHiggsHighPt,numHighBbSep);


				// Detector Plots //

				// MET
				const pat::MET & MET = (*METvec)[0];
				h_["detectorMET"]->Fill(MET.et());

				// ak4 Jets (and HT, MHT & SUSY variables) nb: they are not in pt order
				unsigned int numAk4JetsOver10GeV = 0;
				unsigned int numAk4JetsOver20GeV = 0;
				unsigned int numAk4JetsOver40GeV = 0;
				unsigned int numAk4JetsOver100GeV = 0;
				double ht = 0.0;
				std::vector<double> mhtVec = {0,0}; // first element magnitude, second element phi
				std::vector<int> indexBiasDeltaPhiJets; // use to index jets valid for biasedDeltaPhi and alphaT calculations
				double epsilon = 0.0; // a component of the alphaT calculation
				double leadingAk4JetPt = 0.0;
				double leadingAk4JetEta = 999.99;
				double secondaryAk4JetPt = 0.0;
				double secondaryAk4JetEta = 999.99;

				for (size_t iJ = 0; iJ < ak4Jets->size(); ++iJ){
					const pat::Jet & ak4Jet = (*ak4Jets)[iJ];

					if (ak4Jet.pt() > 10.0) numAk4JetsOver10GeV++;
					if (ak4Jet.pt() > 20.0) numAk4JetsOver20GeV++;
					if (ak4Jet.pt() > 40.0) numAk4JetsOver40GeV++;
					if (ak4Jet.pt() > 100.0) numAk4JetsOver100GeV++;

					if (ak4Jet.pt() > 40.0 && abs(ak4Jet.eta())<3.0){ // corresponds to the alpha_t HT definition
						ht = ht + ak4Jet.pt();
						//function: adds two pt vectors together in the (2d) transverse plane
						//output: first element is magnitude, second element is phi
						mhtVec = addTwoPtVectors(mhtVec[0], mhtVec[1], ak4Jet.pt(), ak4Jet.phi());
						
						indexBiasDeltaPhiJets.push_back(iJ); // use to index jets valid for biasedDeltaPhi and alphaT calculations
						epsilon = epsilon + ak4Jet.et(); // a component of the alphaT calculation
					}

					if (ak4Jet.pt() > leadingAk4JetPt){
						secondaryAk4JetPt = leadingAk4JetPt;
						leadingAk4JetPt = ak4Jet.pt();
						secondaryAk4JetEta = leadingAk4JetEta;
						leadingAk4JetEta = ak4Jet.eta();
					}
					else if (ak4Jet.pt() > secondaryAk4JetPt){
						secondaryAk4JetPt = ak4Jet.pt();
						secondaryAk4JetEta = ak4Jet.eta();
					}
				} // closes loop through the ak4 jets

				// *** biased delta phi calc ***
				std::vector<double> biasedDeltaPhiVec;
				for (size_t i = 0; i < indexBiasDeltaPhiJets.size(); ++i){	
					const pat::Jet & soloJet = (*ak4Jets)[indexBiasDeltaPhiJets[i]];

					std::vector<double> biasedMhtVec = {0,0}; // first element magnitude, second element phi
					for (size_t j = 0; j < indexBiasDeltaPhiJets.size(); ++j){
						const pat::Jet & otherJet = (*ak4Jets)[indexBiasDeltaPhiJets[j]];
						if (j != i) biasedMhtVec = addTwoPtVectors(biasedMhtVec[0], biasedMhtVec[1], otherJet.pt(), -1*otherJet.phi());
					}
					double biasedDeltaPhi = delPhi(soloJet.phi(),biasedMhtVec[1]);
					if (biasedDeltaPhi > M_PI) biasedDeltaPhi -= 2*M_PI;
					if (biasedDeltaPhi < -1 * M_PI) biasedDeltaPhi += 2*M_PI;
					biasedDeltaPhiVec.push_back(fabs(biasedDeltaPhi));
				}
				int minindex = min_element(biasedDeltaPhiVec.begin(), biasedDeltaPhiVec.end()) - biasedDeltaPhiVec.begin();
				double minBiasedDeltaPhi = biasedDeltaPhiVec[minindex];

				// *** alphaT calc ***
				double alphaT = 999.99; // set to a large value so the event default is to pass cut
				if (indexBiasDeltaPhiJets.size() > 1){
					// *** delta epsilon calc *** (part of the general case alpha_t calculation)
					std::vector<int> bitOnOrOff(indexBiasDeltaPhiJets.size(),1); 
					std::vector<int> sameBitCount(indexBiasDeltaPhiJets.size(),0);
					std::vector<double> deltaEpsilonVec;

					for (unsigned int i = 0; i < pow(2,indexBiasDeltaPhiJets.size()) / 2; ++i){

						double pseudoJetOneEnergy = 0.0;
						double pseudoJetTwoEnergy = 0.0;

						for (size_t j = 0; j < indexBiasDeltaPhiJets.size(); ++j){

							if (bitOnOrOff[j] == 1) pseudoJetOneEnergy += (*ak4Jets)[indexBiasDeltaPhiJets[j]].et();
							if (bitOnOrOff[j] == 0) pseudoJetTwoEnergy += (*ak4Jets)[indexBiasDeltaPhiJets[j]].et();

							sameBitCount[j] = sameBitCount[j]+1;
							if (sameBitCount[j] == pow(2,j)){

								sameBitCount[j] = 0;
								if (bitOnOrOff[j] == 0) bitOnOrOff[j] = 1; 
								else if (bitOnOrOff[j] == 1) bitOnOrOff[j] = 0;
							}
						} // closes loop through jet number index

						deltaEpsilonVec.push_back(fabs(pseudoJetOneEnergy - pseudoJetTwoEnergy));
					} // closes loop through 2^(numJets) / 2 iterations for delta epsilon calc

					int minindex2 = min_element(deltaEpsilonVec.begin(), deltaEpsilonVec.end()) - deltaEpsilonVec.begin();
					double minDeltaEpsilon = deltaEpsilonVec[minindex2];

					alphaT = 0.5 * (epsilon - minDeltaEpsilon) / (sqrt(epsilon * epsilon - mhtVec[0] * mhtVec[0]));
				} // closes 'if' more than one valid jet for the calculation
				// *** end of alphaT calc ***
				
				h_["detectorNumAk4JetsOver10GeV"]->Fill(numAk4JetsOver10GeV);
				h_["detectorNumAk4JetsOver20GeV"]->Fill(numAk4JetsOver20GeV);
				h_["detectorNumAk4JetsOver40GeV"]->Fill(numAk4JetsOver40GeV);
				h_["detectorNumAk4JetsOver100GeV"]->Fill(numAk4JetsOver100GeV);
				h_["detectorLeadingAk4JetPt"]->Fill(leadingAk4JetPt);
				if (ak4Jets->size() > 0) h_["detectorLeadingAk4JetEta"]->Fill(leadingAk4JetEta);
				h_["detectorSecondaryAk4JetPt"]->Fill(secondaryAk4JetPt);
				if (ak4Jets->size() > 1) h_["detectorSecondaryAk4JetEta"]->Fill(secondaryAk4JetEta);
				h2_["detectorLeadingAk4JetPt_detectorSecondaryAk4JetPt"]->Fill(secondaryAk4JetPt,leadingAk4JetPt);
				h2_["detectorHT_detectorSecondaryAk4JetPt"]->Fill(secondaryAk4JetPt,ht);
				h2_["detectorHT_detectorLeadingAk4JetPt"]->Fill(leadingAk4JetPt,ht);
				h_["detectorHT"]->Fill(ht);
				h_["detectorMHT"]->Fill(mhtVec[0]);
				h2_["detectorMHT_detectorMET"]->Fill(MET.et(),mhtVec[0]);
				h_["detectorMinBiasedDeltaPhi"]->Fill(minBiasedDeltaPhi);
				h_["detectorAlphaT"]->Fill(alphaT);
				double mhtOverMet = mhtVec[0] / MET.et();
				h_["detectorMHToverMET"]->Fill(mhtOverMet);

				if (alphaT <= 0.52 && minBiasedDeltaPhi <= 0.50) h2_["detectorLogicAlphaTMinBiasedDeltaPhi"]->Fill(0.1, 0.1);
				if (alphaT > 0.52 && minBiasedDeltaPhi <= 0.50) h2_["detectorLogicAlphaTMinBiasedDeltaPhi"]->Fill(1.1, 0.1);
				if (alphaT <= 0.52 && minBiasedDeltaPhi > 0.50) h2_["detectorLogicAlphaTMinBiasedDeltaPhi"]->Fill(0.1, 1.1);
				if (alphaT > 0.52 && minBiasedDeltaPhi > 0.50) h2_["detectorLogicAlphaTMinBiasedDeltaPhi"]->Fill(1.1, 1.1);

				if (mhtOverMet <= 1.25 && mhtVec[0] <= 130) h2_["detectorLogicMHToverMETmht"]->Fill(0.1,0.1);
				if (mhtOverMet > 1.25 && mhtVec[0] <= 130) h2_["detectorLogicMHToverMETmht"]->Fill(1.1,0.1);
				if (mhtOverMet <= 1.25 && mhtVec[0] > 130) h2_["detectorLogicMHToverMETmht"]->Fill(0.1,1.1);
				if (mhtOverMet > 1.25 && mhtVec[0] > 130) h2_["detectorLogicMHToverMETmht"]->Fill(1.1,1.1);

				if (mhtVec[0] > 130 && mhtOverMet < 1.25 && minBiasedDeltaPhi > 0.50 && alphaT > 0.52) h2_["detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50andAlphaTmoreThan0p52"]->Fill(1.1, 0.1);
				else h2_["detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50andAlphaTmoreThan0p52"]->Fill(0.1,0.1);

				// BTAGZZZZ fatJets //
				unsigned int numberLooseDoubleBTagsNoMatching = 0;
				unsigned int numberMediumDoubleBTagsNoMatching = 0;
				unsigned int numberTightDoubleBTagsNoMatching = 0;


				unsigned int numberLooseDoubleBTagsWithMatching = 0;
				unsigned int numberMediumDoubleBTagsWithMatching = 0;
				unsigned int numberTightDoubleBTagsWithMatching = 0;

				
				double dR_fJhiggs0_min = 99999.99; 
				double dR_fJhiggs1_min = 99999.99;
				unsigned int dR_fJhiggs0_min_index = 99999;
				unsigned int dR_fJhiggs1_min_index = 99999;

				for (size_t iFJ = 0; iFJ < fatJets->size(); ++iFJ){
					const pat::Jet & fatJet = (*fatJets)[iFJ];

					if (fatJet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > doubleBtagWP[0]) numberLooseDoubleBTagsNoMatching++;
					if (fatJet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > doubleBtagWP[1]) numberMediumDoubleBTagsNoMatching++; 
					if (fatJet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > doubleBtagWP[2]) numberTightDoubleBTagsNoMatching++;

					double dR_fJhiggs0 = delR( delPhi( fatJet.phi(),higgsVec[0]->phi() ), delEta( fatJet.eta(),higgsVec[0]->eta() ) );
					double dR_fJhiggs1 = delR( delPhi( fatJet.phi(),higgsVec[1]->phi() ), delEta( fatJet.eta(),higgsVec[1]->eta() ) );

					if (dR_fJhiggs0 < dR_fJhiggs1 && dR_fJhiggs0 < dR_fJhiggs0_min){
						dR_fJhiggs0_min = dR_fJhiggs0;
						dR_fJhiggs0_min_index = iFJ;						
					}						

					if (dR_fJhiggs1 < dR_fJhiggs0 && dR_fJhiggs1 < dR_fJhiggs1_min){
						dR_fJhiggs1_min = dR_fJhiggs1;
						dR_fJhiggs1_min_index = iFJ;						
					}	

				} // closes loop through fat Jets

				if (dR_fJhiggs0_min < dRMaxMatch){
					const pat::Jet & fatJet = (*fatJets)[dR_fJhiggs0_min_index];
					if (fatJet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > doubleBtagWP[0]) numberLooseDoubleBTagsWithMatching++;
					if (fatJet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > doubleBtagWP[1]) numberMediumDoubleBTagsWithMatching++; 
					if (fatJet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > doubleBtagWP[2]) numberTightDoubleBTagsWithMatching++;	
				}

				if (dR_fJhiggs1_min < dRMaxMatch){
					const pat::Jet & fatJet = (*fatJets)[dR_fJhiggs1_min_index];
					if (fatJet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > doubleBtagWP[0]) numberLooseDoubleBTagsWithMatching++;
					if (fatJet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > doubleBtagWP[1]) numberMediumDoubleBTagsWithMatching++; 
					if (fatJet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") > doubleBtagWP[2]) numberTightDoubleBTagsWithMatching++;
				}

				h_["fatJetNumberLooseDoubleBTagsNoMatching"]->Fill(numberLooseDoubleBTagsNoMatching);
				h_["fatJetNumberMediumDoubleBTagsNoMatching"]->Fill(numberMediumDoubleBTagsNoMatching);
				h_["fatJetNumberTightDoubleBTagsNoMatching"]->Fill(numberTightDoubleBTagsNoMatching);
				h_["fatJetNumberLooseDoubleBTagsWithMatching"]->Fill(numberLooseDoubleBTagsWithMatching);
				h_["fatJetNumberMediumDoubleBTagsWithMatching"]->Fill(numberMediumDoubleBTagsWithMatching);
				h_["fatJetNumberTightDoubleBTagsWithMatching"]->Fill(numberTightDoubleBTagsWithMatching);

				if (numberLooseDoubleBTagsWithMatching == 2){
					if ( (*fatJets)[dR_fJhiggs0_min_index].pt() > (*fatJets)[dR_fJhiggs1_min_index].pt()){
						h_["leadingMatchedFatJetPtTwoLooseTagsMatch"]->Fill((*fatJets)[dR_fJhiggs0_min_index].pt());
						h_["secondaryMatchedFatJetPtTwoLooseTagsMatch"]->Fill((*fatJets)[dR_fJhiggs1_min_index].pt());
						h2_["leading_secondaryMatchedFatJetPtTwoLooseTagsMatch"]->Fill((*fatJets)[dR_fJhiggs1_min_index].pt(),(*fatJets)[dR_fJhiggs0_min_index].pt());
					}
					else {
						h_["leadingMatchedFatJetPtTwoLooseTagsMatch"]->Fill((*fatJets)[dR_fJhiggs1_min_index].pt());
						h_["secondaryMatchedFatJetPtTwoLooseTagsMatch"]->Fill((*fatJets)[dR_fJhiggs0_min_index].pt());				
						h2_["leading_secondaryMatchedFatJetPtTwoLooseTagsMatch"]->Fill((*fatJets)[dR_fJhiggs0_min_index].pt(),(*fatJets)[dR_fJhiggs1_min_index].pt());	
					}
				}

				if (numberMediumDoubleBTagsWithMatching == 2){
					if ( (*fatJets)[dR_fJhiggs0_min_index].pt() > (*fatJets)[dR_fJhiggs1_min_index].pt()){
						h_["leadingMatchedFatJetPtTwoMediumTagsMatch"]->Fill((*fatJets)[dR_fJhiggs0_min_index].pt());
						h_["secondaryMatchedFatJetPtTwoMediumTagsMatch"]->Fill((*fatJets)[dR_fJhiggs1_min_index].pt());
						h2_["leading_secondaryMatchedFatJetPtTwoMediumTagsMatch"]->Fill((*fatJets)[dR_fJhiggs1_min_index].pt(),(*fatJets)[dR_fJhiggs0_min_index].pt());
					}
					else {
						h_["leadingMatchedFatJetPtTwoMediumTagsMatch"]->Fill((*fatJets)[dR_fJhiggs1_min_index].pt());
						h_["secondaryMatchedFatJetPtTwoMediumTagsMatch"]->Fill((*fatJets)[dR_fJhiggs0_min_index].pt());					
						h2_["leading_secondaryMatchedFatJetPtTwoMediumTagsMatch"]->Fill((*fatJets)[dR_fJhiggs0_min_index].pt(),(*fatJets)[dR_fJhiggs1_min_index].pt());
					}
				}

				if (numberTightDoubleBTagsWithMatching == 2){
					if ( (*fatJets)[dR_fJhiggs0_min_index].pt() > (*fatJets)[dR_fJhiggs1_min_index].pt()){
						h_["leadingMatchedFatJetPtTwoTightTagsMatch"]->Fill((*fatJets)[dR_fJhiggs0_min_index].pt());
						h_["secondaryMatchedFatJetPtTwoTightTagsMatch"]->Fill((*fatJets)[dR_fJhiggs1_min_index].pt());
						h2_["leading_secondaryMatchedFatJetPtTwoTightTagsMatch"]->Fill((*fatJets)[dR_fJhiggs1_min_index].pt(),(*fatJets)[dR_fJhiggs0_min_index].pt());
					}
					else {
						h_["leadingMatchedFatJetPtTwoTightTagsMatch"]->Fill((*fatJets)[dR_fJhiggs1_min_index].pt());
						h_["secondaryMatchedFatJetPtTwoTightTagsMatch"]->Fill((*fatJets)[dR_fJhiggs0_min_index].pt());						
						h2_["leading_secondaryMatchedFatJetPtTwoTightTagsMatch"]->Fill((*fatJets)[dR_fJhiggs0_min_index].pt(),(*fatJets)[dR_fJhiggs1_min_index].pt());
					}
				}

				// ------------------------------------------------------------------------------------------------------------//
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
	// sketchy hack to get the higgsMass, squarkMass, ratio and mass splitting from the root file path name
	size_t endNameContainer = 8888;
	size_t beginNameContainer = 8888;
	std::string forwardSlash = "/";
	for (size_t c = inputFiles_[0].size()-2; c > 0; --c){
		
		if (inputFiles_[0][c] == forwardSlash[0]) endNameContainer = c; 
		if (inputFiles_[0].substr(c,2) == "mH"){
			beginNameContainer = c;
			break;
		}
	}
	std::string titleName = inputFiles_[0].substr(beginNameContainer,endNameContainer-beginNameContainer);
	std::cout << titleName << std::endl;

	if (std::system(Form("test -e %s",outputFile_.c_str())) == 0) // check the file exists
		PlottingMcSignalStudiesCMSSW createPlots(outputFile_.c_str(), titleName.c_str());
	else{
		std::cout << "the following output root file does not exist to make plots from:" << std::endl;
		std::cout << outputFile_ << std::endl;
	}

return 0;
} // closes the 'main' function







void CreateHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_)
{
	// Gen Particle Histograms
    h_["numberOfGluinos"] = new TH1F("numberOfGluinos", ";Number of Gluinos;a.u.", 4, 0, 4);
	
	h_["leadingSquarkPt"] = new TH1F("leadingSquarkPt", ";squark p_{T} (GeV);a.u.", 100, 0, 2500);
	h_["leadingSquarkEta"] = new TH1F("leadingSquarkEta", ";#eta squark;a.u.", 100, -5, 5);
	h_["secondarySquarkPt"] = new TH1F("secondarySquarkPt", ";squark p_{T} (GeV);a.u.", 100, 0, 2500);
	h_["secondarySquarkEta"] = new TH1F("secondarySquarkEta", ";#eta squark;a.u.", 100, -5, 5);
	h2_["leadingSquarkPt_SecondarySquarkPt"] = new TH2F("leadingSquarkPt_SecondarySquarkPt", ";secondary squark p_{T} (GeV);leading squark p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);
	h2_["leadingSquarkEta_SecondarySquarkEta"] = new TH2F("leadingSquarkEta_SecondarySquarkEta", ";#eta secondary squark;#eta leading squark", 400, -5, 5, 400, -5, 5);
	h2_["leadingSquarkPhi_SecondarySquarkPhi"] = new TH2F("leadingSquarkPhi_SecondarySquarkPhi", ";secondary squark Phi;leading squark Phi", 400, -M_PI, M_PI, 400, -M_PI, M_PI);
	h_["leadingSquarkPt_zeroGluinos"] = new TH1F("leadingSquarkPt_zeroGluinos", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["leadingSquarkPt_oneGluinos"] = new TH1F("leadingSquarkPt_oneGluinos", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["leadingSquarkPt_twoGluinos"] = new TH1F("leadingSquarkPt_twoGluinos", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["leadingSquarkEta_zeroGluinos"] = new TH1F("leadingSquarkEta_zeroGluinos", ";#eta squark;a.u.", 50, -5, 5);
	h_["leadingSquarkEta_oneGluinos"] = new TH1F("leadingSquarkEta_oneGluinos", ";#eta squark;a.u.", 50, -5, 5);
	h_["leadingSquarkEta_twoGluinos"] = new TH1F("leadingSquarkEta_twoGluinos", ";#eta squark;a.u.", 50, -5, 5);
	h_["secondarySquarkPt_zeroGluinos"] = new TH1F("secondarySquarkPt_zeroGluinos", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["secondarySquarkPt_oneGluinos"] = new TH1F("secondarySquarkPt_oneGluinos", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["secondarySquarkPt_twoGluinos"] = new TH1F("secondarySquarkPt_twoGluinos", ";squark p_{T} (GeV);a.u.", 50, 0, 2500);
	h_["secondarySquarkEta_zeroGluinos"] = new TH1F("secondarySquarkEta_zeroGluinos", ";#eta squark;a.u.", 50, -5, 5);
	h_["secondarySquarkEta_oneGluinos"] = new TH1F("secondarySquarkEta_oneGluinos", ";#eta squark;a.u.", 50, -5, 5);
	h_["secondarySquarkEta_twoGluinos"] = new TH1F("secondarySquarkEta_twoGluinos", ";#eta squark;a.u.", 50, -5, 5);

	h_["leadingQjetPt"] = new TH1F("leadingQjetPt", ";quark p_{T} (GeV);a.u.", 100, 0, 2500);
	h_["leadingQjetEta"] = new TH1F("leadingQjetEta", ";#eta quark;a.u.", 100, -5, 5);
	h_["secondaryQjetPt"] = new TH1F("secondaryQjetPt", ";quark p_{T} (GeV);a.u.", 100, 0, 2500);
	h_["secondaryQjetEta"] = new TH1F("secondaryQjetEta", ";#eta quark;a.u.", 100, -5, 5);
	h2_["leadingQjetPt_secondaryQjetPt"] = new TH2F("leadingQjetPt_secondaryQjetPt", ";secondary quark p_{T} (GeV);leading quark p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);
	h2_["leadingQjetEta_secondaryQjetEta"] = new TH2F("leadingQjetEta_secondaryQjetEta", ";#eta secondary quark;#eta leading quark", 400, -5, 5, 400, -5, 5);
	h2_["leadingQjetPhi_secondaryQjetPhi"] = new TH2F("leadingQjetPhi_secondaryQjetPhi", ";secondary quark Phi;leading quark Phi", 400, -M_PI, M_PI, 400, -M_PI, M_PI);

	h_["leadingNlspPt"] = new TH1F("leadingNlspPt", ";NLSP p_{T} (GeV);a.u.", 100, 0, 2500);
	h_["leadingNlspEta"] = new TH1F("leadingNlspEta", ";#eta NLSP;a.u.", 100, -5, 5);
	h_["secondaryNlspPt"] = new TH1F("secondaryNlspPt", ";NLSP p_{T} (GeV);a.u.", 100, 0, 2500);
	h_["secondaryNlspEta"] = new TH1F("secondaryNlspEta", ";#eta NLSP;a.u.", 100, -5, 5);

	h_["leadingHiggsPt"] = new TH1F("leadingHiggsPt", ";higgs p_{T} (GeV);a.u.", 100, 0, 2500);
	h_["leadingHiggsEta"] = new TH1F("leadingHiggsEta", "; #eta higgs;a.u.", 100, -5, 5);
	h_["secondaryHiggsPt"] = new TH1F("secondaryHiggsPt", ";higgs p_{T} (GeV);a.u.", 100, 0, 2500);
	h_["secondaryHiggsEta"] = new TH1F("secondaryHiggsEta", ";#eta higgs;a.u.", 100, -5, 5);
	h2_["leadingQjetPt_leadingHiggsPt"] = new TH2F("leadingQjetPt_leadingHiggsPt", ";higgs p_{T} (GeV);quark p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);
	h2_["leadingQjetEta_leadingHiggsEta"] = new TH2F("leadingQjetEta_leadingHiggsEta", ";#eta higgs;#eta quark", 400, -5, 5, 400, -5, 5);
	h2_["leadingQjetPhi_leadingHiggsPhi"] = new TH2F("leadingQjetPhi_leadingHiggsPhi", ";higgs Phi;quark Phi", 400, -M_PI, M_PI, 400, -M_PI, M_PI);
	h2_["secondaryQjetPt_secondaryHiggsPt"] = new TH2F("secondaryQjetPt_secondaryHiggsPt", ";higgs p_{T};quark p_{T}", 400, 0, 2500, 400, 0, 2500);
	h2_["secondaryQjetEta_secondaryHiggsEta"] = new TH2F("secondaryQjetEta_secondaryHiggsEta", ";#eta higgs;#eta quark", 400, -5, 5, 400, -5, 5);
	h2_["secondaryQjetPhi_secondaryHiggsPhi"] = new TH2F("secondaryQjetPhi_secondaryHiggsPhi", ";higgs Phi;quark Phi", 400, -M_PI, M_PI, 400, -M_PI, M_PI);
	h2_["leadingHiggsPt_secondaryHiggsPt"] = new TH2F("leadingHiggsPt_secondaryHiggsPt", ";secondary higgs p_{T} (GeV);leading higgs p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);
	h2_["leadingHiggsEta_secondaryHiggsEta"] = new TH2F("leadingHiggsEta_secondaryHiggsEta", ";#eta secondary higgs;#eta leading higgs", 400, -5, 5, 400, -5, 5);
	h2_["leadingHiggsPhi_secondaryHiggsPhi"] = new TH2F("leadingHiggsPhi_secondaryHiggsPhi", ";secondary higgs Phi;leading higgs Phi", 400, -M_PI, M_PI, 400, -M_PI, M_PI);
	h_["leadingHiggsQjetDphi"] = new TH1F("leadingHiggsQjetDphi", ";higgs Phi - qjet Phi;a.u.", 100, 0, 2*M_PI);
	h_["secondaryHiggsQjetDphi"] = new TH1F("secondaryHiggsQjetDphi", ";higgs Phi - qjet Phi;a.u.", 100, 0, 2*M_PI);
	h_["leadingHiggsSecondaryHiggsDphi"] = new TH1F("leadingHiggsSecondaryHiggsDphi", ";leading higgs Phi - secondary higgs Phi;a.u.", 100, 0, 2*M_PI);
	h_["leadingQjetSecondaryQjetDphi"] = new TH1F("leadingQjetSecondaryQjetDphi", ";leading qjet Phi - secondary qjet Phi;a.u.", 100, 0, 2*M_PI);

	h2_["logicBbDr"] = new TH2F("logicBbDr", ";secondary arm dR_bb > 0.8;leading arm dR_bb > 0.8", 2, 0, 2, 2, 0, 2);
	h2_["logicHiggsPt"] = new TH2F("logicHiggsPt", ";secondary arm higgs p_{T} > 170 (GeV);leading arm higgs p_{T} > 170 (GeV)", 2, 0, 2, 2, 0, 2);
	h2_["logicBbDrHiggsPt"] = new TH2F("logicBbDrHiggsPt", ";higgs p_{T} > 170 (GeV); dR_bb > 0.8", 2, 0, 2, 2, 0, 2);
	h2_["logicNumberHiggsPtNumberBbDr"] = new TH2F("logicNumberHiggsPtNumberBbDr", ";number of higgs with p_{T} > 170 (GeV);number of dR_bb > 0.8", 3, 0, 3, 3, 0, 3);
	
	h_["leadingLspPt"] = new TH1F("leadingLspPt", ";LSP p_{T} (GeV);a.u.", 100, 0, 100);
	h_["leadingLspEta"] = new TH1F("leadingLspEta", ";#eta LSP;a.u.", 100, -5, 5);
	h_["secondaryLspPt"] = new TH1F("secondaryLspPt", ";LSP p_{T} (GeV);a.u.", 100, 0, 100);
	h_["secondaryLspEta"] = new TH1F("secondaryLspEta", ";#eta LSP;a.u.", 100, -5, 5);
	h_["lspMET"] = new TH1F("lspMET", ";LSP E_{T}^{miss} (GeV);a.u.", 100, 0, 100);

	h_["leadingBBbarSeperation"] = new TH1F("leadingBBbarSeperation", ";dR_bb;a.u.", 100, 0, 2.5);
	h_["secondaryBBbarSeperation"] = new TH1F("secondaryBBbarSeperation", ";dR_bb;a.u.", 100, 0, 2.5);
	h2_["leadingBBbarSeperation_secondaryBBbarSeperation"] = new TH2F("leadingBBbarSeperation_secondaryBBbarSeperation", ";secondary dR_bb;leading dR_bb", 400, 0, 2.5, 400, 0, 2.5);	

	// Detector Histograms
	h_["detectorMET"] = new TH1F("detectorMET", ";detector E_{T}^{miss} (GeV);a.u.", 100, 0, 800);
	h_["detectorHT"] = new TH1F("detectorHT", ";detector H_{T} (GeV);a.u.", 100, 0, 7000);
	h2_["detectorMHT_detectorMET"] = new TH2F("detectorMHT_detectorMET", ";detector E_{T}^{miss} (GeV);detector H_{T}^{miss} (GeV)", 400, 0, 800, 400, 0, 800);
	h_["detectorMHT"] = new TH1F("detectorMHT", ";detector H_{T}^{miss} (GeV);a.u.", 100, 0, 800);
	h_["detectorMinBiasedDeltaPhi"] = new TH1F("detectorMinBiasedDeltaPhi", ";Min Biased Delta Phi; a.u.", 100, 0, M_PI);
	h_["detectorAlphaT"] = new TH1F("detectorAlphaT", ";Alpha_T; a.u.", 100, 0, 1.0);
	h_["detectorMHToverMET"] = new TH1F("detectorMHToverMET", ";detector H_{T}^{miss} / E_{T}^{miss}; a.u.", 100, 0, 3.0);
	h_["detectorLeadingAk4JetPt"] = new TH1F("detectorLeadingAk4JetPt", ";AK4 Jet p_{T} (GeV);a.u.", 100, 0, 2500);
	h_["detectorSecondaryAk4JetPt"] = new TH1F("detectorSecondaryAk4JetPt", ";AK4 Jet p_{T} (GeV);a.u.", 100, 0, 2500);
	h_["detectorLeadingAk4JetEta"] = new TH1F("detectorLeadingAk4JetEta", ";#eta AK4 Jet;a.u.", 100, -5, 5);
	h_["detectorSecondaryAk4JetEta"] = new TH1F("detectorSecondaryAk4JetEta", ";#eta AK4 Jet;a.u.", 100, -5, 5);
	h2_["detectorLeadingAk4JetPt_detectorSecondaryAk4JetPt"] = new TH2F("detectorLeadingAk4JetPt_detectorSecondaryAk4JetPt", ";secondary AK4 Jet p_{T} (GeV);leading AK4 Jet p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);
	h2_["detectorHT_detectorSecondaryAk4JetPt"] = new TH2F("detectorHT_detectorSecondaryAk4JetPt", ";secondary AK4 Jet p_{T} (GeV);detector H_{T} (GeV)", 400, 0, 2500, 400, 0, 7000);
	h2_["detectorHT_detectorLeadingAk4JetPt"] = new TH2F("detectorHT_detectorLeadingAk4JetPt", ";leading AK4 Jet p_{T} (GeV);detector H_{T} (GeV)", 400, 0, 2500, 400, 0, 7000);
	h_["detectorNumAk4JetsOver10GeV"] = new TH1F("detectorNumAk4JetsOver10GeV", ";Number of AK4 Jets;a.u.", 25, 0, 25);
	h_["detectorNumAk4JetsOver20GeV"] = new TH1F("detectorNumAk4JetsOver20GeV", ";Number of AK4 Jets;a.u.", 25, 0, 25);
	h_["detectorNumAk4JetsOver40GeV"] = new TH1F("detectorNumAk4JetsOver40GeV", ";Number of AK4 Jets;a.u.", 25, 0, 25);
	h_["detectorNumAk4JetsOver100GeV"] = new TH1F("detectorNumAk4JetsOver100GeV", ";Number of AK4 Jets;a.u.", 25, 0, 25);
	h2_["detectorLogicAlphaTMinBiasedDeltaPhi"] = new TH2F("detectorLogicAlphaTMinBiasedDeltaPhi", ";Alpha T > 0.52;Biased Delta Phi > 0.50", 2, 0, 2, 2, 0, 2);
	h2_["detectorLogicMHToverMETmht"] = new TH2F("detectorLogicMHToverMETmht", ";H_{T}^{miss} / E_{T}^{miss} > 1.25;H_{T}^{miss} > 130 (GeV)", 2, 0, 2, 2, 0, 2);

	h2_["detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50andAlphaTmoreThan0p52"] = new TH2F("detector_MHTmoreThan130andMHToverMETlessThan1p25andBiasedDeltaPhiMoreThan0p50andAlphaTmoreThan0p52", ";pass?;", 2, 0, 2, 1, 0, 1);

	// fatJet Btag histograms
	h_["fatJetNumberLooseDoubleBTagsNoMatching"] = new TH1F("fatJetNumberLooseDoubleBTagsNoMatching", ";Number of Double B Tags;a.u.", 6, 0, 6);
	h_["fatJetNumberMediumDoubleBTagsNoMatching"] = new TH1F("fatJetNumberMediumDoubleBTagsNoMatching", ";Number of Double B Tags;a.u.", 6, 0, 6);
	h_["fatJetNumberTightDoubleBTagsNoMatching"] = new TH1F("fatJetNumberTightDoubleBTagsNoMatching", ";Number of Double B Tags;a.u.", 6, 0, 6);
	h_["fatJetNumberLooseDoubleBTagsWithMatching"] = new TH1F("fatJetNumberLooseDoubleBTagsWithMatching", ";Number of Double B Tags;a.u.", 5, 0, 5);
	h_["fatJetNumberMediumDoubleBTagsWithMatching"] = new TH1F("fatJetNumberMediumDoubleBTagsWithMatching", ";Number of Double B Tags;a.u.", 5, 0, 5);
	h_["fatJetNumberTightDoubleBTagsWithMatching"] = new TH1F("fatJetNumberTightDoubleBTagsWithMatching", ";Number of Double B Tags;a.u.", 5, 0, 5);
	h_["leadingMatchedFatJetPtTwoLooseTagsMatch"] = new TH1F("leadingMatchedFatJetPtTwoLooseTagsMatch", ";doubleBTagJet p_{T} (GeV); a.u.", 50, 0, 2500);
	h_["secondaryMatchedFatJetPtTwoLooseTagsMatch"] = new TH1F("secondaryMatchedFatJetPtTwoLooseTagsMatch", ";doubleBTagJet p_{T} (GeV); a.u.", 50, 0, 2500);
	h_["leadingMatchedFatJetPtTwoMediumTagsMatch"] = new TH1F("leadingMatchedFatJetPtTwoMediumTagsMatch", ";doubleBTagJet p_{T} (GeV); a.u.", 50, 0, 2500);
	h_["secondaryMatchedFatJetPtTwoMediumTagsMatch"] = new TH1F("secondaryMatchedFatJetPtTwoMediumTagsMatch", ";doubleBTagJet p_{T} (GeV); a.u.", 50, 0, 2500);
	h_["leadingMatchedFatJetPtTwoTightTagsMatch"] = new TH1F("leadingMatchedFatJetPtTwoTightTagsMatch", ";doubleBTagJet p_{T} (GeV); a.u.", 50, 0, 2500);
	h_["secondaryMatchedFatJetPtTwoTightTagsMatch"] = new TH1F("secondaryMatchedFatJetPtTwoTightTagsMatch", ";doubleBTagJet p_{T} (GeV); a.u.", 50, 0, 2500);
	h2_["leading_secondaryMatchedFatJetPtTwoLooseTagsMatch"] = new TH2F("leading_secondaryMatchedFatJetPtTwoLooseTagsMatch", ";secondary doubleBTagJet p_{T} (GeV); leading doubleBTagJet p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);
	h2_["leading_secondaryMatchedFatJetPtTwoMediumTagsMatch"] = new TH2F("leading_secondaryMatchedFatJetPtTwoMediumTagsMatch", ";secondary doubleBTagJet p_{T} (GeV); leading doubleBTagJet p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);
	h2_["leading_secondaryMatchedFatJetPtTwoTightTagsMatch"] = new TH2F("leading_secondaryMatchedFatJetPtTwoTightTagsMatch", ";secondary doubleBTagJet p_{T} (GeV); leading doubleBTagJet p_{T} (GeV)", 400, 0, 2500, 400, 0, 2500);

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






// function to get the nmssm cascade particles, returns true if does so sucessfully:)
bool getCascadeParticles(edm::Handle<std::vector<reco::GenParticle>> genParticles, int ievt, std::string filename, unsigned int & gluinoCount, std::vector<reco::GenParticle> & squarkVec, std::vector<const reco::Candidate*> & qjetVec, std::vector<const reco::Candidate*> & nlspVec, std::vector<const reco::Candidate*> & lspVec, std::vector<const reco::Candidate*> & higgsVec, std::vector<const reco::Candidate*> & bVec, std::vector<const reco::Candidate*> & bbarVec)
{
	for (size_t iGen=0; iGen<genParticles->size(); ++iGen){
		const reco::GenParticle & genParticle = (*genParticles)[iGen];
		
		// count the gluinos
		if ( abs(genParticle.pdgId()) == 1000021 && genParticle.numberOfDaughters() == 2){
			
			int pdgId_0 = abs(genParticle.daughter(0)->pdgId());
			int pdgId_1 = abs(genParticle.daughter(1)->pdgId());				
			if (   (pdgId_0 == 1000001 || pdgId_0 == 1000002 || pdgId_0 == 1000003 || pdgId_0 == 1000004 || pdgId_0 == 2000001 || pdgId_0 == 2000002 || pdgId_0 == 2000003 || pdgId_0 == 2000004)
			    || (pdgId_1 == 1000001 || pdgId_1 == 1000002 || pdgId_1 == 1000003 || pdgId_1 == 1000004 || pdgId_1 == 2000001 || pdgId_1 == 2000002 || pdgId_1 == 2000003 || pdgId_1 == 2000004))
				gluinoCount++;
		}
		
		// get the squarks
		if ( (abs(genParticle.pdgId()) == 1000001
		     || abs(genParticle.pdgId()) == 2000001
		     || abs(genParticle.pdgId()) == 1000002
		     || abs(genParticle.pdgId()) == 2000002
		     || abs(genParticle.pdgId()) == 1000003
		     || abs(genParticle.pdgId()) == 2000003
		     || abs(genParticle.pdgId()) == 1000004
		     || abs(genParticle.pdgId()) == 2000004) 
			 && genParticle.numberOfDaughters() == 2){

			int pdgId_0 = abs(genParticle.daughter(0)->pdgId());
			int pdgId_1 = abs(genParticle.daughter(1)->pdgId());

			if ( (pdgId_1 == 1 || pdgId_1 == 2 || pdgId_1 == 3 || pdgId_1 ==4) && pdgId_0 == 1000023 ){
				squarkVec.push_back(genParticle);
				qjetVec.push_back(genParticle.daughter(1));
				nlspVec.push_back(genParticle.daughter(0));		
			} 
			else if ( (pdgId_0 == 1 || pdgId_0 == 2 || pdgId_0 == 3 || pdgId_0 ==4) && pdgId_1 == 1000023 ){
				squarkVec.push_back(genParticle);					
				qjetVec.push_back(genParticle.daughter(0));
				nlspVec.push_back(genParticle.daughter(1));
			}
			else continue; // go look at the next candidate genParticle

			// nlsp decays
			// loop to ensure we get the decaying version of the nlsp
			while (nlspVec.back()->numberOfDaughters() != 0){
				if (nlspVec.back()->numberOfDaughters() == 2) break;
				else nlspVec.back() = nlspVec.back()->daughter(0);
			}

			pdgId_0 = abs(nlspVec.back()->daughter(0)->pdgId());
			pdgId_1 = abs(nlspVec.back()->daughter(1)->pdgId());
			if ( pdgId_1 == 35 && pdgId_0 == 1000022 ){
				higgsVec.push_back(nlspVec.back()->daughter(1));
				lspVec.push_back(nlspVec.back()->daughter(0));			
			} 
			else if ( pdgId_0 == 35 && pdgId_1 == 1000022 ){
				higgsVec.push_back(nlspVec.back()->daughter(0));
				lspVec.push_back(nlspVec.back()->daughter(1));
			}
			else {
				std::cerr << "ERROR: nlsp does not have the correct daughters" << std::endl;
				std::cerr << "File: " << filename << std::endl;
				std::cerr << "Event: " << ievt << std::endl;
				std::cerr << std::endl;
				return false;
			}

			// higgs decays
			if (higgsVec.back()->numberOfDaughters() != 2 ){
				std::cerr << "ERROR: higgs does not have two daughters" << std::endl;
				std::cerr << "File: " << filename << std::endl;
				std::cerr << "Event: " << ievt << std::endl;
				std::cerr << std::endl;
				return false;
			}
			pdgId_0 = higgsVec.back()->daughter(0)->pdgId();
			pdgId_1 = higgsVec.back()->daughter(1)->pdgId();
			if ( pdgId_0 == 5 && pdgId_1 == -5 ){
				bVec.push_back(higgsVec.back()->daughter(0));
				bbarVec.push_back(higgsVec.back()->daughter(1));
			}
			else if ( pdgId_1 == 5 && pdgId_0 == -5 ){
				bVec.push_back(higgsVec.back()->daughter(1));
				bbarVec.push_back(higgsVec.back()->daughter(0));			
			} 
			else {
				std::cerr << "ERROR: higgs does not have the correct daughters" << std::endl;
				std::cerr << "File: " << filename << std::endl;
				std::cerr << "Event: " << ievt << std::endl;
				std::cerr << std::endl;
				return false;
			}

		} // closes 'if' genParticle is a squark with two daughters
	} // closes loop through genParticles

	// check we have the right number of correctly decaying squarks
	if (squarkVec.size() != 2){
		std::cerr << "ERROR: haven't got two correctly decaying squarks" << std::endl;
		std::cerr << "File: " << filename << std::endl;
		std::cerr << "Event: " << ievt << std::endl;
		std::cerr << std::endl;
		return false;	
	}

	// ensure squarkVec entries are in pt order, and that decay arms follow this rule also
	if (squarkVec[1].pt() > squarkVec[0].pt()){
		std::swap(squarkVec[0],squarkVec[1]);
		std::swap(qjetVec[0],qjetVec[1]);
		std::swap(nlspVec[0],nlspVec[1]);
		std::swap(lspVec[0],lspVec[1]);
		std::swap(higgsVec[0],higgsVec[1]);
		std::swap(bVec[0],bVec[1]);
		std::swap(bbarVec[0],bbarVec[1]);
	}

	return true;
} // closes the funciton getCascadeParticles