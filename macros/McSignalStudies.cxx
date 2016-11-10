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
#include <TBranch.h>
#include <TClonesArray.h>

// DELPHES headers
#include "Analysis/Analysis_boostedNmssmHiggs/interface/DelphesClasses.h"
#include "Analysis/Analysis_boostedNmssmHiggs/interface/ExRootTreeReader.h"

// Headers from this package
// #include "Analysis/Analysis_boostedNmssmHiggs/interface/Kinematics.h"
// #include "Analysis/Analysis_boostedNmssmHiggs/interface/PlottingDoubleBTaggerEfficiencyStudies.h"


void CreateHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&);
void WriteHistograms(std::map<std::string,TH1F*>&, std::map<std::string,TH2F*>&, std::string);

void McSignalStudies() 
{
	// already loaded in root session
    // gSystem->Load("libTreePlayer");
    // gSystem->Load("$HOME/MG5_aMC_v2_3_3/Delphes/libDelphes.so");

	// Create histograms, they are accessed by eg: h_["fatJetMass_loose"]->Fill(125.0);
	std::map<std::string, TH1F*> h_;
	std::map<std::string, TH2F*> h2_;
	CreateHistograms(h_, h2_);

	// Running Options
	int maxEvents_ = 50; // -1 for all events
	unsigned int outputEvery_ = 1;
	std::vector<std::string> inputFiles_ = {"~/tag_1_delphes_events.root"};
	std::string outputFile_ = "ABC_v3/output.root";
	bool justDoPlotting_ = false;

	std::string outputDirectory_;
	std::string forwardSlash = "/";
	// strip the directory from the outputfile name
	for (size_t c = outputFile_.size()-1; c > 0; --c){
		if (outputFile_[c] == forwardSlash[0]){
			outputDirectory_ = outputFile_.substr(0, c+1);
			break;
		}
	}
	bool makeDir = !(std::system(Form("mkdir %s",outputDirectory_.c_str())));
	if (justDoPlotting_ == false && makeDir == false){
		std::cout << "The chosen output directory already exists, or the parent directory does not exist:" << std::endl;
		std::cout << "Do not wish to overwrite ROOT file: Exiting..." << std::endl;
		return 1;
	}

	// copy the code used to make the histogram ROOT file into the same directory (this could be out of sync if you edit after compilation)
	if (justDoPlotting_ == false) std::system(Form("cp $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/macros/McSignalStudies.cxx %s",outputDirectory_.c_str()));		


	///////////////////////////////
	// Loop through the input files
	///////////////////////////////
	int ievt=0;
	// but if we are just doing the plotting goto the end
	if (justDoPlotting_) goto plottingLabel;

	for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
	    
	    // Open input file
	    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
	    if( inFile ){
		
			TTree * t = (TTree*)inFile->Get("Delphes");
			ExRootTreeReader *reader = new ExRootTreeReader(t);
			TClonesArray *branchParticle = reader->UseBranch("Particle");
			TClonesArray *branchMET = reader->UseBranch("MissingET");

		    Int_t nevents = (Int_t)t->GetEntries();
		    for (Int_t ievt=0; ievt<nevents; ievt++){  

				// break loop if maximal number of events is reached 
				if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
					
				// event counter
				if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false){
					std::cout << "   File " << iFile+1 << " of " << inputFiles_.size() << ":";
					std::cout << "  processing event: " << ievt << std::endl;
				}

				// load info for this event
		        reader->ReadEntry(ievt);

				// load the gen particle info
				std::vector<GenParticle*> particleVec;
				for (int iPar=0; iPar<branchParticle->GetEntries(); ++iPar){
				    GenParticle * particle = (GenParticle*) branchParticle->At(iPar);
				    particleVec.push_back(particle);
				    // std::cout << particleVec[iPar]->PT << std::endl; // an example of how to read info from this object we created
				}

				// load the met info
				std::vector<MissingET*> metVec;
				for (int iMet=0; iMet<branchMET->GetEntries(); ++iMet){
				    MissingET * met = (MissingET*) branchMET->At(iMet);
				    metVec.push_back(met);
				}


				// A N A L Y S I S
				// ------------------------------------------------------------------------------------------------------------//
				// ------------------------------------------------------------------------------------------------------------//
				
				// count how many gluinos are involved, we know we expect 2 squarks
				unsigned int gluinoCount = 0;

				// the first element is leading arm
				// the second element is for the secondary arm
				// access leading lsp PT like so: particleVec[lspIndices[0]]->PT
				// access secondary lsp PT like so: particleVec[lspIndices[1]]->PT
				std::vector<int> squarkIndices;
				std::vector<int> qjetIndices;
				std::vector<int> nlspIndices;
				std::vector<int> lspIndices;
				std::vector<int> higgsIndices;
				std::vector<int> bIndices;
				std::vector<int> bbarIndices;

				// *** SQUARKS ***
				// loop through gen particle entries to find the squarks (and count the gluinos)
				for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
					if ( abs(particleVec[iPar]->PID) == 1000021 ) gluinoCount++;
					if ( abs(particleVec[iPar]->PID) == 1000001
					|| abs(particleVec[iPar]->PID) == 2000001
					|| abs(particleVec[iPar]->PID) == 1000002
					|| abs(particleVec[iPar]->PID) == 2000002
					|| abs(particleVec[iPar]->PID) == 1000003
					|| abs(particleVec[iPar]->PID) == 2000003
					|| abs(particleVec[iPar]->PID) == 1000004
					|| abs(particleVec[iPar]->PID) == 2000004) squarkIndices.push_back(iPar);
				} // closes loop through gen particle vector
				if (squarkIndices.size() != 2){
					std::cout << "WARNING" << std::endl;
					// could elaborate and write to text or something the problem...
					continue;
				} // closes 'if' we don't have the 2 squarks
				if (particleVec[squarkIndices[1]]->PT > particleVec[squarkIndices[0]]->PT) std::swap(squarkIndices[0],squarkIndices[1]);


				// NB: the squark object do not have daughters...so we have to find cascade particles in a backwards kind of way			

				// *** QJET ***
				// loop through gen particle entries to find the quarks (squark->QUARK+nlsp)
				// nb: currently these are set in madgraph to only be u or d quarks
				// leading arm quark
				for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
					if ( abs(particleVec[iPar]->PID) == 1 || abs(particleVec[iPar]->PID) == 2 ){
						if (particleVec[iPar]->M1 == squarkIndices[0]){
							qjetIndices.push_back(iPar);  
							break;
						}
					}
				} // closes loop through gen particle vector
				// secondary arm quark
				for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
					if ( abs(particleVec[iPar]->PID) == 1 || abs(particleVec[iPar]->PID) == 2 ){
						if (particleVec[iPar]->M1 == squarkIndices[1]){
							qjetIndices.push_back(iPar);  
							break;
						}
					}
				} // closes loop through gen particle vector
				// check there are two entries


				// *** NLSP ***
				// loop through gen particle entries to find the nlsp's (squark->quark+NLSP)
				// leading arm quark
				for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
					if ( abs(particleVec[iPar]->PID) == 1000023 ){
						if (particleVec[iPar]->M1 == squarkIndices[0]){
							nlspIndices.push_back(iPar);  
							break;
						}
					}
				} // closes loop through gen particle vector
				// secondary arm quark
				for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
					if ( abs(particleVec[iPar]->PID) == 1000023 ){
						if (particleVec[iPar]->M1 == squarkIndices[1]){
							nlspIndices.push_back(iPar);  
							break;
						}
					}
				} // closes loop through gen particle vector
				// check there are two entries



				// *** LSP ***
				// loop through gen particle entries to find the lsp's (nlsp->LSP+higgs)
				// leading arm quark
				for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
					if ( abs(particleVec[iPar]->PID) == 1000022 ){
						if (particleVec[iPar]->M1 == nlspIndices[0]){
							lspIndices.push_back(iPar);  
							break;
						}
					}
				} // closes loop through gen particle vector
				// secondary arm quark
				for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
					if ( abs(particleVec[iPar]->PID) == 1000022 ){
						if (particleVec[iPar]->M1 == nlspIndices[1]){
							lspIndices.push_back(iPar);  
							break;
						}
					}
				} // closes loop through gen particle vector
				// check there are two entries
				


				// *** HIGGS ***
				// loop through gen particle entries to find the lsp's (nlsp->LSP+higgs)
				// leading arm quark
				for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
					if ( abs(particleVec[iPar]->PID) == 35 ){
						if (particleVec[iPar]->M1 == nlspIndices[0]){
							higgsIndices.push_back(iPar);  
							break;
						}
					}
				} // closes loop through gen particle vector
				// secondary arm quark
				for (size_t iPar=0; iPar<particleVec.size(); ++iPar){
					if ( abs(particleVec[iPar]->PID) == 35 ){
						if (particleVec[iPar]->M1 == nlspIndices[1]){
							higgsIndices.push_back(iPar);  
							break;
						}
					}
				} // closes loop through gen particle vector
				// check there are two entries



				// do the same for the b's

				// do the same for the bbars's






				// End of Analysis
				// ------------------------------------------------------------------------------------------------------------//
				// ------------------------------------------------------------------------------------------------------------//

			} // closes loop through events for this file
			inFile->Close();
	    } // closes 'if' the file exists

	    // break loop if maximal number of events is reached:
	    // this has to be done twice to stop the file loop as well
	    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;

	} // closes loop through files

	WriteHistograms(h_, h2_, outputFile_.c_str());

	plottingLabel: //DO SOMETHING SIMILAR HERE!!!
	// if (std::system(Form("test -e %s",outputFile_.c_str())) == 0) // check the file exists
	// 	PlottingDoubleBTaggerEfficiencyStudies createPlots(outputFile_.c_str(), doubleBtagWPname, etaBinning, massCut);
	// else{
	// 	std::cout << "the following output root file does not exist to make plots from:" << std::endl;
	// 	std::cout << outputFile_ << std::endl;
	// }

return 0;
} // closes the 'main' function






void CreateHistograms(std::map<std::string,TH1F*> & h_, std::map<std::string,TH2F*> & h2_)
{
	// set the binning for histograms
    std::vector<double> ptBinning;
    for(double binLowerEdge=  0.0; binLowerEdge< 100.0; binLowerEdge+= 100.0) ptBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  100.0; binLowerEdge< 200.0; binLowerEdge+= 50.0) ptBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  200.0; binLowerEdge< 350.0; binLowerEdge+= 25.0) ptBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  350.0; binLowerEdge< 500.0; binLowerEdge+= 50.0) ptBinning.push_back(binLowerEdge);    	
    for(double binLowerEdge=  500.0; binLowerEdge< 600.0; binLowerEdge+= 100.0) ptBinning.push_back(binLowerEdge);
    for(double binLowerEdge=  600.0; binLowerEdge< 800.1; binLowerEdge+= 200.0) ptBinning.push_back(binLowerEdge);

	// create the debugging histograms
	h_["DEBUG_higgsBbDRpreMatching"] = new TH1F("DEBUG_higgsBbDRpreMatching", ";dR_bb;a.u.", 50, 0, 2.50);
	h_["DEBUG_higgsBbDRpreMatching_pt400to415_eta0to2p4"] = new TH1F("DEBUG_higgsBbDRpreMatching_pt400to415_eta0to2p4", ";dR_bb;a.u.", 100, 0, 2.50);
	h_["DEBUG_higgsBbDRpreMatching_pt450to465_eta0to2p4"] = new TH1F("DEBUG_higgsBbDRpreMatching_pt450to465_eta0to2p4", ";dR_bb;a.u.", 100, 0, 2.50);
	h_["DEBUG_higgsBbDRpreMatching_pt500to515_eta0to2p4"] = new TH1F("DEBUG_higgsBbDRpreMatching_pt500to515_eta0to2p4", ";dR_bb;a.u.", 100, 0, 2.50);
	h_["DEBUG_higgsBbDRpreMatching_3p400to415_eta0to2p4"] = new TH1F("DEBUG_higgsBbDRpreMatching_3p400to415_eta0to2p4", ";dR_bb;a.u.", 100, 0, 2.50);
	h_["DEBUG_higgsBbDRpreMatching_3p450to465_eta0to2p4"] = new TH1F("DEBUG_higgsBbDRpreMatching_3p450to465_eta0to2p4", ";dR_bb;a.u.", 100, 0, 2.50);
	h_["DEBUG_higgsBbDRpreMatching_3p500to515_eta0to2p4"] = new TH1F("DEBUG_higgsBbDRpreMatching_3p500to515_eta0to2p4", ";dR_bb;a.u.", 100, 0, 2.50);
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