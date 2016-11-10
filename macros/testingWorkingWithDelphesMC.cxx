#include "TROOT.h"
#include "TFile.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TSystem.h"
#include "TClonesArray.h"

#include <iostream>
#include <math.h>
#include <vector>
#include <string>

#include "/users/jt15104/MG5_aMC_v2_3_3/Delphes/classes/DelphesClasses.h"
#include "/users/jt15104/MG5_aMC_v2_3_3/Delphes/external/ExRootAnalysis/ExRootTreeReader.h"
#include "/users/jt15104/MG5_aMC_v2_3_3/Delphes/external/ExRootAnalysis/ExRootResult.h"

// #include  "/afs/cern.ch/user/t/taylor/MG5_aMC_v2_3_3/Delphes/classes/DelphesClasses.h"
// #include  "/afs/cern.ch/user/t/taylor/MG5_aMC_v2_3_3/Delphes/external/ExRootAnalysis/ExRootTreeReader.h"
// #include  "/afs/cern.ch/user/t/taylor/MG5_aMC_v2_3_3/Delphes/external/ExRootAnalysis/ExRootResult.h"


void testingWorkingWithDelphesMC(){

    gSystem->Load("libTreePlayer");
    gSystem->Load("$HOME/MG5_aMC_v2_3_3/Delphes/libDelphes.so");

// create histos function

  // fill histos...maybe not

      // write histos...function


    // // create the TChain
    // std::cout << "Loading up the TChains..." << std::endl;
    // TChain chain("Delphes");
    // chain.Add("~/tag_1_delphes_events.root");


    TFile * f = new TFile("~/tag_1_delphes_events.root");     //9
      TTree * t = (TTree*)f->Get("Delphes");
  //DECLARE THE VARIABLES WE WILL USE
  Int_t totalEntries = (Int_t)t->GetEntries(); //total number of tree entries


    // setup branches to read from
    ExRootTreeReader *reader = new ExRootTreeReader(t);
    TClonesArray *branchParticle = reader->UseBranch("Particle");
    TClonesArray *branchMET = reader->UseBranch("MissingET");

    // loop through events
    // Int_t nevents = (Int_t)chain.GetEntries();
    Int_t nevents = 10;
    for (Int_t i=0; i<nevents; i++){
      
        reader->ReadEntry(i);

        // load the gen particle info
        std::vector<GenParticle*> particleVec;
        for (size_t iPar=0; iPar<branchParticle->GetEntries(); ++iPar){
            GenParticle * particle = (GenParticle*) branchParticle->At(iPar);
            particleVec.push_back(particle);
            std::cout << particleVec[iPar]->PT << std::endl; // an example of how to read info from this object we created
        }

        // load the met info
        std::vector<MissingET*> metVec;
        for (size_t iMet=0; iMet<branchMET->GetEntries(); ++iMet){
            MissingET * met = (MissingET*) branchMET->At(iMet);
            metVec.push_back(met);
            // std::cout << metVec[iMet]->PT << std::endl;
        }












        if (i % 1000 == 0) std::cout << i << " out of " << nevents << std::endl;
    } // closes loop through events

    
    
 // write the histograms



}//closes function

