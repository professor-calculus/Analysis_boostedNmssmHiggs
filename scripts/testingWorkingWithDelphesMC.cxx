#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TPad.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <string>

#include "TChain.h"


// #include "TSystem.h"
// #ifdef __CLING__
// R__LOAD_LIBRARY(libDelphes)
// #include "classes/DelphesClasses.h"
// #include "external/ExRootAnalysis/ExRootTreeReader.h"
// #include "external/ExRootAnalysis/ExRootResult.h"
// #endif



#include "TSystem.h"




void testingWorkingWithDelphesMC(){

// gSystem->Load("libDelphes");
  gSystem->Load("libTreePlayer");
    gSystem->Load("$HOME/MG5_aMC_v2_3_3/Delphes/libDelphes.so");

// create histos function

  // fill histos...maybe not

      // write histos...function


  // make trees
  std::cout << "Loading up the TChains..." << std::endl;
  TChain * t = new TChain("Delphes");
 
  std::string inputFile01 = "~/tag_1_delphes_events.root";

// use a command line parser...
    t->Add(inputFile01.c_str());



  //DECLARE THE VARIABLES WE WILL USE
  Int_t nevents = (Int_t)t->GetEntries(); //total number of tree entries

  Float_t MET = 0;
  t->SetBranchAddress("MissingET.MET", &MET);




// TBranch *b_met = t->GetBranch("MissingET.MET");
// b_met->SetAddress(&MET);

//t1->SetBranchAddress("ETmis_size", &met);
      // Long64_t nentries = t1->GetEntries();






  // double particle_q = 0;
  // t->SetBranchAddress("Particle.Charge" &particle_q);

  // //neutrino detection variable
  // Float_t MET = 0;                    //missing transverse energy


  // Int_t mu_candidateNumber = 0;        //number of muon candidates in an entry
  // vector<double> * mu_pt = 0;          //array of momenta for each candidate in an entry

  // Int_t mc_candidateNumber = 0;     //no. of monte carlo particles in an event...
  // vector<double> * mc_pt = 0;       //AND the TRUE values. (mc_ subscript)
  // vector<double> * mc_phi = 0;
  // vector<double> * mc_ID = 0;        //identifies each monte carlo particle
  // vector<double> * mc_status = 0;    //array of values representing a certain process

  // t->SetBranchAddress("MET_RefFinal_et", &MET);
  // t->SetBranchAddress("mu_staco_n", &mu_candidateNumber);
  // t->SetBranchAddress("mu_staco_pt", &mu_pt);
  // t->SetBranchAddress("mc_n", &mc_candidateNumber);
  // t->SetBranchAddress("mc_pt", &mc_pt);
  // t->SetBranchAddress("mc_phi", &mc_phi);
  // t->SetBranchAddress("mc_pdgId", &mc_ID);
  // t->SetBranchAddress("mc_status", &mc_status);


///////////////////////
//LOOP THROUGH EVENTS//
///////////////////////
     for (Int_t i=0; i<nevents; i++){       //this is for full length runs
      
        t->GetEntry(i);               //loads all data for a given event


        std::cout << MET << std::endl;        


        // //TRUTH
        // //finds the true Z-->nunu neutrinos
        // mk=0;
        // for (m=0; m<mc_candidateNumber && mk<2; m++){

        //   if ( mc_status->at(m) == 3  &&  (  abs(mc_ID->at(m))==12 ||  abs(mc_ID->at(m))==14  ||   abs(mc_ID->at(m))==16 ) ){
        //     //ID 12,14,16 are the neutrinos
        //     mk++;
        //     if (mk==1) {m1=m;}    //sets the labels to the neutrinos
        //     if (mk==2) {m2=m;}
        //   }
        // }
        // //Calculations&Plots with the true nunu pair
        // mc_xptZ = mc_pt->at(m1) * cos(mc_phi->at(m1)) + mc_pt->at(m2) * cos(mc_phi->at(m2));
        // mc_yptZ = mc_pt->at(m1) * sin(mc_phi->at(m1)) + mc_pt->at(m2) * sin(mc_phi->at(m2));
        // mc_ptZ = sqrt( pow(mc_xptZ,2) + pow(mc_yptZ,2) );


      if (i % 1000 == 0){    //to see how it is progressing
      cout << i << " out of " << nevents << "\n";}     
     }//close loop through events (for a given file)

    
    
 
   // //create a root file to save all the histograms to. (had to use a different file pointer)
   //  TFile g( outputFilename.c_str() , "new");

   //  //save the histograms to the stated file
   //  h2->Write();



}//closes function

