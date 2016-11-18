// CPP headers
#include <string>
#include <vector>
#include <iostream>

// ROOT headers
#include <TROOT.h>

// script to submit McSignalStudies many times
// run with
// $ root -q -b -l submit_McSignalStudies.cxx
// change the input and output vectors as you wish:)

void submit_McSignalStudies(){

	std::vector<std::vector<std::string>> vecInputFiles_; // can add more than one file if you do so wish
	// vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/testing/mH125p0_mSusy1000p0_ratio0p96_splitting2p0_10000events/dirA/dirB/dirC/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH70p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH70p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH70p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH70p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH70p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH70p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH80p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH80p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH80p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH80p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH80p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH80p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH90p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH90p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH90p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH90p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH90p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH90p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH100p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH100p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH100p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH100p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH100p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH100p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH110p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH110p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH110p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH110p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH110p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH110p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH120p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH120p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH120p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH120p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH120p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH120p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH130p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH130p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH130p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH130p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH130p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});
	vecInputFiles_.push_back({"/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/mH130p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/Events/run_01/tag_1_delphes_events.root"});

	std::vector<std::string> vecOutputFile_; // follow the mHXXpX_mSusyXXXXpX_ratioXpXX_splittingXpX_XXXXXevents/output.root format please:)
	// vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/test2/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH70p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH70p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH70p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH70p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH70p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH70p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH80p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH80p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH80p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH80p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH80p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH80p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH90p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH90p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH90p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH90p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH90p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH90p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH100p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH100p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH100p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH100p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH100p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH100p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH110p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH110p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH110p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH110p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH110p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH110p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH120p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH120p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH120p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH120p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH120p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH120p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH130p0_mSusy1000p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH130p0_mSusy1200p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH130p0_mSusy1400p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH130p0_mSusy1600p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH130p0_mSusy1800p0_ratio0p95_splitting1p0_10000events/output.root");
	vecOutputFile_.push_back("/users/jt15104/local_Analysis_boostedNmssmHiggs/output_McSignalStudies/nmssmCascadeAnalysis_v01/mH130p0_mSusy2000p0_ratio0p95_splitting1p0_10000events/output.root");


	if (vecInputFiles_.size() != vecOutputFile_.size()){
		std::cout << "The size of the two inputs vectors is not the same. Exiting..." << std::endl;
		return;
	}


	for (size_t i=0; i<vecOutputFile_.size(); ++i){
		std::cout << "******************************" << std::endl;
		std::cout << "*** Doing MG setup " << i+1 << " of " << vecOutputFile_.size() << " ***" << std::endl;
		std::cout << "*** Project Name: " << vecInputFiles_[i][0] << " ***" << std::endl;
		std::cout << "******************************" << std::endl;

		std::string submissionString = "{";
		for (size_t j=0; j<vecInputFiles_[i].size(); ++j){
			submissionString = submissionString + "\"";
			submissionString = submissionString + vecInputFiles_[i][j];
			submissionString = submissionString + "\",";
		}
		submissionString = submissionString + "}, \"" + vecOutputFile_[i] + "\"";
		// std::cout << submissionString << std::endl;

		std::system(Form("root -q -l -b 'McSignalStudies.cxx(%s)'", submissionString.c_str()));

	}

} // closes the function "McSignalStudies"
