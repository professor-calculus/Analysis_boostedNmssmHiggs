Information on how to create/process data.

You will first need a valid grid certificate
$ voms-proxy-init -voms cms --valid 168:00

You will also need to source crab3
$ source /cvmfs/cms.cern.ch/crab3/crab.sh (shortcut $loadCrab)
If you are using an old version of CMSSW you will want to use the following
$ source /cvmfs/cms.cern.ch/crab3/crab_standalone.sh (shortcut $loadCrabOldSkuul)

Crab jobs will go to '/hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/@dataset/@outputDatasetTag/@dateStamp/0000/'
You should copy them to a suitable place on '/hdfs/user/jt15104/Analysis_boostedNmssmHiggs/' using hadoop commands

All the python configuration files direct you to what needs changing on a run by run basis.



ONE. Processing AOD MC into patTuples with fat double B tags.
Note that you may have to use a specific version of CMSSW to get this to work (eg CMSSW_8_0_20)
$ crab submit -c crab3config_patTuple_addBTagging.py
This calls on the CMSSW config file patTuple_addBTagging_cfg.py



TWO. Create signal MC.
This is composed of three stages, starting from the point where we have signal MadEvent LHE files
Note that you may have to use a specific version of CMSSW to get this to work (eg CMSSW_8_0_3_patch1)
(The stages are copied from the commands used by central production for recent susy MC. (Details in mcProductionWorkFlow.txt)

Old
1. Put LHE through pythia8 and simulates the CMS detector
$ crab submit -c crab3config_processMc01_mgLheToGenSim.py

2. First step going from GenSim to AOD (with genParticles)
$ crab submit -c crab3config_processMc02_genSimToAOD_step1of2.py

3. Second step going from GenSim to AOD (with genParticles)
$ crab submit -c crab3config_processMc03_genSimToAOD_step2of2.py

New
This is now all controlled my McProduction_Controller.py
it easily does submission, resubmission, and status checks
and passes the jobs through the necessary stages easily.
Select the correct settings then...
$ python McProduction_Controller.py (from the crab3 directory)
The one main pain in the bum is that the main MC stages use a different CMSSW version to the patTuple stage
So there are options related to this you must check are correct
(don't worry, the code won't let you submit from the wrong version)
