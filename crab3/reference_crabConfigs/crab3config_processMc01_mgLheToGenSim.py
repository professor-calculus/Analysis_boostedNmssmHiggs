from CRABClient.UserUtilities import config
import os
config = config()

# NOTE THAT THIS NOW EXISTS PRIMARILY AS A RECORD
# SHOULD BE USING mcProduction_Controller.py FOR THIS WORKFLOW

# submit, from crab directory, with $ crab submit -c crab3config_processMc01_mgLheToGenSim.py
# but first...

# 1. ensure that you have unzipped the necessary madGraph files using python/unzipMadgraphLhe.py

# 2. EDIT THE FOLLOWING
madGraphProjectName = 'mH110p0_mSusy1800p0_ratio0p95_splitting1p0_25000events'
partOneUniqueName = 'processMc01_nmssmSignalV01_TEST10'

# 3. MAKE SURE THAT THE FOLLOWING ARE UP TO DATE
madGraphProjectStore = '/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/paramCard_type02/'
pathWithinMadgraphProject = "Events/run_01/unweighted_events.lhe"

# 4. MAKE SURE unitsPerJob (and totalUnits) ARE SET CORRECTLY
config.Data.unitsPerJob = 5 # number of events per job
config.Data.totalUnits = 100 # -1 for all

config.Data.splitting = 'EventBased' 
config.Data.inputDBS = 'global'
config.Data.publication = True # means we can access the data in step 2
config.Data.outputPrimaryDataset = 'nmssmSignalCascadeV01_13TeV' 

partOneUniqueNameMGP = partOneUniqueName + "_" + madGraphProjectName # combine info for a truly unique record (must have less than 100 characters)
config.General.requestName = partOneUniqueNameMGP # name of the crab job project (eg on dashboard)
config.Data.outputDatasetTag = partOneUniqueNameMGP # name for storage directory beneath /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/CRAB_PrivateMC/

config.General.workArea = 'crab_projects' # name of directory where crab project info is stored
config.General.transferOutputs = True
config.General.transferLogs    = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'processMc01_mgLheToGenSim_cfg.py' # the CMSSW config file used
config.JobType.inputFiles = [os.path.join(madGraphProjectStore,madGraphProjectName,pathWithinMadgraphProject)] # copies this file over to the workernode

config.Site.storageSite = 'T2_UK_SGrid_Bristol'


### ## # EXTRA INFO # ## ###
#
# this crab job turns MG_LHE files into CMSSW_GENSIM ROOT files
# part 1/3 of making CMSSW_RECO data from MG output
#
# the public location of the ouput will go here:
# /<outputPrimaryDataset>/taylor-<partOneUniqueNameMGP>-<randomHash>/USER
#
# output ROOT files will go here:
# /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/<outputPrimaryDataset>/<partOneUniqueNameMGP>/<dateStamp>/0000/
#
# the names of the output ROOT files will be:
# nmssmSignal_GENSIM_<jobNum>.root
# as stated by the config file processMc01_mgLheToGenSim_cfg.py
############################