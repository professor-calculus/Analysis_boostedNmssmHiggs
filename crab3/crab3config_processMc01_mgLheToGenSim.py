from CRABClient.UserUtilities import config
import os
config = config()

# submit, from crab directory, with $ crab submit -c crab3config_MC.py
# but first...

# 1. ensure that you have unzipped the necessary madGraph files using python/unzipMadgraphLhe.py

# 2. EDIT THE FOLLOWING
madGraphProjectName = 'mH110p0_mSusy1800p0_ratio0p95_splitting1p0_25000events'
uniqueName = 'processMc01_nmssmSignalV01_REPOcommitHash'

# 3. MAKE SURE THAT THE FOLLOWING ARE UP TO DATE
madGraphProjectStore = '/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/paramCard_type02/'
pathWithinMadgraphProject = "Events/run_01/events.lhe"

# 4. MAKE SURE unitsPerJob (and totalUnits) ARE SET CORRECTLY
config.Data.unitsPerJob = 1 # number of events per job
config.Data.totalUnits = 2 # -1 for all

config.Data.splitting = 'EventBased' 
config.Data.inputDBS = 'global'
config.Data.publication = False

uniqueNameMGP = uniqueName + "_" + madGraphProjectName # combine info for a truly unique record
config.General.requestName = uniqueNameMGP # name of the crab job project (eg on dashboard)
config.Data.outputDatasetTag = uniqueNameMGP # name for storage directory beneath /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/CRAB_PrivateMC/

config.General.workArea = 'crab_projects' # name of directory where crab project info is stored
config.General.transferOutputs = True
config.General.transferLogs    = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'processMc01_mgLheToGenSim_cfg.py'
config.JobType.inputFiles = [os.path.join(madGraphProjectStore,madGraphProjectName,pathWithinMadgraphProject)]

config.Site.storageSite = 'T2_UK_SGrid_Bristol'


### ## # EXTRA INFO # ## ###
#
# this crab job turn MG_LHE files into CMSSW_GENSIM ROOT files
# part 1/3 of making CMSSW_RECO data from MG output
#
# output ROOT files will go here:
# /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/CRAB_PrivateMC/<uniqueNameMGP>/<dateStamp>/0000/
#
# the names of the output ROOT files will be:
# nmssmSignalV01_GENSIM_<jobNum>.root
# as stated by the config file processMc01_mgLheToGenSim_cfg.py
############################