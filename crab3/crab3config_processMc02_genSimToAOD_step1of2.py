from CRABClient.UserUtilities import config
import os
config = config()

# NOTE THAT THIS NOW EXISTS PRIMARILY AS A RECORD
# SHOULD BE USING mcProduction_Controller.py FOR THIS WORKFLOW

# submit, from crab directory, with $ crab submit -c crab3config_processMc02_genSimToAOD_step1of2.py
# but first...

# 1. EDIT THE FOLLOWING #
madGraphProjectName = 'mH110p0_mSusy1800p0_ratio0p95_splitting1p0_25000events'
partTwoUniqueName = 'processMc02_nmssmSignalV01_F5'
# partOneUniqueName = 'processMc01_nmssmSignalV01_REPOcommitHash' # for the input dataset MIGHT BE MORE COMPLEX!!!


partTwoUniqueNameMGP = partTwoUniqueName + "_" + madGraphProjectName # combine info for a truly unique record
# partOneUniqueNameMGP = partOneUniqueName + "_" + madGraphProjectName

config.General.requestName = partTwoUniqueNameMGP # name of the crab job project (eg on dashboard)
config.Data.outputDatasetTag = partTwoUniqueNameMGP # name for storage directory beneath /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/<dataset>/

config.General.workArea = 'crab_projects' # name of directory where crab project info is stored
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'processMc02_genSimToAOD_step1of2_cfg.py'

config.Data.inputDataset = '/nmssmCascade/taylor-processMc01_nmssmSignalV01_TEST9_mH110p0_mSusy1800p0_ratio0p95_splitting1p0_25000events-fc4d6f674e85c70df8b656268b1e35aa/USER' # will need to make more general
config.Data.unitsPerJob = 1 # to correspond to a single unit output in partOne of the process
config.Data.totalUnits = -1 # '-1' does them all
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.publication = True

config.Site.storageSite = 'T2_UK_SGrid_Bristol'


### ## # EXTRA INFO # ## ###
#
# this crab job turns CMSSW_GENSIM files into INTERMEDIATE2AOD ROOT files
# part 2/3 of making CMSSW_RECO data from MG output
#
# the public location of the ouput will go here:
# /<outputPrimaryDatasetFROMPART1>/taylor-<partTwoUniqueNameMGP>-<randomHash>/USER
#
# output ROOT files will go here:
# /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/CRAB_PrivateMC/<partTwoUniqueNameMGP>/<dateStamp>/0000/
#
# the names of the output ROOT files will be:
# nmssmSignal_AODSIMstep1of2_<jobNum>.root
# as stated by the config file processMc01_mgLheToGenSim_cfg.py
############################