from CRABClient.UserUtilities import config
import os
config = config()

# NOTE THAT THIS NOW EXISTS PRIMARILY AS A RECORD
# SHOULD BE USING mcProduction_Controller.py FOR THIS WORKFLOW

# submit, from crab directory, with $ crab submit -c crab3config_processMc02_genSimToAOD_step1of2.py
# but first...

# 1. EDIT THE FOLLOWING #
madGraphProjectName = 'mH110p0_mSusy1800p0_ratio0p95_splitting1p0_25000events'
partThreeUniqueName = 'processMc03_nmssmSignalV01_F5'

partThreeUniqueNameMGP = partThreeUniqueName + "_" + madGraphProjectName # combine info for a truly unique record

config.General.requestName = partThreeUniqueNameMGP # name of the crab job project (eg on dashboard)
config.Data.outputDatasetTag = partThreeUniqueNameMGP # name for storage directory beneath /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/<dataset>/

config.General.workArea = 'crab_projects' # name of directory where crab project info is stored
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'processMc03_genSimToAOD_step2of2_cfg.py'

config.Data.inputDataset = '/nmssmCascade/?/USER' # will need to make more general
config.Data.unitsPerJob = 1 # to correspond to a single unit output in partOne of the process
config.Data.totalUnits = -1 # '-1' does them all
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.publication = True

config.Site.storageSite = 'T2_UK_SGrid_Bristol'


### ## # EXTRA INFO # ## ###
#
# this crab job turns CMSSW_GENSIM files into INTERMEDIATE2AOD ROOT files
# part 3/3 of making CMSSW_RECO data from MG output
#
# the public location of the ouput will go here:
# /<outputPrimaryDatasetFROMPART1>/taylor-<partThreeUniqueNameMGP>-<randomHash>/USER
#
# output ROOT files will go here:
# /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/CRAB_PrivateMC/<partThreeUniqueNameMGP>/<dateStamp>/0000/
#
# the names of the output ROOT files will be:
# nmssmSignal_AODSIMstep2of2_<jobNum>.root
# as stated by the config file processMc01_mgLheToGenSim_cfg.py
############################