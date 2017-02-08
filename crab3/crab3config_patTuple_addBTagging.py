from CRABClient.UserUtilities import config
config = config()

# submit, from crab directory, with $ crab submit -c crab3config_patTuple_addBTagging.py

# EDIT the following
# uniqueName, config.Data.inputDataset, config.Data.unitsPerJob, (totalUnits)


uniqueName = 'ttbar_CMSSW8021_90f89e0426' # name for this ntuple production (include commit hash so it is unique and can be traced) 
config.General.requestName = uniqueName # name of the crab job project (eg on dashboard)
config.Data.outputDatasetTag = uniqueName # name for storage directory beneath /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/<dataset>/

config.General.workArea = 'crab_projects' # name of directory where crab project info is stored
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'patTuple_addBTagging_cfg.py'

#################################################################################################
# DATASETS --------------------------------------------------------------------------------------
# for testing the X->bb tagger on single boosted higgs (use 1 unitsPerJob, 1unit here ~=9k events)
# config.Data.inputDataset = '/GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph/RunIISpring16reHLT80-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v3/AODSIM'

# a dummy qcd background sample (use 1 unitsPerJob, 1unit here ~= 7k events) (there are ~10M events)
# config.Data.inputDataset = '/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISpring16DR80-FlatPU20to70HcalNZSRAW_withHLT_80X_mcRun2_asymptotic_v14-v1/AODSIM'

# a dummy ttbar background sample (use 1 unitsPerJob, 1unit here ~= 10.5k events) (there are ~10M events)
config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM'

#################################################################################################

config.Data.unitsPerJob = 1 # depends on the dataset ~1-10 (will list recommendation for datasets previously used)
config.Data.totalUnits = -1 # '-1' does them all
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.publication = False

config.Site.storageSite = 'T2_UK_SGrid_Bristol'