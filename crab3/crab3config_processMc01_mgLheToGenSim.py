
from CRABClient.UserUtilities import config
config = config()


uniqueName = 'processMc01_nmssmSignalV01_TEST' # name for this production process (include commit hash so it is unique and can be traced) 
config.General.requestName = uniqueName # name of the crab job project (eg on dashboard)
config.Data.outputDatasetTag = uniqueName # name for storage directory beneath /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/<dataset>/

config.General.workArea = 'crab_projects' # name of directory where crab project info is stored
config.General.transferOutputs = True
config.General.transferLogs    = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'processMc01_mgLheToGenSim_cfg.py'
config.JobType.inputFiles = '/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/paramCard_type02/mH110p0_mSusy1800p0_ratio0p95_splitting1p0_25000events/Events/run_01/events.lhe'


# config.Data.inputDataset         = '/NNPDF30_13TeV_Pseudoscalar_100_1-LHE/sbreeze-ICL-POWHEG_DMS_NNPDF30_13TeV_Pseudoscalar_100_1-LHE-d1e6a7a52c9e635a587ef1faf77a1500/USER' # insert dataset
config.Data.inputDBS             = 'global' # global, phys03
config.Data.publication = False

config.Data.splitting            = 'EventBased' 
config.Data.unitsPerJob          = 10 # number of events per job
config.Data.totalUnits           = 40 # -1 for all


config.Site.storageSite = 'T2_UK_SGrid_Bristol'
