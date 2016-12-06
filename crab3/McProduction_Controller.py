import os
import sys

# select the correct options!!!
# make sure you are in the Analysis_boostedNmssmHiggs/crab3 repo
# run with $ python McProduction_Controller.py

#################################################################
#################################################################
#################################################################
#################################################################
###### @ U S E R @ O P T I O N S @ ##############################


mode = 'submit'
# mode = 'resubmit'
# mode = 'checkStatus'

whichPartOfProcess = 'processMc01' # turns madgraph LHE into cmssw GENSIM
# whichPartOfProcess = 'processMc02' # step one of GENSIM into AOD
# whichPartOfProcess = 'processMc03' # step two of GENSIM into AOD
# whichPartOfProcess = 'patTupleAddBTag' # formats the AOD into patTuple form


#-----------------------------------------------
##### INFO constant workflow INFO ##############
madGraphProjects = [
                    'mH70p0_mSusy1000p0_ratio0p95_splitting1p0_25000events',
                    'mH70p0_mSusy1400p0_ratio0p95_splitting1p0_25000events',
                    'mH70p0_mSusy1800p0_ratio0p95_splitting1p0_25000events',
                    'mH70p0_mSusy2200p0_ratio0p95_splitting1p0_25000events',
                   ]
outputPrimaryDatasetIntro = 'nmssmSignalCascadeV01_13TeV'

storageSite = 'T2_UK_SGrid_Bristol'
#-----------------------------------------------

#-------------------------------------------
##### INFO 'processMc01' INFO ##############
editionNamePro01 = "ed01"

eventsPerJob = 20
totalNumberOfEvents = 100 # -1 to select them all
localMadGraphProjectStore = '/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/paramCard_type02/'
pathWithinMadgraphProject = 'Events/run_01/events.lhe' # ensure that you have unzipped these files using python/unzipMadgraphLhe.py
#-------------------------------------------

#------------------------------------------- Note that running submission of processMc02
##### INFO 'processMc02' INFO ############## requires valid editionNamePro01
editionNamePro02 = "ed01"

filesPerJobPro02 = 1
totalNumberOfFilesPro02 = -1 # -1 to select them all
#-------------------------------------------

#------------------------------------------- Note that running submission of processMc03
##### INFO 'processMc03' INFO ############## requires valid editionNamePro02
editionNamePro03 = "ed01"

filesPerJobPro03 = 1
totalNumberOfFilesPro03 = -1 # -1 to select them all
#-------------------------------------------

#------------------------------------------- Note that running submission of patTupleAddBTag
##### INFO 'patTupleAddBTag' INFO ########## requires valid editionNamePro03
editionNamePAT = "ed01"

filesPerJobPAT = 1
totalNumberOfFilesPAT = -1 # -1 to select them all
#-------------------------------------------


#################################################################
#################################################################
### ## # HELPFUL INFO # ## ###
#
# outputPrimaryDataset = outputPrimaryDatasetIntro + '_' + madGraphProjectsStripOffEvents[i]
# partOneUniqueName = outputPrimaryDatasetIntro + '_processMc01_' + editionNamePro01 + '_' + madGraphProjectsStripOffEvents[i]
# partTwoUniqueName = outputPrimaryDatasetIntro + '_processMc02_' + editionNamePro02 + '_' + madGraphProjectsStripOffEvents[i]
# partThreeUniqueName = outputPrimaryDatasetIntro + '_processMc03_' + editionNamePro03 + '_' + madGraphProjectsStripOffEvents[i]
# patTupleUniqueName = outputPrimaryDatasetIntro + '_patTupleAddBTag_' + editionNamePAT + '_' + madGraphProjectsStripOffEvents[i]
# (the unique names must be less than 100 characters)
#
# the public location of the ouput will go here:
# 1. /<outputPrimaryDataset>/taylor-<partOneUniqueName>-<randomHash>/USER
# 2. /<outputPrimaryDataset>/taylor-<partTwoUniqueName>-<randomHash>/USER
# 3. /<outputPrimaryDataset>/taylor-<partThreeUniqueName>-<randomHash>/USER
# 4. /<outputPrimaryDataset>/taylor-<patTupleUniqueName>-<randomHash>/USER
#
# crab projects will go here:
# 1. ./crab_projects/crab_<partOneUniqueName>/
# 2. ./crab_projects/crab_<partTwoUniqueName>/
# 3. ./crab_projects/crab_<partThreeUniqueName>/
# 4. ./crab_projects/crab_<patTupleUniqueName>/
#
# output ROOT files will go here:
# 1. /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/<outputPrimaryDataset>/<partOneUniqueName>/<dateStamp>/0000/
# 2. /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/<outputPrimaryDataset>/<partTwoUniqueName>/<dateStamp>/0000/
# 3. /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/<outputPrimaryDataset>/<partThreeUniqueName>/<dateStamp>/0000/
# 4. /hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/<outputPrimaryDataset>/<patTupleUniqueName>/<dateStamp>/0000/
#
# the names of the output ROOT files will be:
# 1. nmssmSignal_GENSIM_<jobNum>.root (as stated by the config file processMc01_mgLheToGenSim_cfg.py)
# 2. nmssmSignal_AODSIMstep1of2_<jobNum>.root (as stated by the config file processMc02_genSimToAOD_step1of2_cfg.py)
# 3. nmssmSignal_AODSIMstep2of2_<jobNum>.root (as stated by the config file processMc03_genSimToAOD_step1of2_cfg.py)
# 4. bTagPatTuple_<jobNum>.root (as stated by the config file patTuple_addBTagging_cfg.py)
###############################

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################

# check that you submit from the correct directory
cmsswBase = os.popen("echo $CMSSW_BASE", "r").readline()
dirShouldBe = cmsswBase.rstrip() + "/src/Analysis/Analysis_boostedNmssmHiggs/crab3"
if os.getcwd() != dirShouldBe:
	print "We are not in the correct directory to run this script"
	print "Get yourself in $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/crab3"
	print "You Fool"
	sys.exit()

# load crab3
os.system("source /cvmfs/cms.cern.ch/crab3/crab_standalone.sh")

# load GRID certificate for a week
# os.system("voms-proxy-init -voms cms --valid 168:00")

# strip events section off madGraph info
madGraphProjectsStripOffEvents = []
for madGraphProject in madGraphProjects:
	for c in range(len(madGraphProject)-1, 0, -1):
		if madGraphProject[c] == '_':
			madGraphProjectsStripOffEvents.append(madGraphProject[:c])
			break




#-----------------------------------------------------------#
#----------------------processMc01--------------------------#
#-----------------------------------------------------------#

if mode == 'submit' and whichPartOfProcess == 'processMc01':

	for i in range(0,len(madGraphProjects)):

		outputPrimaryDataset = outputPrimaryDatasetIntro + '_' + madGraphProjectsStripOffEvents[i]
		partOneUniqueName = outputPrimaryDatasetIntro + '_processMc01_' + editionNamePro01 + '_' + madGraphProjectsStripOffEvents[i]
		inputFiles = os.path.join(localMadGraphProjectStore,madGraphProjects[i],pathWithinMadgraphProject)
	
		# create the tempory crab config file to submit
		f = open("temp_crab3config_processMc01.py", 'w')
		f.write("from CRABClient.UserUtilities import config\n")
		f.write("config = config()\n")
		f.write("config.Data.unitsPerJob = %d\n" % eventsPerJob)
		f.write("config.Data.totalUnits = %d\n" % totalNumberOfEvents)
		f.write("config.Data.splitting = 'EventBased'\n")
		f.write("config.Data.inputDBS = 'global'\n")
		f.write("config.Data.publication = True\n")
		f.write("config.Data.outputPrimaryDataset = '%s'\n" % outputPrimaryDataset)
		f.write("config.Data.outputDatasetTag = '%s'\n" % partOneUniqueName)
		f.write("config.General.requestName = '%s'\n" % partOneUniqueName)
		f.write("config.General.workArea = 'crab_projects'\n")
		f.write("config.General.transferOutputs = True\n")
		f.write("config.General.transferLogs    = True\n")
		f.write("config.JobType.pluginName = 'PrivateMC'\n")
		f.write("config.JobType.psetName = 'processMc01_mgLheToGenSim_cfg.py'\n")
		f.write("config.JobType.inputFiles = ['%s']\n" % inputFiles)
		f.write("config.Site.storageSite = '%s'\n" % storageSite)
		f.close()
		print ""
		# os.system("cat temp_crab3config_processMc01.py") # for testing
		os.system("crab submit -c temp_crab3config_processMc01.py") # for the real deal
		os.system("rm temp_crab3config_processMc01.py")
		print ""



if mode == 'checkStatus' and whichPartOfProcess == 'processMc01':
	print ""
	print "*** CHECKING STATUS FOR PROCESSMC01 ***"
	print "NB if no info printed the task is most likely not bootstrapped yet"
	for i in range(0,len(madGraphProjects)):
		print ""
		print ""
		partOneUniqueName = outputPrimaryDatasetIntro + '_processMc01_' + editionNamePro01 + '_' + madGraphProjectsStripOffEvents[i]
		statusLines = os.popen("crab status -d crab_projects/crab_%s" % partOneUniqueName, "r").readlines()
		print partOneUniqueName
		foundJobLine = False
		for line in statusLines:
			if foundJobLine == False and line[:12] == "Jobs status:":
				print line.rstrip()
				foundJobLine = True
				continue
			if foundJobLine == True:
				if line == "\n":
					break
				print line.rstrip()
		foundPubLine = False
		for line in statusLines:
			if foundPubLine == False and line[:19] == "Publication status:":
				print line.rstrip()
				foundPubLine = True
				continue
			if foundPubLine == True:
				if line == "\n":
					break
				print line.rstrip()
	print ""



if mode == 'resubmit' and whichPartOfProcess == 'processMc01':

	for i in range(0,len(madGraphProjects)):
		partOneUniqueName = outputPrimaryDatasetIntro + '_processMc01_' + editionNamePro01 + '_' + madGraphProjectsStripOffEvents[i]
		os.system("crab resubmit -d crab_projects/crab_%s" % partOneUniqueName)



		
#-----------------------------------------------------------#
#----------------------processMc02--------------------------#
#-----------------------------------------------------------#

if mode == 'submit' and whichPartOfProcess == 'processMc02':

	for i in range(0,len(madGraphProjects)):
		partOneUniqueName = outputPrimaryDatasetIntro + '_processMc01_' + editionNamePro01 + '_' + madGraphProjectsStripOffEvents[i]
		partTwoUniqueName = outputPrimaryDatasetIntro + '_processMc02_' + editionNamePro02 + '_' + madGraphProjectsStripOffEvents[i]

		# get the name of the DAS input dataset name
		inputDataset = []
		statusLines = os.popen("crab status -d crab_projects/crab_%s" % partOneUniqueName, "r").readlines()
		for line in statusLines:
			if line[:15] == "Output dataset:":
				inputDataset = line.rstrip()
				for c in range(15,len(inputDataset)):
					if inputDataset[c] == '/':
						inputDataset = inputDataset[c:]
						# print inputDataset # for debugging
						break
				break

		# create the tempory crab config file to submit
		f = open("temp_crab3config_processMc02.py", 'w')
		f.write("from CRABClient.UserUtilities import config\n")
		f.write("config = config()\n")
		f.write("config.Data.inputDataset = '%s'\n" % inputDataset)
		f.write("config.Data.unitsPerJob = %d\n" % filesPerJobPro02)
		f.write("config.Data.totalUnits = %d\n" % totalNumberOfFilesPro02)
		f.write("config.Data.inputDBS = 'phys03'\n")
		f.write("config.Data.splitting = 'FileBased'\n")
		f.write("config.Data.publication = True\n")
		f.write("config.Data.outputDatasetTag = '%s'\n" % partTwoUniqueName)
		f.write("config.General.requestName = '%s'\n" % partTwoUniqueName)
		f.write("config.General.workArea = 'crab_projects'\n")
		f.write("config.General.transferOutputs = True\n")
		f.write("config.General.transferLogs = True\n")
		f.write("config.JobType.pluginName = 'Analysis'\n")
		f.write("config.JobType.psetName = 'processMc02_genSimToAOD_step1of2_cfg.py'\n")
		f.write("config.Site.storageSite = '%s'\n" % storageSite)
		f.close()
		print ""
		# os.system("cat temp_crab3config_processMc02.py") # for testing
		os.system("crab submit -c temp_crab3config_processMc02.py") # for the real deal
		os.system("rm temp_crab3config_processMc02.py")
		print ""



if mode == 'checkStatus' and whichPartOfProcess == 'processMc02':
	print ""
	print "*** CHECKING STATUS FOR PROCESSMC02 ***"
	print "NB if no info printed the task is most likely not bootstrapped yet"
	for i in range(0,len(madGraphProjects)):
		print ""
		print ""
		partTwoUniqueName = outputPrimaryDatasetIntro + '_processMc02_' + editionNamePro02 + '_' + madGraphProjectsStripOffEvents[i]
		statusLines = os.popen("crab status -d crab_projects/crab_%s" % partTwoUniqueName, "r").readlines()
		print partTwoUniqueName
		foundJobLine = False
		for line in statusLines:
			if foundJobLine == False and line[:12] == "Jobs status:":
				print line.rstrip()
				foundJobLine = True
				continue
			if foundJobLine == True:
				if line == "\n":
					break
				print line.rstrip()
		foundPubLine = False
		for line in statusLines:
			if foundPubLine == False and line[:19] == "Publication status:":
				print line.rstrip()
				foundPubLine = True
				continue
			if foundPubLine == True:
				if line == "\n":
					break
				print line.rstrip()
	print ""



if mode == 'resubmit' and whichPartOfProcess == 'processMc02':

	for i in range(0,len(madGraphProjects)):
		partTwoUniqueName = outputPrimaryDatasetIntro + '_processMc02_' + editionNamePro02 + '_' + madGraphProjectsStripOffEvents[i]
		os.system("crab resubmit -d crab_projects/crab_%s" % partTwoUniqueName)



		
#-----------------------------------------------------------#
#----------------------processMc03--------------------------#
#-----------------------------------------------------------#

if mode == 'submit' and whichPartOfProcess == 'processMc03':

	for i in range(0,len(madGraphProjects)):
		partTwoUniqueName = outputPrimaryDatasetIntro + '_processMc02_' + editionNamePro02 + '_' + madGraphProjectsStripOffEvents[i]
		partThreeUniqueName = outputPrimaryDatasetIntro + '_processMc03_' + editionNamePro03 + '_' + madGraphProjectsStripOffEvents[i]

		# get the name of the DAS input dataset name
		inputDataset = []
		statusLines = os.popen("crab status -d crab_projects/crab_%s" % partTwoUniqueName, "r").readlines()
		for line in statusLines:
			if line[:15] == "Output dataset:":
				inputDataset = line.rstrip()
				for c in range(15,len(inputDataset)):
					if inputDataset[c] == '/':
						inputDataset = inputDataset[c:]
						# print inputDataset # for debugging
						break
				break

		# create the tempory crab config file to submit
		f = open("temp_crab3config_processMc03.py", 'w')
		f.write("from CRABClient.UserUtilities import config\n")
		f.write("config = config()\n")
		f.write("config.Data.inputDataset = '%s'\n" % inputDataset)
		f.write("config.Data.unitsPerJob = %d\n" % filesPerJobPro03)
		f.write("config.Data.totalUnits = %d\n" % totalNumberOfFilesPro03)
		f.write("config.Data.inputDBS = 'phys03'\n")
		f.write("config.Data.splitting = 'FileBased'\n")
		f.write("config.Data.publication = True\n")
		f.write("config.Data.outputDatasetTag = '%s'\n" % partThreeUniqueName)
		f.write("config.General.requestName = '%s'\n" % partThreeUniqueName)
		f.write("config.General.workArea = 'crab_projects'\n")
		f.write("config.General.transferOutputs = True\n")
		f.write("config.General.transferLogs = True\n")
		f.write("config.JobType.pluginName = 'Analysis'\n")
		f.write("config.JobType.psetName = 'processMc03_genSimToAOD_step2of2_cfg.py'\n")
		f.write("config.Site.storageSite = '%s'\n" % storageSite)
		f.close()
		print ""
		# os.system("cat temp_crab3config_processMc03.py") # for testing
		os.system("crab submit -c temp_crab3config_processMc03.py") # for the real deal
		os.system("rm temp_crab3config_processMc03.py")
		print ""



if mode == 'checkStatus' and whichPartOfProcess == 'processMc03':
	print ""
	print "*** CHECKING STATUS FOR PROCESSMC03 ***"
	print "NB if no info printed the task is most likely not bootstrapped yet"
	for i in range(0,len(madGraphProjects)):
		print ""
		print ""
		partThreeUniqueName = outputPrimaryDatasetIntro + '_processMc03_' + editionNamePro03 + '_' + madGraphProjectsStripOffEvents[i]
		statusLines = os.popen("crab status -d crab_projects/crab_%s" % partThreeUniqueName, "r").readlines()
		print partThreeUniqueName
		foundJobLine = False
		for line in statusLines:
			if foundJobLine == False and line[:12] == "Jobs status:":
				print line.rstrip()
				foundJobLine = True
				continue
			if foundJobLine == True:
				if line == "\n":
					break
				print line.rstrip()
		foundPubLine = False
		for line in statusLines:
			if foundPubLine == False and line[:19] == "Publication status:":
				print line.rstrip()
				foundPubLine = True
				continue
			if foundPubLine == True:
				if line == "\n":
					break
				print line.rstrip()
	print ""



if mode == 'resubmit' and whichPartOfProcess == 'processMc03':

	for i in range(0,len(madGraphProjects)):
		partThreeUniqueName = outputPrimaryDatasetIntro + '_processMc03_' + editionNamePro03 + '_' + madGraphProjectsStripOffEvents[i]
		os.system("crab resubmit -d crab_projects/crab_%s" % partThreeUniqueName)




#-----------------------------------------------------------#
#----------------------patTupleAddBTag----------------------#
#-----------------------------------------------------------#

if mode == 'submit' and whichPartOfProcess == 'patTupleAddBTag':

	for i in range(0,len(madGraphProjects)):
		partThreeUniqueName = outputPrimaryDatasetIntro + '_processMc03_' + editionNamePro03 + '_' + madGraphProjectsStripOffEvents[i]
		patTupleUniqueName = outputPrimaryDatasetIntro + '_patTupleAddBTag_' + editionNamePAT + '_' + madGraphProjectsStripOffEvents[i]

		# get the name of the DAS input dataset name
		inputDataset = []
		statusLines = os.popen("crab status -d crab_projects/crab_%s" % partThreeUniqueName, "r").readlines()
		for line in statusLines:
			if line[:15] == "Output dataset:":
				inputDataset = line.rstrip()
				for c in range(15,len(inputDataset)):
					if inputDataset[c] == '/':
						inputDataset = inputDataset[c:]
						# print inputDataset # for debugging
						break
				break

		# create the tempory crab config file to submit
		f = open("temp_crab3config_patTuple.py", 'w')
		f.write("from CRABClient.UserUtilities import config\n")
		f.write("config = config()\n")
		f.write("config.Data.inputDataset = '%s'\n" % inputDataset)
		f.write("config.Data.unitsPerJob = %d\n" % filesPerJobPAT)
		f.write("config.Data.totalUnits = %d\n" % totalNumberOfFilesPAT)
		f.write("config.Data.inputDBS = 'phys03'\n")
		f.write("config.Data.splitting = 'FileBased'\n")
		f.write("config.Data.publication = False\n")
		f.write("config.Data.outputDatasetTag = '%s'\n" % patTupleUniqueName)
		f.write("config.General.requestName = '%s'\n" % patTupleUniqueName)
		f.write("config.General.workArea = 'crab_projects'\n")
		f.write("config.General.transferOutputs = True\n")
		f.write("config.General.transferLogs = True\n")
		f.write("config.JobType.pluginName = 'Analysis'\n")
		f.write("config.JobType.psetName = 'patTuple_addBTagging_cfg.py'\n")
		f.write("config.Site.storageSite = '%s'\n" % storageSite)
		f.close()
		print ""
		os.system("cat temp_crab3config_patTuple.py") # for testing
		# os.system("crab submit -c temp_crab3config_patTuple.py") # for the real deal
		os.system("rm temp_crab3config_patTuple.py")
		print ""



if mode == 'checkStatus' and whichPartOfProcess == 'patTupleAddBTag':
	print ""
	print "*** CHECKING STATUS FOR PROCESSMC03 ***"
	print "NB if no info printed the task is most likely not bootstrapped yet"
	for i in range(0,len(madGraphProjects)):
		print ""
		print ""
		patTupleUniqueName = outputPrimaryDatasetIntro + '_patTupleAddBTag_' + editionNamePAT + '_' + madGraphProjectsStripOffEvents[i]
		statusLines = os.popen("crab status -d crab_projects/crab_%s" % patTupleUniqueName, "r").readlines()
		print partThreeUniqueName
		foundJobLine = False
		for line in statusLines:
			if foundJobLine == False and line[:12] == "Jobs status:":
				print line.rstrip()
				foundJobLine = True
				continue
			if foundJobLine == True:
				if line == "\n":
					break
				print line.rstrip()
		foundPubLine = False
		for line in statusLines:
			if foundPubLine == False and line[:19] == "Publication status:":
				print line.rstrip()
				foundPubLine = True
				continue
			if foundPubLine == True:
				if line == "\n":
					break
				print line.rstrip()
	print ""



if mode == 'resubmit' and whichPartOfProcess == 'patTupleAddBTag':

	for i in range(0,len(madGraphProjects)):
		patTupleUniqueName = outputPrimaryDatasetIntro + '_patTupleAddBTag_' + editionNamePAT + '_' + madGraphProjectsStripOffEvents[i]
		os.system("crab resubmit -d crab_projects/crab_%s" % patTupleUniqueName)


#################################################################
#################################################################
#################################################################
#################################################################
#################################################################