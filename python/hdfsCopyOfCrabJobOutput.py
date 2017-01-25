import os
import sys
import re

# makes for easier copying of root files from crab3 output on /hdfs to your local area
#############################
### run with
### $ python $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/python/hdfsCopyOfCrabJobOutput.py

# currently assumes < 1000 root files so that the 0000 dir is always selected

####################
####################
### USER OPTIONS ###

saveBaseDir = "/hdfs/user/jt15104/Analysis_boostedNmssmHiggs/patTuples/CMSSW_8_0_21/signalSamples/"

crabOutputBaseDir = "/hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/"

crabOutputProjectList = [
                        "nmssmSignalCascadeV05_13TeV_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/nmssmSignalCascadeV05_13TeV_patTupleAddBTag_ed8021v1_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/",
                        ]

####################
####################
####################

# 1. check the saveBaseDir exists
if not os.path.isdir("%s" % (saveBaseDir)):
	print "saveBaseDir does not exist...Exiting..."
	sys.exit()

# 2. check the crabOutputBaseDir exists
if not os.path.isdir("%s" % (crabOutputBaseDir)):
	print "crabOutputBaseDir does not exist...Exiting..."
	sys.exit()

# 3. loop through the crabOutputProjectList
for crabOutputProject in crabOutputProjectList:

	# # 3a. check the directory we will make does not exist
	newDirectoryForRootFiles = os.path.join(saveBaseDir,crabOutputProject)
	if os.path.isdir("%s" % newDirectoryForRootFiles):
		print "directory already exists: " + newDirectoryForRootFiles
		print "NOT doing any copying for this project"
		print ""
		continue

	# 3b. check the directory we want to copy from exists
	dirToCopyFrom = os.path.join(crabOutputBaseDir,crabOutputProject)
	if not os.path.isdir("%s" % dirToCopyFrom):
		print "directory does not exist: " + newDirectoryForRootFiles
		print "NOT doing any copying for this project"
		print ""
		continue

	# 3c. create the 'newDirectoryForRootFiles'
	newDirectoryForRootFiles = newDirectoryForRootFiles[5:]
	# print newDirectoryForRootFiles
	os.system("hadoop fs -mkdir -p %s" % newDirectoryForRootFiles)


	# 3d. get the full path to the ROOT files
	nextPartInPath = os.popen("ls %s | head -1" % dirToCopyFrom, "r").readline() # most recent date stamp should be a larger number on crab, hence use of head
	nextPartInPath = nextPartInPath.rstrip()
	dirToCopyFrom = os.path.join(dirToCopyFrom,nextPartInPath,"0000")
	# print dirToCopyFrom

	# 3e. copy the root files from crab location to newly created hdfs location
	print "*************"
	print "*************"
	print "Copying Over:"
	print dirToCopyFrom
	print "To:"
	print newDirectoryForRootFiles
	print "*************"
	print "*************"
	crabRootFiles = dirToCopyFrom[5:] + "/*.root"
	# print crabRootFiles
	os.system("hadoop fs -cp %s %s" % (crabRootFiles,newDirectoryForRootFiles))
	print ""