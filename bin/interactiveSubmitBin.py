import os

#############################
### run with
### $ python $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/bin/interactiveSubmitBin.py

# remember to $ scram b first!!!
# script exists to allow easy running over multiple files
# does the full script with plotting (if you wish to just to the plotting again just run from the command line, don't use this)
##########################################################
##########################################################
### ## # USER INPUTS # ## ###

# the executable you wish to use
executable = "$CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/McSignalStudiesCMSSW/McSignalStudiesCMSSW"

# input files to run over, can use wildcards (*,?)
filesToUse = "/hdfs/user/jt15104/Analysis_boostedNmssmHiggs/patTuples/CMSSW_8_0_20/signalSamples/nmssmSignalCascadeV05_13TeV_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/nmssmSignalCascadeV05_13TeV_patTupleAddBTag_ed12_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/bTagPatTuple_88?.root"

# the output dir you wish to create to hold the output
outputDir = "/users/jt15104/CMSSW_8_0_20/src/Analysis/Analysis_boostedNmssmHiggs/testingABC/"

#############################
##########################################################
##########################################################

commandToRun = executable + " outputfile=" + outputDir + " inputfiles="

listOfRootFiles = os.popen("ls %s" % filesToUse, "r").readlines()
for rootFile in listOfRootFiles:
	rootFile = rootFile.rstrip()
	commandToRun = commandToRun + rootFile + ","

# print commandToRun
os.system("%s" % commandToRun)