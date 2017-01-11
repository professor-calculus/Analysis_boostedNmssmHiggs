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
# executable = "$CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/McSignalStudiesCMSSW/McSignalStudiesCMSSW"
executable = "$CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies/DoubleBTaggerEfficiencyStudies"

# input files to run over (the equiv of what you would get back form 'ls', thus can use wildcards (*,?) and multiple arguments)
filesToUse = "/hdfs/user/jt15104/Analysis_boostedNmssmHiggs/patTuples/CMSSW_8_0_20/signalSamples/nmssmSignalCascadeV05_13TeV_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/nmssmSignalCascadeV05_13TeV_patTupleAddBTag_ed12_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/bTagPatTuple_1*.root /hdfs/user/jt15104/Analysis_boostedNmssmHiggs/patTuples/CMSSW_8_0_20/signalSamples/nmssmSignalCascadeV05_13TeV_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/nmssmSignalCascadeV05_13TeV_patTupleAddBTag_ed12_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/bTagPatTuple_2*.root /hdfs/user/jt15104/Analysis_boostedNmssmHiggs/patTuples/CMSSW_8_0_20/signalSamples/nmssmSignalCascadeV05_13TeV_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/nmssmSignalCascadeV05_13TeV_patTupleAddBTag_ed12_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/bTagPatTuple_3*.root"

# the NEW output DIR & file in which you wish to create to hold the output
# outputFile = "/users/jt15104/CMSSW_8_0_20/src/Analysis/Analysis_boostedNmssmHiggs/testingABC/test.root"
outputFile = "/users/jt15104/local_Analysis_boostedNmssmHiggs/output_DoubleBTaggerEfficiencyStudies/CMSSW_8_0_20/nmssmSignalCascadeV05_13TeV_mH70p0_mSusy1000p0_ratio0p99_splitting0p5/histos.root"

#############################
##########################################################
##########################################################

commandToRun = executable + " outputfile=" + outputFile + " inputfiles="

listOfRootFiles = os.popen("ls %s" % filesToUse, "r").readlines()
for rootFile in listOfRootFiles:
	rootFile = rootFile.rstrip()
	commandToRun = commandToRun + rootFile + ","

# print commandToRun
# print len(commandToRun)
os.system("%s" % commandToRun)