import os

#############################
### run with
### $ python $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/bin/interactiveSubmitBin.py

# remember to $ scram b first!!!
# script exists to allow easy running over multiple files
# does the full script with plotting (if you wish to just to the plotting again just run from the command line, don't use this)
# NOTE THAT YOU SHOULD BE USING batch/submit_htcondorExeJob.py FOR FULL SCALE JOBS
##########################################################
##########################################################
### ## # USER INPUTS # ## ###

# the executable you wish to use
# executable = "$CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/McSignalStudiesCMSSW/McSignalStudiesCMSSW"
# executable = "$CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies/DoubleBTaggerEfficiencyStudies"
executable = "$CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/InvestigateEventSelection/InvestigateEventSelection"

# input files to run over (the equiv of what you would get back form 'ls', thus can use wildcards (*,?) and multiple arguments)
filesToUse = "/hdfs/user/jt15104/Analysis_boostedNmssmHiggs/patTuples/CMSSW_8_0_21/dummyBackground/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/qcdSpring16_CMSSW8021_90f89e0426/bTagPatTuple_*.root"

# the NEW output DIR & file in which you wish to create to hold the output
# outputFile = "$CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/testingScript/histos.root"
outputFile = "/users/jt15104/local_Analysis_boostedNmssmHiggs/output_InvestigateEventSelection/100K_events/dummyBackground/CMSSW_8_0_21/qcd/histos.root"

# max number of events to run over (-1 does them all)
maxevents = 100000

#############################
##########################################################
##########################################################

commandToRun = executable + " outputfile=" + outputFile + " maxevents=" + str(maxevents) + " inputfiles=" 

listOfRootFiles = os.popen("ls %s" % filesToUse, "r").readlines()
# listOfRootFiles = os.popen("ls %s | head -485" % filesToUse, "r").readlines() # this is for cases where otherwise the command is too long for unix shell
for rootFile in listOfRootFiles:
	rootFile = rootFile.rstrip()
	commandToRun = commandToRun + rootFile + ","

# print commandToRun
# print len(commandToRun)
os.system("%s" % commandToRun)