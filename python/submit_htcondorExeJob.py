import os
import sys
from htcondorExeJob import *


# set the user variables and then run with
# $ python $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/python/submit_htcondorExeJob.py
# FROM THE DIRECTORY WHERE YOU WISH TO SAVE THE ROOT FILES!


############################
############################
### ## # USER INPUT # ## ###
############################
############################
executable = "tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies/DoubleBTaggerEfficiencyStudies" # wrt CMSSW location
code = "src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies.cpp"
inputMotherDir = "/hdfs/user/jt15104/Analysis_boostedNmssmHiggs/patTuples/CMSSW_8_0_20/boostedHiggsSamples/GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph/boostedHiggsMC_25e777/" # should be on hdfs
logLocation = "/storage/jt15104/jobLog_Analysis_boostedNmssmHiggs/testing/"

# Note that you can put multiple files in a single element of vecInputFiles by just comma seperating them with no spaces (within the same set of "")
vecInputFiles = []
vecInputFiles.append(os.path.join(inputMotherDir,"bTagPatTuple_1.root"))
vecInputFiles.append(os.path.join(inputMotherDir,"bTagPatTuple_2.root"))
vecInputFiles.append(os.path.join(inputMotherDir,"bTagPatTuple_3.root"))
vecInputFiles.append(os.path.join(inputMotherDir,"bTagPatTuple_5.root"))
vecInputFiles.append(os.path.join(inputMotherDir,"bTagPatTuple_6.root"))

vecOutputFiles = []
vecOutputFiles.append("histos_1.root")
vecOutputFiles.append("histos_2.root")
vecOutputFiles.append("histos_3.root")
vecOutputFiles.append("histos_5.root")
vecOutputFiles.append("histos_6.root")

############################
############################
############################
############################
############################


if len(vecInputFiles) != len(vecOutputFiles):
	print "Length of input vector does not equal length of output vector"
	print "Exiting..."
	sys.exit()

os.system("mkdir -p %s" % logLocation)
executable = os.path.join(os.popen("echo $CMSSW_BASE").read().rstrip(),executable)
code = os.path.join(os.popen("echo $CMSSW_BASE").read().rstrip(),code)
os.system("cp %s ." % code)
os.system("cp $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/python/submit_htcondorExeJob.py .")

file = open("log.txt", "w")
for i in range (0,len(vecInputFiles)):

	job = htcondorExeJob()
	job.add_Executable(executable)
	job.add_Arguments("inputfiles=%s outputfile=%s runOnDice" % (vecInputFiles[i],vecOutputFiles[i])) # note the runOnDice argument required at the end
	job.add_Output("%s/htcondor.job.$(cluster).$(process).out" % logLocation)
	job.add_Log("%s/htcondor.job.$(cluster).$(process).log" % logLocation)
	job.add_Error("%s/htcondor.job.$(cluster).$(process).err" % logLocation)
	job.add_RequestMem("1000")
	job.add_Input("") # leave empty if you do not need them
	job.add_OtherInputFiles("") # leave empty if you do not need them
	# job.add_Input(ntuple.root)
	# job.add_OtherInputFiles("header1.h, header2.h")
	tempStr = job.submitJob()
	file.write(tempStr + "\n")
	file.write("input: " + vecInputFiles[i] + "\n")
	file.write("output: " + vecOutputFiles[i] + "\n\n")

file.close()
