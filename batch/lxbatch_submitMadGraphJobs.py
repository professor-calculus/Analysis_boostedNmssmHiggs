import os

#############################
### run with
### $ python $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/batch/lxbatch_submitMadGraphJobs.py
### (best to run from within a jobs log dir)
### makes .txt file madgraph instructions and then uses these to submit jobs to lxbatch

#############################
### ## # USER INPUTS # ## ###
### where is the mad-graph project stored
projectPath = "/afs/cern.ch/work/t/taylor/public/madGraphProjects/nmssmCascadeAnalysis_v01"

### where the parameter cards are stored (must run madGraphParamCardGenerator.py first)
paramCardDir = "/afs/cern.ch/user/t/taylor/MG5_aMC_v2_3_3/nmssmCascadeParamCards" 

### parameter options
higgsMassScan = [70.0, 80.0]
susyMassScan = [800.0, 1000.0]
massRatio = 0.95 # massHiggs / massNLSP
massSplitting = 1.0 # massNLSP - massHiggs - massLSP
numberEvents = 1000
#############################

os.system("mkdir %s/productionInstructions/" % projectPath)

for higgsMass in higgsMassScan:
	for susyMass in susyMassScan:

		temp = "mH" + str(higgsMass) + "_mSusy" + str(susyMass) + "_ratio" + str(massRatio) + "_splitting" + str(massSplitting)
		temp = temp.replace(".", "p")
		paramCard = paramCardDir + "/" + temp + ".param_card"
		temp = temp + "_" + str(numberEvents) + "events"
		instructionText = projectPath + "/productionInstructions/" + temp + ".txt"
		
		f = open(instructionText, 'w')
		f.write("launch %s -n %s;\n" % (projectPath,temp))
		f.write("pythia=ON\n")
		f.write("delphes=ON\n")
		f.write("done\n")
		f.write("%s/Cards/delphes_card_CMS.dat\n" % projectPath)
		f.write("%s\n" % paramCard)
		f.write("set nevents %d\n" % numberEvents)
		f.write("done\n")
		f.close()

		os.system('bsub -q 8nh "sh $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/batch/lxbatch_madGraphJob.sh %s"' % instructionText)
