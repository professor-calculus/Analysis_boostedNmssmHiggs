import os

#############################
### run with
### $ python $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/batch/lxbatch_submitMadGraphJobs.py
### (best to run from within a jobs log dir)
### copies the base madgraph project (have a baseMGR and a testMGR project dir to being with)
### creates the .txt file instructions for madGraph
### submit jobs to lxbatch
#############################
### ## # USER INPUTS # ## ###
#
### where is the mad-graph project stored
projectLocation = "/afs/cern.ch/work/t/taylor/public/madGraphProjects/nmssmCascadeAnalysis_v01"
#
### where the parameter cards are stored (must run madGraphParamCardGenerator.py first to create these!)
paramCardDir = "/afs/cern.ch/user/t/taylor/MG5_aMC_v2_3_3/nmssmCascadeParamCards" 
#
### parameter options
# higgsMassScan = [70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0]
higgsMassScan = [70.0]
# susyMassScan = [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0]
susyMassScan = [1000.0]
massRatio = 0.95 # massHiggs / massNLSP
massSplitting = 1.0 # massNLSP - massHiggs - massLSP
numberEvents = 500
#############################

# check that the base directory does indeed exist
if os.path.isdir("%s/baseMGR/" % projectLocation):

	for higgsMass in higgsMassScan:
		for susyMass in susyMassScan:

			temp = "mH" + str(higgsMass) + "_mSusy" + str(susyMass) + "_ratio" + str(massRatio) + "_splitting" + str(massSplitting)
			temp = temp.replace(".", "p")
			paramCard = paramCardDir + "/" + temp + ".param_card"
			temp = temp + "_" + str(numberEvents) + "events"
			MGProject = projectLocation + "/" + temp 
			# check the dir does not exist already
			if not os.path.isdir("%s" % (MGProject)):

				# make the working directory by copying the base directory
				os.system("cp -r %s/baseMGR/ %s" % (projectLocation,MGProject))
				# create the instructions for the madgraph job (links to the corresponding .param_card)
				instructionText = MGProject + "/instructionsToRun.txt"		
				f = open(instructionText, 'w')
				f.write("launch %s;\n" % (MGProject))
				f.write("pythia=ON\n")
				f.write("delphes=ON\n")
				f.write("done\n")
				f.write("%s/Cards/delphes_card_CMS.dat\n" % MGProject)
				f.write("%s\n" % paramCard)
				f.write("set nevents %d\n" % numberEvents)
				f.write("done\n")
				f.close()

				os.system('bsub -q 8nh "sh $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/batch/lxbatch_madGraphJob.sh %s %s"' % (instructionText, MGProject))
				# print('bsub -q 8nh "sh $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/batch/lxbatch_madGraphJob.sh %s %s"' % (instructionText, MGProject))
				print("job %s submitted to lxbatch" % temp)

			else:
				print("job %s not submitted to lxbatch, the directory already exists:" % temp)

else:
	print("base project %s/baseMGR/ does not exist: Exiting..." % projectLocation)