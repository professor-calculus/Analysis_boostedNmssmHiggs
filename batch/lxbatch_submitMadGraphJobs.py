import os

#############################
### run with
### $ python $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/batch/lxbatch_submitMadGraphJobs.py
### (best to run from within a jobs log dir)
### copies the base madgraph project (have a baseMGR dir to being with)
### creates the .txt file instructions for madGraph
### submit jobs to lxbatch


#############################
#############################
### ## # USER INPUTS # ## ###
#
### where is the mad-graph project stored
projectLocation = "/afs/cern.ch/work/t/taylor/public/madGraphProjects/nmssmCascadeAnalysis_v05/paramCard_type03"
#
### where the parameter cards are stored (must run madGraphParamCardGenerator.py first to create these!)
paramCardDir = "/afs/cern.ch/user/t/taylor/MG5_aMC_v2_3_3/nmssmCascadeParamCards/type_03" 
#
### parameter options
higgsMassScan = [70.0, 90.0, 110.0, 125.0]
# higgsMassScan = [70.0]
susyMassScan = [1000.0, 1400.0, 1800.0, 2200.0]
# susyMassScan = [1000.0]
massRatio = 0.99 # massHiggs / massNLSP
massSplitting = 0.5 # massNLSP - massHiggs - massLSP
numberEvents = 600000

# usePythiaAndDelplhes = True
usePythiaAndDelplhes = False

# Only need True if you are Matching in MadGraphs built in pythia
# willWeNeedMatching = True
willWeNeedMatching = False
#############################
#############################
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
				f.write("launch %s;\n" % MGProject)
				if usePythiaAndDelplhes == True:
					f.write("pythia=ON\n")
					f.write("delphes=ON\n")
				f.write("done\n")
				if usePythiaAndDelplhes == True:
					f.write("%s/Cards/delphes_card_CMS.dat\n" % MGProject)
				f.write("%s\n" % paramCard)
				f.write("set nevents %d\n" % numberEvents)
				if willWeNeedMatching == True:
					f.write("set ickkw 1\n")
				f.write("done\n")
				f.close()

				os.system('bsub -q 8nh "sh $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/batch/lxbatch_madGraphJob.sh %s %s %s"' % (instructionText, projectLocation, temp))
				# print('bsub -q 8nh "sh $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/batch/lxbatch_madGraphJob.sh %s %s %s"' % (instructionText, projectLocation, temp))
				print("job %s submitted to lxbatch" % temp)

			else:
				print("job %s not submitted to lxbatch, the directory already exists:" % temp)

else:
	print("base project %s/baseMGR/ does not exist: Exiting..." % projectLocation)
