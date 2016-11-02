import os

# run with
# $ python $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/python/mcCampaign.py > $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/python/mcCampaignCommands.txt
# copy and paste the output into madgrpah command line (can do it in one big go)

### ## # USER INPUTS # ## ###
### where you wish to store your projects (too big for /users)
projectPath = "/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01" 
# projectPath = "/users/jt15104/MG5_aMC_v2_4_2/nmssmCascadeAnalysis_v01" # FILES TOO BIG FOR THIS
# projectPath = "/afs/cern.ch/work/t/taylor/public/madGraphProjects/nmssmCascadeAnalysis_v01"

### where madgraph is installed
madGraphPath = "/users/jt15104/MG5_aMC_v2_4_2/"
# madGraphPath = "/afs/cern.ch/user/t/taylor/MG5_aMC_v2_4_3/"

createNewProject = True # True if want to generate the new physics
# createNewProject = False # if using project that already exists
# higgsMassScan = [70.0, 80.0, 90.0, 100.0, 120.0, 130.0]
higgsMassScan = [85.0]
# susyMassScan = [800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0]
susyMassScan = [1000.0]
massRatio = 0.95 # massHiggs / massNLSP
massSplitting = 1.0 # massNLSP - massHiggs - massLSP
numberEvents = 1
#############################

if createNewProject:
	print("import model nmssmCascadeAnalysis;")
	print("define squark = ul ur dl dr cl cr sl sr ul~ ur~ dl~ dr~ cl~ cr~ sl~ sr~;")
	print("generate p p > squark squark;")
	print("add process p p > squark go;")
	print("add process p p > go go;")
	print("output %s;" % projectPath)

for higgsMass in higgsMassScan:
	nlspMass = higgsMass / massRatio
	lspMass = nlspMass - higgsMass - massSplitting

	for susyMass in susyMassScan:
		gluinoMass = 1.01 * susyMass

		outputName = "mH" + str(higgsMass) + "_mSusy" + str(susyMass) + "_ratio" + str(massRatio) + "_splitting" + str(massSplitting) + "_" + str(numberEvents) + "events"
		outputName = outputName.replace(".", "p")

		print("launch %s -n %s;" % (projectPath, outputName))
		print("pythia=ON")
		print("delphes=ON")
		print("done")
		print("%s/Cards/delphes_card_CMS.dat" % projectPath)
		print("%smodels/nmssmCascadeAnalysis/P1.param_card" % madGraphPath)
		print("set mass mneu1 %f" % lspMass)
		print("set mass mneu2 %f" % nlspMass)
		print("set mass mh01 %f" % higgsMass)
		print("set mass msd1 %f" % susyMass)
		print("set mass msd2 %f" % susyMass)
		print("set mass msd4 %f" % susyMass)
		print("set mass msd5 %f" % susyMass)
		print("set mass msu1 %f" % susyMass)
		print("set mass msu2 %f" % susyMass)
		print("set mass msu4 %f" % susyMass)
		print("set mass msu5 %f" % susyMass)
		print("set mass mgo %f" % gluinoMass)
		print("set nevents %d" % numberEvents)
		print("done")