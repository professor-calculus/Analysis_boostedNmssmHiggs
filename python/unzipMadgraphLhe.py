import os
import sys
import re

#############################
### run with
### $ python $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/python/unzipMadgraphLhe.py
### 
### unzips MG .lhe files (so that we can then put them through CMSSW)
#############################

#############################
#############################
### ## # USER INPUTS # ## ###

motherDir = "/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/paramCard_type01/"
pathWithinMadgraphProject = "Events/run_01/events.lhe.gz"

#############################
#############################
#############################

madGraphProjects = os.popen("ls %s" % motherDir, "r").readlines()

for madGraphProject in madGraphProjects:
	fileToUnzip = os.path.join(motherDir,madGraphProject.rstrip(),pathWithinMadgraphProject)
	# print fileToUnzip
	os.system("gunzip %s" % fileToUnzip)