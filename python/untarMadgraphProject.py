import os
import sys
import re

#############################
### run with
### $ python $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/python/untarMadgraphProject.py
### 
### untars the MG projects
### MUST SET IT OFF FROM THE DIRECTORY WHERE THE FILES ARE!!!
#############################

#############################
#############################
### ## # USER INPUTS # ## ###

motherDir = "/storage/jt15104/madGraphProjects/nmssmCascadeAnalysis_v01/paramCard_type02/"

#############################
#############################
#############################

madGraphProjects = os.popen("ls %s" % motherDir, "r").readlines()

for madGraphProject in madGraphProjects:

	if madGraphProject.rstrip()[-7:] == ".tar.gz":
		projectToUnzip = os.path.join(motherDir,madGraphProject.rstrip())
		# print projectToUnzip
		os.system("tar -xvf %s" % projectToUnzip)