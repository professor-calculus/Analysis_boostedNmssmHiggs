Information on how to run the tools in this repository
This repo is for nmssm cascade analysis.


Creating Data:

Monte Carlo.
Submit the script 'crab3/crab3config_MC.py' to the grid
$ crab submit -c crab3config_MC.py (need a valid grid certificate and to have loaded crab)
This script calls on 'python/patTuple_addBTagging_cfg.py'

'crab3/crab3config_MC.py'
uniqueName: to describe the job, good to include repo commit hash
DATASETS: config.Data.inputDataset and config.Data.unitsPerJob
(keep info about past datasets stored for reference)

'python/patTuple_addBTagging_cfg.py'
can control the detailed specifics of the crab jobs, eg what information to retain.




bin:
compile everything with $scram b -j8
this then give us executables in CMSSW_8_0_20/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/...
run with
$ executable options
(it runs with respect to where you are when you execute)

All the code utilises the FWLite command parser.
Mandatory Options:
inputfiles=file.root
inputfiles=file1.root,file2.root,file3.root (needs to be comma seperated without a space)

outputfile=outputfile.root
this should be the full path and include a new directory (which the code creates) where the file will go.
The new directory should include useful info about what you were doing.

Other options (don't need to explicitly write):
maxevents=-1 can set the max events (-1 is all of them)
outputevery=1000 how often we print to screen

orderedsecondaryfiles=0 (default)
orderedsecondaryfiles=1 this does the plotting part only
ie you have the histograms in the root file already and we run on them to remake the plots, for example changing the aesthetics
it runs on the .root file stated by outputfile=outputfile.root (so in this case this is the input) 

eg.1 STANDARD
$ ~/CMSSW_8_0_20/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies/DoubleBTaggerEfficiencyStudies inputfiles=XYZ_1.root,XYZ_2.root outputfile=ABC/out.root

eg.2 REMAKE THE PLOTS
$ ~/CMSSW_8_0_20/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies/DoubleBTaggerEfficiencyStudies outputfile=ABC/out.root



'bin/DoubleBTaggerEfficiencyStudies.cpp'
To run on boosted higgs mc samples. It looks at the effiency and properties of the boosted h->bb tagger.
it uses the class 'interface/PlottingDoubleBTaggerEfficiencyStudies.h' to create the plots

