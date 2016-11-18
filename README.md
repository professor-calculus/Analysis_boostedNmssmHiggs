Information on how to run the tools in this repository
This repo is for the nmssm cascade analysis.





Creating Data From the GRID:

Monte Carlo.
Submit the script 'crab3/crab3config_MC.py' to the grid. This is done like so:
$ crab submit -c crab3config_MC.py (need a valid grid certificate and to have loaded crab)
-grid certificate: $ voms-proxy-init -voms cms --valid 168:00
-crab: $ source /cvmfs/cms.cern.ch/crab3/crab.sh (shortcut $loadCrab)
The script 'crab3/crab3config_MC.py' calls on 'python/patTuple_addBTagging_cfg.py'

'crab3/crab3config_MC.py' info.
uniqueName: to describe the job, good to include repo commit hash so can trace back what you did.
datasets: config.Data.inputDataset and config.Data.unitsPerJob
(keep info about past datasets stored for reference)

'python/patTuple_addBTagging_cfg.py' info.
can control the detailed specifics of the crab jobs, e.g. what information to retain.
Currently returns a patTuple, maybe this is too larger type of file for the future...

For Bristol:
Crab jobs will go to 'dpm/phy.bris.ac.uk/home/cms/store/user/taylor/'
You should copy them to '/hdfs/user/jt15104/Analysis_boostedNmssmHiggs/CMSSW_8_X_YZ/sampleType/...' using haddop commands







bin:

Compile everything with $scram b -j8
This then gives us executables in CMSSW_8_0_20/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/...
run with
$ executable options
(it runs with respect to where you are when you execute)
(the executable is just the script name without the .cpp extension)
(for the options see the specific script)

All the code utilises the FWLite command parser.
***Mandatory Options***:

-Input of single files
inputfiles=file.root
-Input of multiple files
inputfiles=file1.root,file2.root,file3.root (needs to be comma seperated without a space)

outputfile=outputfile.root
this should be the full path and include a new directory (which the code creates) where the file will go.
The new directory should include useful info about what you were doing. The directory should also be unique
eg "output_DoubleBTaggerEfficiencyStudies_boostedHiggsToBB/histos.root"
Good idea to use the name of the dir where you got the ntuples.
Link clearly to what you used from /hdfs output_DoubleBTaggerEfficiencyStudies/GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph_boostedHiggsMC_25e777/
(typically use full paths to inputs and outputs)

***Other options (don't need to explicitly write)***:

maxevents=-1 can set the max events (-1 is all of them)

outputevery=1000 how often we print to screen

orderedsecondaryfiles=0 (default)
orderedsecondaryfiles=1 this does the plotting part only
ie you have the histograms in the root file already and we run on them to remake the plots, for example changing the aesthetics
it runs on the .root file stated by outputfile=outputfile.root (so in this case this is the input) 

eg.1 STANDARD
$ ~/CMSSW_8_0_20/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies/DoubleBTaggerEfficiencyStudies inputfiles=XYZ_1.root,XYZ_2.root outputfile=ABC/out.root

eg.2 REMAKE THE PLOTS
$ ~/CMSSW_8_0_20/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies/DoubleBTaggerEfficiencyStudies outputfile=ABC/histos.root orderedsecondaryfiles=1

nb: this code makes a copy of itself and puts it in the output directory
but if you have made changes after compilation this will be out of sync, watch out!
it also makes a copy of the parser values in a seperate text file (this will be in sync)
-together, along with non-lazy directory name, you should have decent records of what you have done:)

*** Explicit descriptions ***
'bin/DoubleBTaggerEfficiencyStudies.cpp'
To run on a boosted higgs mc samples. It looks at the effiency and properties of the boosted h->bb tagger.
it uses the class 'interface/PlottingDoubleBTaggerEfficiencyStudies.h' to create the plots

e.g. of running. (the bad practise here is that the patTuples have not been copied to '/hdfs/user/jt15104/')
$CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/DoubleBTaggerEfficiencyStudies/DoubleBTaggerEfficiencyStudies
inputFiles=/hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph/boostedHiggsMC_25e777/161011_112815/0000/bTagPatTuple_1.root,/hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph/boostedHiggsMC_25e777/161011_112815/0000/bTagPatTuple_2.root,/hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph/boostedHiggsMC_25e777/161011_112815/0000/bTagPatTuple_3.root,/hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph/boostedHiggsMC_25e777/161011_112815/0000/bTagPatTuple_5.root,/hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph/boostedHiggsMC_25e777/161011_112815/0000/bTagPatTuple_6.root
outputfile=~/local_Analysis_boostedNmssmHiggs/output_DoubleBTaggerEfficiencyStudies_boostedHiggsMC_25e777/histos.root
orderedsecondaryfiles=0

For now target all outputs to '~/local_Analysis_boostedNmssmHiggs/'. This may eventually take up to much space on /users







macros:
To be run interactively with ROOT.

massParamsToKinematics.cxx info.
Makes many plots to do with the kinematics of the decay chain by just doing calculations, not simulation.
Also produces .tex files which will automatically make pdfs with all the plots when run in latex locally.
There is a simple 'USER INPUT' section, and that is all you need to change.
$ root -l -b -q massParamsToKinematics.cxx


McSignalStudies.cxx info.
Analyses the output of the MadGraph/Pythia/Delphes production of the signal.
Has the option of just re-doing the plots and not having to run over everything again.
Takes the following arguements. The input files as a string vector and output files as a string.
(Link the output to where you got the ntuples on /storage)

-to submit on mass use submit_McSignalStudies.cxx
-to analysize different outputs together use combine_McSignalStudies.cxx







python:

madGraphParamCardGenerator.py info.
use to create param cards needed by my workflow to create MadGraph simulations of the data.

patTuple_addBTagging_cfg.py info.
Already described above.

submit_htcondorExeJob.py info.
WORK IN PROGRESS, NOT CURRENTLY FUNCTIONAL.
it will be used for submitting jobs to htcondor








batch:

lxbatch_submitMadGraphJobs.py info.
Use to submit Madgraph jobs onto lxplus.
Will scan over mass of higgs and mass susy for a given mass ratio and mass splitting.
(The corresponding param cards must have been created first, see 'madGraphParamCardGenerator.py')
-it calls upon lxbatch_madGraphJob.sh
The output will be tarred up and ready to copy over to soolins /storage for safe keeping
from soolin...
$scp taylor@lxplus.cern.ch:/afs/cern.ch/work/t/taylor/public/madGraphProjects/nmssmCascadeAnalysis_v01/*.tar.gz .
tar -xvf *.tar.gz (wildcard might not work so well…python script…bash line)
rm *.tar.gz

htcondor_executable.job info.
NEEDS FINISHING...to do with batch jobs on dice
