Information on running batch code.

lxbatch_submitMadGraphJobs.py info.
Use to submit Madgraph jobs onto lxplus.
Will scan over mass of higgs and mass susy for a given mass ratio and mass splitting.
(The corresponding param cards must have been created first, see 'madGraphParamCardGenerator.py')
It calls upon lxbatch_madGraphJob.sh
The output will be tarred up and ready to copy over to soolins /storage for safe keeping
from the correct soolin storage location...
$scp taylor@lxplus.cern.ch:/afs/cern.ch/work/t/taylor/public/madGraphProjects/nmssmCascadeAnalysis_v01/*.tar.gz .

htcondor_executable.job info.
A htcondor .job template, used by htcondorExeJob.py
htcondorExeJob.py info.
A class for submitting jobs to DICE, used by submit_htcondorExeJob.py
submit_htcondorExeJob.py info.
Choose an executable and input and output files. It will run them on DICE (as long as the cpp is set up for it)
Will need to hadd the outputs
$ hadd combinedHistos.root *.root
Then run the plotter on this output.