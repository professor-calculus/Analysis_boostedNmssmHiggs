#!/bin/bash
# run with $ bsub -q 8nh "sh lxbatch_madGraphJob.sh /path/To/madGraphInstructions.txt"
cd /afs/cern.ch/user/t/taylor/CMSSW_8_0_20/src/
eval `scramv1 runtime -sh`
cd /afs/cern.ch/user/t/taylor/MG5_aMC_v2_3_3/
./bin/mg5_aMC $1
tar -zcf $2.tar.gz $2 # tar it up at the end for easy transfer to soolin

# need to keep CMSSW version and MadGraph verison up-to-date