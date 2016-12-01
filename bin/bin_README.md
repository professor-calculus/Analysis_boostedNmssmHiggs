Information on how to run the cpp.

Compile everything with $scram b -j8
This then gives us executables in $CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin/...
You should add this address to your path with
$ export PATH=$CMSSW_BASE/tmp/slc6_amd64_gcc530/src/Analysis/Analysis_boostedNmssmHiggs/bin:$PATH

run with
$ executableName options
Note:
-it runs with respect to where you are when you execute
-the executable is just the script name without the .cpp extension
-for the options see the specific script



All the code utilises the FWLite command parser.
@Mandatory Options:

-Input of single files
inputfiles=file.root
-Input of multiple files
inputfiles=file1.root,file2.root,file3.root (needs to be comma seperated without a space)

outputfile=outputfile.root
*This should be the full path and include a new directory (which the code creates) where the file will go.
The directory should be clear and unique.
eg.
~/local_Analysis_boostedNmssmHiggs/output_DoubleBTaggerEfficiencyStudies/CMSSW_8_0_20/GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph_boostedHiggsMC_25e777/histos.root

@Other options (don't need to explicitly write them):

maxevents=-1 can set the max events (-1 is all of them)
outputevery=1000 how often we print to screen
orderedsecondaryfiles=0 (default)
orderedsecondaryfiles=1 this does the plotting part only
ie.
you have the histograms in the root file already and we run on them to remake the plots, for example changing the aesthetics
it runs on the .root file stated by outputfile=outputfile.root (so in this case this is the input) 

eg.1 STANDARD
$ DoubleBTaggerEfficiencyStudies inputfiles=XYZ_1.root,XYZ_2.root outputfile=ABC/out.root

eg.2 REMAKE THE PLOTS
$ DoubleBTaggerEfficiencyStudies outputfile=ABC/histos.root orderedsecondaryfiles=1

NOTE:
This code makes a copy of itself and puts it in the output directory, but if you have made changes after compilation this will be out of phase, watch out!
It also makes a copy of the parser values in a seperate text file
Together, along with non-lazy directory name, you should have decent records of what you have done:)



To run on DICE
write 'runOnDice' as your last argument
It uses different logic in the scripts so that they can be used on DICE.
use batch/submit_htcondorExeJob.py to submit the jobs.




*Code Descriptions
'bin/DoubleBTaggerEfficiencyStudies.cpp'
To run on a boosted higgs mc samples. It looks at the effiency and properties of the boosted h->bb tagger.
it uses the class 'interface/PlottingDoubleBTaggerEfficiencyStudies.h' to create the plots