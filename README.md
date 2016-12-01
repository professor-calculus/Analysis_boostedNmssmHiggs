Information on how to run the tools in this repository.
This repo is for the nmssm cascade analysis.

Setup within a CMSSW project.
$ cd CMSSW_X_Y_Z/src
$ cmsenv
$ git clone git@github.com:joseph-taylor/Analysis_boostedNmssmHiggs.git Analysis/Analysis_boostedNmssmHiggs
$ scram b -j8

crab3 directory:
Creating/Processing data
see crab3/crab3_README.md

bin directory:
All the cpp code
see bin/bin_README.md

batch directory:
How to use batch programs.


macros:
To be run interactively with ROOT.
$ root -l -b -q macro.cxx

massParamsToKinematics.cxx info.
Makes many plots to do with the kinematics of the decay chain by just doing calculations, not simulation.
Also produces .tex files which will automatically make pdfs with all the plots when run in latex locally.
There is a simple 'USER INPUT' section, and that is all you need to change.

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

untarMadgraphProject.py info.
use to untar all the projects when copied over from lxplus.

unzipMadgraphLhe.py info.
use to unzip the events.lhe file in the madgraph project.