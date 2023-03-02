# Monopole Analysis

Structure for this repo:

[MonoAnalyzerPhoton/interface](https://github.com/sun51027/Monopole-Analysis/tree/main/MonoAnalyzerPhoton/interface)
:header file to be included in the analysis code.

[MonoAnalyzerPhoton/src](https://github.com/sun51027/Monopole-Analysis/tree/main/MonoAnalyzerPhoton/src)
:source code to be compiled and run the analysis

## MC (signal only) analysis

+ [MonoAnalyzerPhoton](https://github.com/sun51027/Monopole-Analysis/blob/main/MonoAnalyzerPhoton/src/MonoAnalyzerPhoton.cc)

This source file only analyze the signal MC. The input file are in Phat eos space:

`/eos/cms/store/group/offcomp_upgrade-sw/srimanob/monopole/13TeV/Legacy-RECO-v2/`

**note that you must change the path of inputFile in the code [line](https://github.com/sun51027/Monopole-Analysis/blob/main/MonoAnalyzerPhoton/src/MonoAnalyzerPhoton.cc#L393)**

before you run, create two directories:

`mkdir output
 cd output
 mkdir csv_file` 

### How to run

We use macro (root) to compile/run the analysis. 

`root -l 'MonoAnalyzerPhoton.cc("year","mass",(bool)matching_option,(int)analysis_type)' `

Four arguments are needed: (string)year, (string)mass, (bool)matching_option, (int)analysis_type.

If there is no neccessary to study the truth-matching and systematic uncertainty, the last two options use 0, e.g.

`root -l 'MonoAnalyzerPhoton.cc("year","mass",0,0)' `

year = 2018, 2017, 2016, 2016APV

mass = 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500

### output file:

There are two output files: 

1. Signaleff_*.csv : all cutflow counting and signal efficiency written in csv file. 

2. MonoPhotonAnalysis_*.root : all cutflow plots and N-1 plots written in ROOT file.

### other 

Don't care about the "warning" message.

Contact Lin if you have any questions: lshih@cern.ch



