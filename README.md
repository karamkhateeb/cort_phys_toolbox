# README #

### What is this repository for? ###

This repository contains code for predicting photothrombotic lesion profiles based on illumination parameters along with code for analyzing the effect of illumination parameters on lesion sizes. It also contains code for analyzing neurophysiological analysis following lesioning and/or stimulation with our semi-transparent electrocorticographic array.

### Power Analysis ###
Code in the PowerAnalysis directory demonstrates the analysis of results from monkeys D, E, F, and G.

Important MATLAB Functions/Scripts:
- extract_all_data_function is used to extract signals saved using Ripple's Trellis software
- SaveSignals.m contains code for loading, downsampling, and filtering signals into frequency bands
- SignalPower is used to calculate signal power
- findArtifacts is used to find the time indices of signals for which the signal exceeds a specified threshold to identify motion/touch/etc. artifacts during recording
- VisBadChannels is used to determine which channels are bad
- PowerAnalysisAndFigs is the script used to conduct power analysis and generate figures for manuscript. To successfully run these scripts, will need to download boundedline.m and suplabel.m, which can be downloaded from MATLAB file exchange.
- StimAnalysisAndFigs is the script used to conduct power analysis on stimulation data and generate figures for manuscript. To successfully run these scripts, will need to download boundedline.m and suplabel.m, which can be downloaded from MATLAB file exchange.

### Histo Analysis ###
Code in the HistoAnalysis directory demonstrates the 3D reconstruction and volume estimation of photothrombotic lesions for monkeys B, C, D, and E. There is also a generic script that can be altered to perform histological analysis on any photothrombotic lesion (PT_3D_Reconstruct).

### Lesion Quantification Analysis ###
Code in this directory is used to determine the effects of light illumination parameters on photothrombotic lesion sizes as measured histologically or through OCTA imaging. It also contains data for the light intensities and aperture diameters tested.

### Light Simulation ###
The LightSim directory contains code for modeling light propagation through cortical tissue and predicting lesion profiles based on illumination parameters.

Note: Throughout the scripts contained in this repository, the following terms are used interchangeable:
- PT1 = Monkey A
- PT2 = Monkey B
- PT3 = Monkey C
- PT4 = Monkey D
- PT5 = Monkey E
- PT6 = Monkey F
- PT7 = Monkey G


* Repo owners and admins: Karam Khateeb (kkhateeb@uw.edu), Julien Bloch (julienb@uw.edu), Jasmine Zhou (jzhou33@uw.edu), Azadeh Yazdan-Shahmorad (azadehy@uw.edu)
