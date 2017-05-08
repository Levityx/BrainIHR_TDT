# BrainIHR_TDT

## Description 

This Matlab package was developped to interface Matlab with the IHR's TDT system, 
in order to record single units with 128 channels. 
Since this is optimised for the lab's hardware, it will not work as is in any other place 
and is provided as a model. 

To start -assuming the package +IHR_TDT is accessible- run

```>> guideSetupSystem3```

which will launch a GUI that manages AlbanTDT and SweepTDT objects in the background. 
Calibrating and recording data (up to 128 electrodes) will generates the pictures found in the folder pics/. In this folder, the report was automatically generated, and the tdt.mat file contains a variable

### Author

Alban, June 2016
