# Next-Stage-Plastic-Beaching
Building on the model setup for [Onink et al. (2021)](https://doi.org/10.1088/1748-9326/abecbd),
this respository contains an updated model setup for running, analyzing and visualizing 
Lagrangian particle simulations. 

Almost all code is written in python and uses the parcels (**P**robably **A** 
**R**eally **C**omputationally **E**fficient **L**agrangian **S**imulator) package, which is 
described in detail in [Lange & van Sebille (2017)](https://doi.org/10.5194/gmd-10-4175-2017)
 and [Delandmeter & van Sebille (2019)](https://doi.org/10.5194/gmd-12-3571-2019). For 
detailed explanations regarding the working of parcels, I refer you to the parcels website
[here]([parcels](http://oceanparcels.org/)), which has a number of tutorials and plenty of
documentation/references of papers applying parcels in a variety of settings.

Each subdirectory within this repository contains a README file that explains the basic
function of all the files. More detailed comments regarding the inner working of the code
is within the python files themselves.

The directories contained here are:
- `bin`: This directory contains a number of bash scripts that allow submission of jobs on
the [ubelix](https://ubelix.unibe.ch/) HPC cluster at the University of Bern.
  
- `old_model_setup`: This directory contains an old model setup. It is included as a
reference, but its use is not encouraged.
  
- `src`: This contains all the model source code.

All code was written by [Victor Onink](https://github.com/VictorOnink). The structuring of 
the model code was done in close collaboration with 
[Inger van Boeijen](https://github.com/IngerMathilde).
