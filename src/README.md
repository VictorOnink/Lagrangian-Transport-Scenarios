# src

This directory contains the source code for all components of the model except job submission
(files for this are found in `bin`).

Contained within this directory are the following files:

- `main.py`: This is essentially the master python file which actually runs the simulation, analysis and visualization
  jobs.

- `settings.py`: This is a settings file that sets a number of variables that are used throughout the model code. This
  includes settings file paths (e.g. for input data, parcels output, post-processing output and generated figures),
  loading in the run parameters set within the job submission bash files and a number of model parameters (e.g.
  gravitational acceleration, horizontal diffusion parameter, integration timesteps, etc.)

We also have a number of subdirectories:

- `advection_scenarios`: This contains code for setting the input data files that are used in running the parcels
  simulations.

- `Analysis`: This contains standardized post-processing code for analyzing the raw parcels output.

- `factories`: This contains a number of factories that are used to run procedures such as creating the fieldset,
carrying out post-processing and defining the simulation scenario.
  
- `scenarios`: This contains all the various model scenarios, where a scenario is a particular model setup. For example,
  there are separate scenarios for various types of plastic beaching 
  (see [Onink et al. (2021)](https://doi.org/10.1088/1748-9326/abecbd)), size dependent transport and fragmentation
  implementations, where each scenario differs in which physical processes are included. A new model setup would be set 
  up by creating a new scenario.
  
- `utils`: This contains a number of utility functions that are used throughout the model code. This ranges from 
  standardized parcels kernels (e.g. for wind mixing) to functions for removing a file.
  
- `visualization`: This contains all the code for creating figures visualizing the post-processed parcels output.