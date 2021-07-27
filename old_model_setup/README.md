# old_model_setup
This contains the old model_setup prior to the restructuring, and it is essentially the same
code as was used for [Onink et al. (2021)](https://doi.org/10.1088/1748-9326/abecbd). However,
I would not recommend using the current code, as some parts are outdated, and others are now
set up much more efficiently in the new setup within the `src` directory. 

The files within are:
- `DataSetup.py`: this contains all code for creating the model fieldset
- `ModelMasterFile.py`: This loads the parameters from the bash script and runs the actual
simulation (note: the bash script for this model setup is an outdated version of 
  `bin/simulation_submission.sh` and is not provided in this repository)
  
- `ParticleBeaching.py`: This contains all the beaching kernels
- `ParticleClasses.py`: This contains the various types of particle classes
- `ParticleSetup.py`: This is used to create the particle sets, whether that involves loading
initial input arrays or restarting from a previous file.
  
- `ParticleTransport.py`: This contains all the parcels kernels for particle transport.

Again, this code is a very old setup and its use beyond acting as a reference is not
recommended.
