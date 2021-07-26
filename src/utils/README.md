# utils

This directory contains various files containing a range of utility functions, where a utility function is just any
function I find is useful to use throughout the code. This ranges from simple functions to check if a file exists
to Parcels kernels for standard processes such as advection and diffusion.

Contained within this directory are the following files:
- `analysis_utils.py`: This contains a number of standardized functions (such as a histogram function) used in the post-
processing of the Parcels output file.
- `BaseParticle.py`: This contains the baseparticle class used as a basis for creating all other particle classes in the
parcels simulations
- `file_utils.py`: This contains functions for basic manipulation of files, such as checking if a file exists, setting
directory paths and saving/loading pickle objects.
- `physics_utils.py`: This contains parcels kernels representing various physical processes, such as for computing
advective and diffusive transport.
- `run_utils.py`: This contains functions for running the parcels simulations, such as adding a variable to a
particle class or setting the seed value for the random number generator.