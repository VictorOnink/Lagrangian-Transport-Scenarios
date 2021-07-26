# Scenarios

This directory contains all the scenario files for the model setup. Each scenario should be seen as a specific model setup. The `BaseScenario` within `base_scenario.py` specifies
the basic scenario building blocks. A number of functions are mandatory for any new scenario, such as the `file_names` function that sets the names of the output files. There are
also a number of functions that are the same for all scenarios, such as `run()` to actually execute the parcels simulation.

The basic mandatory components for any scenario are:
- `var_list`: This is a list containing the various variables that need to be added to the `parcels.ParticleSet`.
- `create_fieldset()`: This is a function that calls the fieldset factory and creates the fieldset object with all necessary data loaded in. The fieldset factory by default assumes
that data is **not** required, so to load specific data or create a certain fieldset constant, the boolean variable for that specific data must be set to True (e.g. 
setting `fieldset_factory.FieldSetFactory().create_fieldset(border_current=True)` would add the border current data to the fieldset)
- `get_pset()`: This function creates the `ParticleSet` object for the simulation.
- `get_pclass()`: This creates the particle class for the simulation (so where the variables tracked for each particle are defined).
- `file_names()`: This sets the names of the particle files.
- `beaching_kernel()`: This sets the beaching behavior at the coastal boundaries.
- `get_particle_behavior()`: This creates the particle `Kernel`, which defines the particle behavior and indicates which physical processes are taken into account.

There are also a number of functions which are the same across all scenarios (at the moment at least):
- `get_restart_variables()`: If the simulation is a restart from a preceeding simulation, this function loops through all the variables in `self.var_list` and loads these variables
from the old file.
- `get_var_dict()`: This function creates a dictionary that contains arrays to initialize the particleset, whether `restart == 0` or not.
- `run()`: This is really the main function that runs everything, loading the fieldset, creating the particleset and executing the parcels simulation. 
- `return_full_run_directory()`: This function is used more in the analysis, but it returns a dictionary that contains all the filenames within a simulation (so looped over all 
run and restart files)
