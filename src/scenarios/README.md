# Scenarios

This directory contains all the scenario files for the model setup. Each scenario should be seen as a specific model setup. The `BaseScenario` within `base_scenario.py` specifies
the basic scenario building blocks. A number of functions are mandatory for any new scenario, such as the `file_names` function that sets the names of the output files. There are
also a number of functions that are the same for all scenarios, such as `run()` to actually execute the parcels simulation.

The basic mandatory components for any scenario are:
- `set_prefix()`: Sets the prefix of the scenario to label the output files.
- `set_time_steps()`: This sets the integration, output and release periodicity timesteps
- `set_var_list()`: This sets a list containing the names of the various variables that need to be added to the `parcels.ParticleSet`.
- `create_fieldset()`: This is a function that calls the fieldset factory and creates the fieldset object with all necessary data loaded in. The fieldset factory by default assumes
that data is **not** required, so to load specific data or create a certain fieldset constant, the boolean variable for that specific data must be set to True (e.g. 
setting `fieldset_factory.FieldSetFactory().create_fieldset(border_current=True)` would add the border current data to the fieldset)
- `get_pset()`: This function creates the `ParticleSet` object for the simulation.
- `get_pclass()`: This creates the particle class for the simulation (so where the variables tracked for each particle are defined).
- `file_names()`: This sets the names of the particle files.
- `beaching_kernel()`: This sets the beaching behavior at the coastal boundaries.
- `get_particle_behavior()`: This creates the particle `Kernel`, which defines the particle behavior and indicates which physical processes are taken into account.

There are also a number of functions which are the same across all scenarios:
- `get_restart_variables()`: If the simulation is a restart from a preceeding simulation, this function loops through all the variables in `self.var_list` and loads these variables
from the old file.
- `get_var_dict()`: This function creates a dictionary that contains arrays to initialize the particleset, whether `restart == 0` or not.
- `run()`: This is really the main function that runs everything, loading the fieldset, creating the particleset and executing the parcels simulation. 
- `return_full_run_directory()`: This function is used more in the analysis, but it returns a dictionary that contains all the filenames within a simulation (so looped over all 
run and restart files)


Now, we have a number of different scenarios built in at the moment:
- `FragmentationKaandorp.py`: This is a fragmentation scenario, seeking to adapt the transition matrix microplastic fragmentation model from [Kaandorp et al. (2021)](https://doi.org/10.1088/1748-9326/abe9ea).
- 'FragmentationKaandorpPartial.py': Very similar to `FragmentationKaandorp.py`, but with a different approach to splitting particles.
- `SizeTransport.py`: This is a scenario for studying the effect of initial particle size on the 3D transport behavior of microplastic particles.
- `Turrel_Beaching_scenario.py`: This is a scenario where plastic resuspension is dependent on wind direction and water level. This is based on [Turrell (2018)](https://doi.org/10.1016/j.marpolbul.2018.10.024) and [Turrell (2020)](https://doi.org/10.1016/j.marpolbul.2020.111600), but hasn't been extensively tested.
- `advection_diffusion_only_scenario.py`: This is a beaching scenario where a particle beaches if it is advected onto a land cell. Based on the simple beaching implementations often used in the field, this has not been included in any final paper.
- `coastal_proximity.py`: This is a beaching scenario where a particle beaches when it is within the coastal zone for more than a predetermined period of time. This has not been included in any final paper
- `stochastic_scenario.py`: This is a beaching/resuspension scenario where beaching/resuspension are implemented stochastically when the particle is within the coastal zone. This is the main scenario described in [Onink et al. (2021)](https://doi.org/10.1088/1748-9326/abecbd).
- `shore_dependent_resuspension_scenario.py`: This beaching/resuspension scenario is very similar to the `stochastic_scenario.py` scenario, but here the resuspension rate is spatially varying based on the sandiness of the coastline. This scenario was also described in [Onink et al. (2021)](https://doi.org/10.1088/1748-9326/abecbd).
- `BlueCloudBackwards.py`: This is a 2D transport scenario set up to do back-in-time simulations for the Blue Cloud Hackathon
- `BlueCloudForwards.py`: The equivalent to `BlueCloudBackwards.py`, but then for forward-in-time simulations for the Blue Cloud Hackathon
- `BlueCloud.py`: This is a more general version of the `BlueCloudBackwards` and `BlueCloudForwards` scenarios set for for forward and backward in time simulations