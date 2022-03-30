# factories

This directory contains a number of files, each of which contains code to run a certain procedure or create a specific object. In principle, setting it up like this should make
it easier to add new features. For example, if you have a new scenario, then all that would need to be added to the factory would be an `if` statement for the new scenario which
would load the corresponding scenario file from the scenario directory. Similarly, adding new fields to the fieldset would just require creating a new variable to the 
`create_fieldset()` function. 

The directory has the following files:
- `analysys_factory.py`: This file contains a class for running post-processing of the parcels output. The function `AnalysisFactory.create_procedure()` has various 
boolean variables, which trigger running the corresponding post-processing procedure.
- `fieldset_factory.py`: This is the class to create the fieldset for the parcels simulations. Depending on the boolean functions that are activated, the fieldset loads the 
field data files specified in `file_dict`.
- `pset_variable_factory.py`: This is the class to initialize the particle variables such as lon, lat, depth, etc. in the first year of a particle simulation (for subsequent years
we just load the variable data from the preceeding model year).
- `scenario_factory.py`: This class creates the model scenario that is used to run the parcels model (and which is also used in the analysis/visualization steps to load the file
names). 
- `visualization_factory.py`: This contains the class to run the visualization code. 
