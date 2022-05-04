# Lagrangian Transport Scenarios
## Contents
1. [Introduction](#introduction)
2. [Citing the LTS framework](#citing)
3. [General overview](#overview)
4. [Usage guide](#usage)
5. [System requirements](#requirements)
6. [Acknowledgements](#acknowledgements)


## Introduction <a name="introduction"></a>
The Lagrangian Transport Scenarios (LTS) repository contains the framework to run, process and visualise oceanic Lagrangian particle-tracking simulations. The LTS framework was developed to standardize and simplify running a wide array of Lagrangian ocean scenarios, where each scenario contains a set of kernels that describe the behavior of particles which are advected with a specified set of oceanographic and meteorological data.

Almost all the code within this repository is written in python, with bash scripts being used to submit and run jobs. All the Lagrangian simulation code uses the Parcels (**P**robably **A** **R**eally **C**omputationally **E**fficient **L**agrangian **S**imulator) package, which is 
described in detail in [Lange & van Sebille (2017)](https://doi.org/10.5194/gmd-10-4175-2017) and [Delandmeter & van Sebille (2019)](https://doi.org/10.5194/gmd-12-3571-2019). For explanations of Parcels, please refer to the the [OceanParcels](http://oceanparcels.org/) website, which has a number of tutorials and references to various use-cases of Parcels.

## Citing the LTS framework <a name="citing"></a>
The LTS framework is largely based on the work for the [Onink et al. (2021)](https://doi.org/10.1088/1748-9326/abecbd) paper, which described stochastic plastic beaching and resuspension parametrizations within large-scale ocean simulations. However, various other scenarios and kernels have been described in other publications outlined below.

- The stochastic beaching and resuspension kernels within the `Stochastic`  and `ShoreDependentResuspension` scenarios are described in [Onink et al. (2021)](https://doi.org/10.1088/1748-9326/abecbd).
- The KPP vertical turbulent mixing kernel is described in [Onink et al. (2022)](https://doi.org/10.5194/gmd-15-1995-2022), with the model code to simulate 1D vertical transport of buoyant particles being available [here](https://github.com/VictorOnink/Wind-Mixing-Diffusion).
- The `SizeTransport` and `FragmentationKaandorpPartial` scenarios are described in Onink et al. (in prep) "The influence of particle size and fragmentation on large-scale microplastic transport in the Mediterranean Sea"

For any model components described above, please cite the relevant paper. For citing the LTS framework itself, please cite: *zenodo link*.

## General overview <a name="overview"></a>
The LTS framework is set up for three basic tasks: run the Lagrangian simulations, analyze the output and visualize the post-processed Lagrangian simulations. For each of these tasks, there is a dedicated submission bash script in the `bin` directory which is used to submit jobs onto the server. Currently, all these script are written for the [ubelix](https://ubelix.unibe.ch/) HPC cluster at the University of Bern, which uses a version of the Slurm workload manager. In principle the LTS framework can be applied to run servers that use different workload managers, but this would require modification of the submission scripts.

To speed up computation, the particle tracking simulations are parallelised into `runs`, where each run contains a subselection of the total number of particles in the simulation. This subselection is allocated at the beginning of the simulation according to the starting locations of the particles and remain fixed throughout the entire duration of the simulation, and the number of runs is dependent on the input scenario. Each run is also further subdivided into according to the `restart` variable. All multi-year simulations are subdivided into 1-year sections, where `restart = 0` indicates that this is a new simulation starting from the initial input files. The second simulation year would then be indicated by `restart = 1`, where each subsequent simulation year only starts when the simulation for the previous simulation year is complete (e.g. `restart = 2` only starts when `restart = 1` ends, which in turn only started once `restart = 0` is complete). Breaking up the simulation into these subsimulations according to the `restart` and `run` variables speeds up the computational time to run these simulations, and implies that an error in a later simulation year doesn't require rerunning the entire simulation.

The `bin/simulation_submission.sh` script also allows the user to set a number of variables, such as beaching timescales, fragmentation rates, ensemble members, etc. These are saved in the local `.env` file, which are then loaded into the `src/settings.py` file. This governs the scenario parameters, and the submission script is set up for easily submit a large number of jobs for e.g. a range of beaching timescale values. 

The `bin/analysis_submission.sh` is used to submit analysis jobs that take the unprocessed Parcels output and produce post-processing output files. A number of standardized analysis functions are built into the LTS framework, such as for calculating annual mean horizontal and vertical concentration distributions, calculating bayesian statistics and simple timeseries of beaching/adrift particle fractions. For almost all these analysis functions, the analysis consists of two steps. For the first step, a job is submitted that carries out the specified analysis procedure for each individual `restart`/`run` file and produces an intermediate output file. The second step then takes all these intermediate output files and compiles them into one output file which contains the analysis results for the simulation as a whole. For model scenarios with a N `restart`/`run` files, this two step procedure can increase the speed of the analysis by a factor ~N. Almost all analysis functions are standardized to work with any basic scenario, but a number are scenario specific (e.g. `src/Analysis/parcels_to_sizespectrum.py`).

Finally, the `bin/visualization_submission.sh` script is used to generate figures and animations of the post-processed data from the analysis stage. A number of standard functions for e.g. creating figure with a (M, N) subplot array or cartopy maps are provided, but in general each scenario will have its own set of code to generate figures. 

## Usage guide <a name="usage"></a>
In general, any general questions about the LTS framework (for either understanding it or applying it to another setting) can be posted in the [Issues](https://github.com/VictorOnink/Lagrangian-Transport-Scenarios/issues) section of this repository. However, a few basic use cases are described below.
<details> 
<summary> Running simulations backwards in time
</summary>
<p>
  
The LTS framework is set up to run both forward-in-time and backward-in-time simulations for any model scenario. Within the bash scripts, the `BACKWARD` variable specifies if a simulation is forward or backward, where `BACKWARD = 0` indicates forward in time and `BACKWARD = 1` indicates backwards in time.
</p>
</details>

<details> 
<summary> Running LTS on a different server
</summary>
<p>
  
In order to run LTS on a new server, it is first important to modify the job submission bash scripts in the `bin` directory to state `SERVER='NewServerName'`. Depending on the exact setup of this new server, it is likely also necessary to modify the job submission scripts to account for a different workload manager. In case the new server uses Slurm, these modifications might be minor (e.g. updating the partition and qos names to whatever names are used on the server). For other workload managers, this could require significantly more effort, and if necessary perhaps ask your local systems administrator for assistance.

Once the bash scripts are done, go to `src/settings.py` and modify the `DATA_DIREC` variable to include the path to the model data on `'NewServerName'`. This then automatically updates other directory paths such as `DATA_INPUT_DIREC`, `DATA_OUTPUT_DIREC` and `FIGURE_OUTPUT_DIREC`. However, `SCRATCH_DIREC` will need to be updated manually since the scratch directory is likely not a subdirectory of `DATA_DIREC`.

In principle, this is all that is required to set up LTS to run on a new server. Just make sure that all the necessary oceanographic/meteorological data is available on the server following the paths set in `src/advection_scenarios/advection_files.py` and you should be good to go.
  
</p>
</details>

<details>
<summary> Introducing new model variables
</summary>
<p>
  
In order to create a new variable that can be accessed throughout the LTS framework, you simply need to add it to `src/settings.py`, where the convention is to have all `src/settings.py` variables be completely capitalised. If the variable is e.g. a physical constant such as the density of air or some other quantity that doesn't generally vary between different scenarios, then it is easiest just to give it a set value within the settings file. However, if it is a model parameter (e.g. beaching timescale, fragmentation timescale, etc.), then it is preferable to set the value using the `load_env_variable()` function. This loads the variable value from the `.env` file, where the value can be set within the bash submission scripts. 
  
</p>
</details>

<details>
<summary> Creating a new model scenario
</summary>
<p>
  
In order to create a new scenario, a number of steps need to be taken:
1. Within the bash scripts, set `SCENARIO='NewScenario'`. The `src/settings.py` file will load the `SCENARIO` variable factor, and then know that you are running a `NewScenario` simulation. If you need to specify new model variables within the bash scripts, then do so as well.
2. Not much needs to happen within `src/settings.py`, except to update it for any new scenario variables you defined within the bash scripts before and to update the log section at the end of file to print the particular model variables relevant for `NewScenario`
3. Within `src/scenarios`, you need to create a new scenario file for `NewScenario`. Each scenario has its own python Class, which contains all the necessary functions to run a new simulation. This scenario class is actually a derived class from `BaseScenario`, which is defined in `src/scenarios/base_scenario.py`. As such, there are a number of mandatory functions that need to be defined in the `NewScenario` class in order for the scenario to work, such as functions defining the particle behavior and the names of the output files. Please refer to `src/scenarios/example_scenario.py` or any of the other scenario files for examples in how to set up a scenario. The `BaseScenario` class contains methods such as `BaseScenario.run()` to actually carry out the final model run, but if you need a modified form of these functions then you can define them again in the derived class to overwrite these functions in the base class.
4. Once you have created the scenario file, you must import the scenario class in `src/factories/scenario_factory.py` and add the new scenario to the `ScenarioFactory.create_scenario()` method.
Once these steps have been taken, everything should be ready to go to run any `NewScenario` simulation. In principle, the analysis functions should be ready as well, although it might be necessary to modify some of the functions if you have any model parameters that you want to loop through aside from just `run` or `restart` (see `src/Analysis/parcels_to_concentration.py` as an example of how changes were made for the `SizeTransport` and `FragmentationKaandorpPartial` scenarios.) If you want to also build in visualization capabilities, then there is an additional step:
6. Add a `NewScenario` directory to `src/visualization`, within which you create `NewScenario_Figures.py`. This file will then have a function `run()`, which will contain all the code for creating figures. 
7. Go to `src/factories/visualization_factory.py` and add `visualization.NewScenario.NewScenario_Figures.run()` to the `VisualizationFactory.run()` method. This will finalize visualization capabilities for `NewScenario`.
  
</p>
</details>

<details>
<summary> Creating a new advection scenario
</summary>
<p>

The advection scenario refers to the set of data including the ocean currents, temperature/salinity, mixed layer depth, wind, etc. that is used to drive the particle transport simulation. Given that the ocean currents are generally the main driver of particle transport, the advection scenario is generally named after the ocean current data (e.g. `CMEMS_MEDITERRANEAN` refers to the advection scenario where the ocean currents are from the [CMEMS Mediterranean Sea Physics Reanalysis](https://doi.org/10.25423/CMCC/MEDSEA_MULTIYEAR_PHY_006_004_E3R1)).

Within the submission scripts and `src/settings.py`, the only steps that need to be taken is that the `ADVECTION_DATA` variables within the bash scripts must be updated to the new advection scenario name `NewAdvection`. This is then automatically set within the `src/settings.py` file.

The `FieldSet` object that contains all the field data paths within the Parcels simulations is created with the `FieldSetFactory.create_fieldset()` method defined in `src/factories/fieldset_factory.py`.  Within the specific scenario class, the `FieldSetFactory.create_fieldset()` method is called, where the variables specify what data is supposed to be loaded (e.g. `FieldSetFactory.create_fieldset(wind=True)` will create a `FieldSet` object containing wind data). Now, the `FieldSetFactory.create_fieldset()` knows the paths to all the necessary data by loading the `file_dict` variable, which is defined within `src/advection_scenarios/advection_files.py`. 

`src/advection_scenarios/advection_files.py` contains the `AdvectionFiles` class, where the `AdvectionFiles.file_names()` method returns the `file_dict` dictionary that contains all the paths to the data necessary to create the `FieldSet` object and run the simulation. Each advection scenario has its own set of paths, and for all the variables the `file_dict` contains the file names, the variable names and the data dimensions. Please check the other advection scenarios for naming conventions, as otherwise the `FieldSetFactory.create_fieldset()` will not be able to read in the data. It is not necessary that each advection scenario specifies paths for all types of oceanographic data, but you do have to make sure you specify paths for all the data you call with the `FieldSetFactory.create_fieldset()` method (e.g. if you run `FieldSetFactory.create_fieldset(wind=True, MLD=False)`, you need to have paths specified for the wind data, but not for the mixed layer depth data).

Aside from file paths for oceanographic data specific to an advection scenario, there are also a number of input files that need to be created in order for e.g. beaching and diffusion kernels to run. In principle code is provided to automatically calculate distance to shore, land ID and horizontal grid size files. However, all these codes were written under the assumption that the ocean current data for the zonal and meridional currents are provided on an A grid (so no staggered grids). In order for the code to run on C grids, the following files would need to be updated:
- `src/advection_scenarios/create_boundary_current.py`
- `src/advection_scenarios/create_distance_to_shore.py`
- `src/advection_scenarios/create_grid_spacing.py`
- `src/advection_scenarios/create_land_ID.py`
- `src/advection_scenarios/create_tidal_Kz_files.py`, where this is only necessary if you are running simulations including vertical tidal mixing.

</p>
</details>

<details>
<summary> Creating a new input scenario
</summary>
<p>
  
A number of input scenarios are currently specified within the LTS framework:
- `Jambeck`: This weighs inputs according to mismanaged waste estimates from [Jambeck et al. (2015)](https://doi.org/10.1126/science.1260352) and coastal [population densities](https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals-rev11). A fixed fraction of all estimated mismanaged waste within 50 km of the ocean is estimated to enter the ocean, and the input distribution and particle weights are scaled according to the estimated waste input
- `Lebreton`, `LebretonDivision` and `LebretonKaandorpInit` are based on riverine plastic inputs from [Lebreton et al. (2017)](https://doi.org/10.1038/ncomms15611), where the Lebreton inputs are used for the spatial distribution and/or the weights of the particles
- `Point_Release`: Releases `PARTICLE_NUMBER` number of particles at a single point, where this point is indicated by the `RELEASE_SITE` variable within the submission scripts and the `SITE_LONLAT` within `src/settings.py`. 
- `Uniform`: Releases particles throughout the model domain with a fixed `RELEASE_GRID` degree spacing. 
To create the input files for each of these input scenarios, the `create_input_files.create_files()` method is run within `src/advection_scenarios/advection_files.py`. This method assures that all the inputs are placed on ocean cells instead of on land. However, all this code is again written assuming that the zonal and meridional currents are provided on an A grid. In order to run simulations with C grid data, `create_input_files.create_files()` would need to be modified.

In order to create an entirely new input scenario, the `INPUT` variable within the submission scripts need to be set to `INPUT=NewInput`. Then, `INPUT_DIREC_DICT` within `src/settings.py` needs to be updated with the path to the directory that will contain the input files.

Subsequently, if the new input scenario scales the particle release according to input estimates such as riverine inputs, the `INPUT_MAX` and `INPUT_MIN` variables within `src/settings.py` specify the maximum and minimum weight a single parcels particle represents. For example, if  `INPUT_MAX = 1` and `INPUT_MIN = 0.2` and we have a site with 2.1 input units, then we'd release 2 particles at this site with weight 1, whereas the remaining 0.1 input units is not large enough to warrant creating another particle. This means that some inputs are neglected, but this can dramatically reduce the number of particles (and thus computational cost) for a simulation. By setting the `INPUT_MAX` and `INPUT_MIN` parameters the user can determine which sources to include or not. Then, `INPUT_DIV` sets the number of particles within each `run`, where e.g. a simulation with 25,000 particles and `INPUT_DIV = 10 000` would mean the simulation is broken up into three runs with 10 000, 10 000 and 5 000 particles. Again, it is up to the discretion of the user to determine a suitable  `INPUT_DIV` value. 

Then, `create_input_files.create_files()` needs to be updated to account for the new input scenario. For this, we refer you to the code for the other input scenarios for examples on how to set up this code.

Once the number of runs has been determined, it is important to update `RUNRANGE` in the submission scripts and `src/settings.py`, such that the code knows how many runs to submit jobs for, and how many output files need to be considered in the final analysis. 
  
</p>
</details>

<details>
<summary> Creating a new analysis function
</summary>
<p>
  
In order to create a new general analysis function, follow the subsequent steps:
1. Add a variable `NEWANALYSIS` to the `bin/analysis_submission.sh` script, where `NEWANALYSIS = 0` indicates not to run this analysis step and `NEWANALYSIS = 1` does indicate to run it.
2. Add `NEWANALYSIS` to `src/settings.py`, which loads the `NEWANALYSIS` variable from the `.env` file.
3. Add a file for the analysis procedure to the `src/Analysis` directory. The typical naming convention is `parcels_to_NEWANALYSIS.py`, which will give a general indication of what type of analysis is being conducted in this run. 
4. If possible, we suggest setting up the `parcels_to_NEWANALYSIS.run()` class so that it follows the two step parallelization procedure outlined earlier, as this will generally speed up the code. For examples on how to approach this, please refer to the other codes. In addition, we suggest making the code as general as possible so that it might be applied to other scenarios as well. However, this might not always be possible and if not, we suggest adding a statement such as `assert self.scenario_name in ['ScenarioName'], "The NEWANALYSIS function is not set up for {}".format(self.scenario_name)` to prevent the code from running for scenarios for which it is not intended.
5. Once the analysis code has been written, go to `src/analysis_factory.py` and add the new analysis code to the `AnalysisFactory.run_analysis_procedure()` method. 
  
</p>
</details>

## System requirements <a name="requirements"></a>
To install the necessary packages for the LTS framework, run `conda env create -f environment_LTS.yml`. This will install all the necessary packages and dependencies to run all the LTS code.

## Acknowledgements <a name="acknowledgements"></a>
The following people have been involved in the development of the LTS framework:
- [Victor Onink](https://github.com/VictorOnink) has done the full implementation and had a leading role in the development of all the scenarios within the LTS framework.
- [Charlotte Laufkötter](https://github.com/blauhai) has been involved in the conceptual and theoretical development of all the scenarios.
- [Inger van Boeijen](https://github.com/IngerMathilde) helped in the initial structuring of the LTS framework.
- [Erik van Sebille](https://github.com/erikvansebille) was involved in the conceptual and theoretical development of all the scenarios.
- [Christian Kehl](https://github.com/CKehl) and [Philippe Delandmeter](https://github.com/delandmeterp) contributed to the development of all Parcels components of the LTS framework. 
- [Cleo Jongedijk](https://github.com/cjongedijk) assisted in the conceptual development of the beaching and resuspension parametrizations in the `Stochastic` scenario.
- [Mikael Kaandorp](https://github.com/mikaelk) was involved in the development of the `FragmentationKaandorpPartial` scenario.
