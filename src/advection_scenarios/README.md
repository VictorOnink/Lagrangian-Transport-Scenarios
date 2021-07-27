# advection_scenarios
This directory contains code related to running the model with different ocean reanalysis data. This involves both
defining the file paths for the data that is loaded into the `fieldset` as well as creating netCDF files containing
fields such as anti-beaching currents, a land ID field or a field defining the distance to the nearest model coastline
for each ocean gridcell. Everything is written as generally as possible such that adding a new advection scenario should
be relatively straightforward, but additional debugging would likely be necessary.

## Main data sources
### Ocean circulation data
Currently, the model is set up to run with three overarching sets of reanalysis data for the ocean circulation data, and
these are the ones we use to define the advection scenario.
- `HYCOM_GLOBAL`: This uses the GOFS 3.1 41-layer HYCOM NCODA GLOBAL 1/12° Reanalysis as described [here](https://www.hycom.org/dataserver/gofs-3pt1/reanalysis)
- `HYCOM_CARIBBEAN`: This uses the HYCOM + NCODA Gulf of Mexico 1/25° Reanalysis described [here](https://www.hycom.org/data/gomu0pt04/expt-50pt1)
- `CMEMS_MEDITERRANEAN`: This uses the Mediterranean Sea Physics Reanalysis from CMEMS described [here](https://doi.org/10.25423/CMCC/MEDSEA_MULTIYEAR_PHY_006_004_E3R1)

### Wave Data
These are the main sources of wave related data (e.g. Stokes drift, wave periods)
- `CMEMS_MEDITERRANEAN`: For this advection scenario, we use the CMEMS sea wave reanalysis described [here](https://doi.org/10.25423/CMCC/MEDSEA_HINDCAST_WAV_006_012)
- `HYCOM_GLOBAL` and `HYCOM_CARIBBEAN`: For these scenarios we use the [WaveWatch III reanalysis](https://polar.ncep.noaa.gov/waves/hindcasts/)

### Wind data
For the wind data, all advection scenarios currently use data from the ERA5 reanalysis, which can be found [here](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5)

## Code structure
The advection scenarios are defined in `advection_files.py`, within the `AdvectionFiles` class. Running the
`AdvectionFiles.file_names()` function (which is called at the beginning of each scenario) creates the `file_dict`
dictionary, which contains all the paths for the files necessary to create the fieldset object for a given scenario.
Meanwhile, there are also a number of files which are created specific to an advection scenario, such as a file containing
the border current fields. When filling the `file_dict` dictionary, the code checks if these files already exists. In
case they do not, the necessary code is run to create the file (but this is only if it doesn't exist, as creating some
of these files can take a while).

Aside from `advection_files.py`, this directory contains the following files:
- `create_boundary_current.py`: This file is used to create the anti-beaching boundary current that prevents particles
from getting stuck at land/ocean boundaries.
- `create_distance_to_shore.py`: This file computes the distance to the nearest model coastline, either for land cells
(used for creating the input files) or for ocean cells (used to determine whether particles are in the coastal beaching
zone).
- `create_grid_spacing.py`: This field calculates the horizontal grid spacing.
- `create_land_id.py`: This creates a boolean grid, where 0 indicates ocean and 1 indicates land cells.
- `create_MLD_files.py`: This uses temperature and salinity data to compute the MLD according to [de Boyer Montegut et al. (2004)](https://doi.org/10.1029/2004JC002378).
- `create_tidal_Kz_files.py`: This computes the vertical Kz diffusion coefficient using data from [de Lavergne et al. (2020)](https://doi.org/10.1029/2020MS002065).

## Input scenarios
Finally, there are the input scenarios, where the files containing initial particle longitudes, latitudes and weights
are created with the `create_input_files.py` file. Currently, we have the following input scenarios possible:
- `Jambeck`: This uses the nationwide estimates of mismanaged plastic waste per capita from [Jambeck et al. (2015)](https://doi.org/10.1126/science.1260352)
together with population density estimates for 2010.
- `Lebreton`: This uses estimates of annual riverine plastic inputs from [Lebreton et al. (2017)](https://doi.org/10.1038/ncomms15611)
- `Point_Release`: This releases a particle at a single point, where the initial coordinates are set in `src/settings.py`
- `Uniform`: This releases particles on a uniformly spaced lon/lat grid over the entire basin (with the exception of
land).