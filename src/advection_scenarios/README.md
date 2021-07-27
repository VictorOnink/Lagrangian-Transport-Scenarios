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

