# Analysis

This directory contains a number of files, each containing a standardized function for post-processing of the parcels
output.

The files (and corresponding functions) contained within this directory are:
- `parcels_to_basicstatistics.py`: Saves some basic statistics for a given simulation, such as max/min/mean/std values
of various particle characteristics.
- `parcels_to_concentrations.py`: Computes the horizontal concentration, both averaged over the entire length of the
simulation and for each simulation year.
- `parcels_to_concentration_monthly.py`: Computes the horizontal concentration, for each month in a simulation.
- `parcels_to_max_distance.py`: Computes the maximum distance each particle is removed from the model coastline throughout
a simulation. Note: this analysis should only be done if the particle simulation tracks the distance to shore throughout
the parcels simulation.
- `parcels_to_timeseries.py`: Computes the beached/floating/seabed fraction over time.
- `parcels_to_timeslices.py`: Splits up the parcels files into timeslices, where each file corresponds to all the
particles in the simulation at a given timepoint. This allows for easier loading of particles data in e.g. making
animations, as for large simulations loading all particle trajectory data is too much for the computer RAM.
- `parcels_to_vertical_concentration.py`: Computes vertical concentration profiles.
- `parcels_to_spatial_vertical_profiles.py`: Calculates vertical profiles for all particles within a certain gridded area.
- `parcels_to_separation_distance.py`: Computes distance between particles over time
- `parcels_to_particle_number.py`: Calculates the particle number for a given `LAMBDA_FRAG` value for parcels output
in the `KaandorpFragmentationPartial` scenario.
- `parcels_bayesian.py`: Calculates bayesian origin statistics for release clusters.
- `FragmentationKaandorp_boxmodel.py`: Contains code to run the box model of fragmentation and size dependent transport from Kaandorp et al. (2020).
- `parcels_to_lonlat_average.py`: Calculates lon lat averages of the particle concentrations, averaged over each simulation year.
- `CMEMS_mediterranean_mean_MLD.py`: Calculates monthly and seasonal mean MLD's for the `CMEMS_MEDITERRANEAN` data.
- `FragmentationKaandorp_boxmodel.py`: This contains an adaptation of the box model from  [Kaandorp et al. (2021)](https://doi.org/10.1088/1748-9326/abe9ea),
originally from [here](https://github.com/OceanParcels/ContinuousCascadingFragmentation/blob/main/box_model_Mediterranean.py).
- `FragmentationKaandorpPartial_mass_loss.py`: Calculates the mass loss due to microplastic fragmentation in the `FragmentationKaandorpPartial` scenario
- `FragmentationKaandorpPartial_Statistics.py`: Contains code to calculate some basic statistics for the `FragmentationKaandorpPartial` scenario
- `SizeTransport_Statistics.py`: Contains code to calculate some basic statistics for the `SizeTransport` scenario


We also include code for standardizing field data so that it can be easily plotted
- `FragmentationKaandorpPartial_fielddata.py`: standardizing field data of microplastic size distributions
