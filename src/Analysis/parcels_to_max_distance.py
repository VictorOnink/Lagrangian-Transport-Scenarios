import settings as settings
import utils

from netCDF4 import Dataset
import numpy as np
from scipy import io
import progressbar
import os


def parcels_to_max_distance(file_dict: dict, lon_min: float = -180, lon_max: float = 180, lat_min: float = -90,
                            lat_max: float = 90):
    domain = [lon_min, lon_max, lat_min, lat_max]
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/'

    # loop through the runs
    for run in progressbar.progressbar(range(settings.RUN_RANGE)):
        # Loop through the restart files
        for restart in range(settings.SIM_LENGTH):
            # Load the lon, lat, time, beach and weight data
            parcels_file = file_dict[run][restart]
            dataset = Dataset(parcels_file)
            lon, lat = dataset.variables['lon'][:, :-1], dataset.variables['lat'][:, :-1]
            distance = dataset.variables['distance'][:, :-1]
            # Get the maximum distance over the course of each trajectory
            max_distance_file = np.nanmax(distance, axis=1, keepdims=True)
            if restart == 0:
                # If this is the beginning of the run, we want to see if the initial positions of the particles fall
                # within the area we have defined as the domain
                within_domain = utils.particles_in_domain(domain=domain, lon=lon[:, :1], lat=lat[:, :1])
            # Making a cumulative array for the run
            if restart == 0:
                max_distance_file = max_distance_file[within_domain]
                max_distance_run = max_distance_file.reshape((len(max_distance_file), 1))
            else:
                max_distance_file = max_distance_file[within_domain]
                max_distance_run = np.hstack((max_distance_run, max_distance_file.reshape((len(max_distance_file), 1))))
        if run == 0:
            max_distance_simulation = np.nanmax(max_distance_run, axis=1, keepdims=True)
        else:
            max_distance_simulation = np.vstack(
                (max_distance_simulation, np.nanmax(max_distance_run, axis=1, keepdims=True)))
    # Get the output dictionary
    output_dict = {'max_dist': max_distance_simulation}
    # Saving the computed concentration, where we have a few default names for the prefix
    if lon_min == -180 and lon_max == 180 and lat_min == -90 and lat_max == 90:
        prefix = 'maximum_distance-global'
    else:
        prefix = 'maximum_distance-lon_min={}-lon_max={}-lat_min={}-lat_max={}'.format(lon_min, lon_max, lat_min,
                                                                                       lat_max)

    output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
    io.savemat(output_name, output_dict)
    utils.print_statement("The maximum distance has been saved")
