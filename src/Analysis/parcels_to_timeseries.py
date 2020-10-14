import settings as settings
import utils

from netCDF4 import Dataset
import numpy as np
from scipy import io
import progressbar
import os


def parcels_to_timeseries(file_dict: dict, lon_min: int = 26.7, lon_max: int = 42.4, lat_min: int = 41,
                          lat_max: int = 47):
    domain = [lon_min, lon_max, lat_min, lat_max]
    output_direc = utils._get_output_directory(server=settings.SERVER) + 'timeseries/'
    # Get the time axis
    for restart in range(settings.SIM_LENGTH):
        if restart == 0:
            parcels_file = file_dict[0][restart]
            dataset = Dataset(parcels_file)
            time_list = dataset.variables['time'][0, :-1]
        else:
            time_list = np.append(time_list, dataset.variables['time'][0, :-1])
    beached_weight = np.zeros(time_list.shape)
    total_weight = np.zeros(time_list.shape)
    os.system('echo "Start running through the loops of the files')
    # loop through the runs
    for run in progressbar.progressbar(range(settings.RUN_RANGE)):
        # Loop through the restart files
        for restart in range(settings.SIM_LENGTH):
            # Load the lon, lat, time, beach and weight data
            parcels_file = file_dict[run][restart]
            dataset = Dataset(parcels_file)
            time = dataset.variables['time'][:, :-1]
            lon, lat = dataset.variables['lon'][:, :-1], dataset.variables['lat'][:, :-1]
            weight = dataset.variables['weights'][:, :-1] * settings.BUOYANT
            beach = dataset.variables['beach'][:, :-1]
            # Get the domain specific timeseries
            weight_sel, beach_sel, time_sel = _particles_in_subarea(lon=lon, lat=lat, weight=weight, beach=beach,
                                                                    time=time, domain=domain)
            # Loop through all time_list values, which cover the entire length of the simulation and not solely this
            # specific year saved in the file
            for t, t_value in enumerate(time_list):
                # If any of the times of the arrays match that of the time_list, then...
                if np.nansum(np.array(time_sel == t_value) * 1) > 0:
                    # We calculate the total mass of particles that are beached
                    weight_beached = np.nansum(weight_sel[(time_sel == t_value) & (beach_sel == 1)])
                    if weight_beached != np.nan:
                        beached_weight[t] += weight_beached
                    # and the total mass of particles in the simulation at that time, excluding all particles where the
                    # beach indicator is 2, since that is only for particles that were previously deleted.
                    weight_total = np.nansum(weight_sel[(time_sel == t_value) & (beach_sel != 2)])
                    total_weight[t] +=
    # Get the output dictionary
    output_dict = {'beached': beached_weight, 'total': total_weight}
    # Saving the computed concentration, where we have a few default names for the prefix
    if lon_min == -180 and lon_max == 180 and lat_min == -90 and lat_max == 90:
        prefix = 'timeseries-global'
    else:
        prefix = 'timeseries-lon_min={}-lon_max={}-lat_min={}-lat_max={}'.format(lon_min, lon_max, lat_min, lat_max)

    output_name = output_direc + utils._analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
    io.savemat(output_name, output_dict)
    os.system('echo "The timeseries has been saved"')


def _particles_in_subarea(lon, lat, weight, beach, time, domain):
    lon_min, lon_max, lat_min, lat_max = domain
    # Select only particles that are within the domain that we are interested in
    select = (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)
    # And now we return the time, weight and beach arrays of the particles that are within
    # the domain of interest
    return weight[select], beach[select], time[select]
