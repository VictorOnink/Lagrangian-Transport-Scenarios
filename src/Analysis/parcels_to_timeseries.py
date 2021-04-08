import settings as settings
import utils

from netCDF4 import Dataset
import numpy as np
from scipy import io
import progressbar
import os


def parcels_to_timeseries(file_dict: dict, lon_min: float = -180, lon_max: float = 180, lat_min: float = -90,
                          lat_max: float = 90):
    domain = [lon_min, lon_max, lat_min, lat_max]
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/'
    # Get the time axis
    for restart in range(settings.SIM_LENGTH):
        if restart == 0:
            parcels_file = file_dict[0][restart]
            dataset = Dataset(parcels_file)
            time_list = dataset.variables['time'][0, :-1]
        else:
            time_list = np.append(time_list, dataset.variables['time'][0, :-1])

    # Initializing the arrays of the timeseries
    beached_weight = np.zeros(time_list.shape)
    total_weight = np.zeros(time_list.shape)
    coastal_weight = np.zeros(time_list.shape)

    os.system('echo "Start running through the restart and run files"')
    # loop through the runs
    for run in progressbar.progressbar(range(settings.RUN_RANGE)):
        # Loop through the restart files
        for restart in range(settings.SIM_LENGTH):
            # Load the lon, lat, time, beach, distance and weight data
            parcels_file = file_dict[run][restart]
            dataset = Dataset(parcels_file)
            time = dataset.variables['time'][:, :-1]
            lon, lat = dataset.variables['lon'][:, :-1], dataset.variables['lat'][:, :-1]
            weight = dataset.variables['weights'][:, :-1] * settings.BUOYANT
            beach = dataset.variables['beach'][:, :-1]
            distance = dataset.variables['distance'][:, :-1]
            # Get the domain specific timeseries
            weight_sel, beach_sel, distance_sel, time_sel = utils._particles_in_domain(domain=domain, lon=lon, lat=lat,
                                                                                       weight=weight,
                                                                                       beach=beach, time=time,
                                                                                       distance=distance,
                                                                                       return_sub_array=True)
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
                    if weight_total != np.nan:
                        total_weight[t] += weight_beached
                    # and the total amount of plastic that is afloat in the coastal/beaching zone
                    weight_coastal = np.nansum(weight_sel[(time_sel == t_value) & (beach_sel == 0) & (distance_sel <= settings.COAST_D)])
                    if weight_coastal != np.nan:
                        coastal_weight[t] += weight_coastal
    # Get the output dictionary
    output_dict = {'beached': beached_weight, 'coastal': coastal_weight, 'total': total_weight}
    # Saving the computed concentration, where we have a few default names for the prefix
    if lon_min == -180 and lon_max == 180 and lat_min == -90 and lat_max == 90:
        prefix = 'timeseries-global'
    else:
        prefix = 'timeseries-lon_min={}-lon_max={}-lat_min={}-lat_max={}'.format(lon_min, lon_max, lat_min, lat_max)

    output_name = output_direc + utils._analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
    io.savemat(output_name, output_dict)
    os.system('echo "The timeseries has been saved"')
