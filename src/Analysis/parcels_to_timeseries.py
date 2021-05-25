import settings as settings
import utils

from netCDF4 import Dataset
import numpy as np
from scipy import io
import progressbar
import os


def parcels_to_timeseries(file_dict: dict, lon_min: float = -180, lon_max: float = 180, lat_min: float = -90,
                          lat_max: float = 90):
    domain = lon_min, lon_max, lat_min, lat_max
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/{}/'.format(
        settings.SCENARIO_NAME)
    utils.check_direc_exist(output_direc)
    # Get the time axis
    time_list = np.array([], dtype=np.int)
    for restart in range(settings.SIM_LENGTH):
        parcels_file = file_dict[0][restart]
        parcels_dataset = Dataset(parcels_file)
        time_list = np.append(time_list, parcels_dataset.variables['time'][0, :-1])

    # Initializing the arrays of the timeseries
    beach_state_dict = {'beach': np.zeros(time_list.shape, dtype=float),
                        'afloat': np.zeros(time_list.shape, dtype=float),
                        'seabed': np.zeros(time_list.shape, dtype=float),
                        'removed': np.zeros(time_list.shape, dtype=float),
                        'total': np.zeros(time_list.shape, dtype=float),
                        'time': time_list}
    beach_label_dict = {'beach': 1, 'afloat': 0, 'seabed': 3, 'removed': 2}

    os.system('echo "Start running through the restart and run files"')
    # loop through the runs
    for run in range(settings.RUN_RANGE):
        # Loop through the restart files
        for restart in range(settings.SIM_LENGTH):
            # Load the lon, lat, time, beach and weight data
            parcels_file = file_dict[run][restart]
            parcels_dataset = Dataset(parcels_file)
            full_data_dict = {}
            for variable in ['lon', 'lat', 'beach', 'time']:
                full_data_dict[variable] = parcels_dataset.variables[variable][:, :-1]
            if 'weights' in parcels_dataset.variables.keys():
                full_data_dict['weights'] = parcels_dataset.variables['weights'][:, :-1] * settings.BUOYANT
            else:
                full_data_dict['weights'] = np.ones(full_data_dict['lon'].shape, dtype=float)

            # Just get the particles within the domain, which we do by setting all values not within the domain to nan
            # These will therefore not be taken into account in the calculations of total counts/weights
            within_domain = utils._particles_in_domain(domain=domain, lon=full_data_dict['lon'],
                                                       lat=full_data_dict['lat'])
            for variable in full_data_dict.keys():
                full_data_dict[variable][within_domain == False] = np.nan

            # Now, looping through all the time steps and adding up the total weight/counts of particles within each of
            # the beach state domains at each time step
            for index, time_value in enumerate(time_list):
                time_selection = full_data_dict['time'] == time_value
                if np.nansum(full_data_dict['time'] == time_value) > 0:
                    time_dict = {}
                    for variable in ['beach', 'weights']:
                        time_dict[variable] = full_data_dict[variable][time_selection]
                    for beach_state in beach_label_dict.keys():
                        beach_state_dict[beach_state][index] += np.nansum(time_dict['weights'][time_dict['beach'] == beach_label_dict[beach_state]])
                    beach_state_dict['total'][index] += np.nansum(full_data_dict['time'] == time_value)
    prefix = 'timeseries'

    output_name = output_direc + utils._analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
    utils.save_obj(output_name, beach_state_dict)
    os.system('echo "The timeseries has been saved"')
