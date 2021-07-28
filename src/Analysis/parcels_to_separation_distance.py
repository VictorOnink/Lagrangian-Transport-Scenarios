from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
from progressbar import ProgressBar
import settings
import utils


def parcels_to_separation_distance(file_dict: dict, scenario):
    """
    This code is to calculate the average separation distance over time between two sets of particles
    This is currently being written to work just for the SizeTransport scenario, but future versions might be able to
    adapt that to be more flexible
    :param file_dict: dictionary containing all file runs
    :param scenario: this is the scenario object, from which we can get the file_names function
    :return:
    """
    # Getting the directory saving the output files
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'separation_distance/{}/'.format(
        settings.SCENARIO_NAME)
    utils.check_direc_exist(output_direc)

    # Setting the particle sizes we loop through and creating a dictionary to save the average separation distance
    particle_size = np.array([1]) * settings.SIZE_FACTOR
    time_steps = len(Dataset(scenario.file_names(new=True, run=0, restart=0)).variables['time'][0, :])
    output_dict = dict.fromkeys(['STD', 'MEAN', 'MEDIAN'])
    for key in output_dict.keys():
        output_dict[key] = dict.fromkeys(range(time_steps))
        for time in output_dict[key].keys():
            output_dict[key][time]= {}
            for size in particle_size:
                output_dict[key][time][utils.init_size_key(size)] = 0

    # Setting the file name
    prefix = 'separation_distance'
    output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)

    # Starting to loop through the runs
    pbar = ProgressBar()
    for run in pbar(range(settings.RUN_RANGE)):
        # We only look at the first year of the simulation at the moment
        for restart in range(1):
            dataset = Dataset(file_dict[run][restart])
            lon_reference, lat_reference = dataset.variables['lon'][:], dataset.variables['lat'][:]
            for size in particle_size:
                comparison_dataset = Dataset(scenario.file_names(new=True, run=run, restart=restart, init_size=size))
                lon_comparison, lat_comparison = comparison_dataset.variables['lon'][:], comparison_dataset.variables['lat'][:]
                # Looping through all the time steps
                for time in range(1, time_steps):
                    key_size = utils.init_size_key(size)
                    distance = utils.distance_between_points(lon_reference[:, time], lat_reference[:, time],
                                                             lon_comparison[:, time], lat_comparison[:, time])

                    def median(array):
                        sorted_array = np.sort(array)
                        if time == 176:
                            print(len(array))
                            print(len(sorted_array))
                            print(array[len(array) // 2])
                            print(sorted_array[len(array) // 2])
                        return sorted_array[len(array) // 2]

                    print('{} {} {}'.format(time, np.mean(distance), median(distance)))
                    # if time == 176:
                    #     print(distance)
    #                 output_dict['MEAN'][time][key_size] = np.nanmean(distance)
    #                 output_dict['MEDIAN'][time][key_size] = np.nanmedian(distance)
    #                 output_dict['STD'][time][key_size] = np.nanstd(distance)
    # print('final mean {} median {} std {}'.format(output_dict['MEAN'][720][key_size],
    #                                               output_dict['MEDIAN'][720][key_size],
    #                                               output_dict['STD'][720][key_size]))

