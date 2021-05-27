import settings as settings
import utils
from advection_scenarios import advection_files
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import progressbar
import os
from copy import deepcopy


def parcels_to_basicstatistics(file_dict: dict):
    """
    This is a general function that we can use to calculate max, min and mean values of particle values. It won't be
    equally meaningfun for each variable, but these properties should be relatively cheap to compute
    :param file_dict:
    :return:
    """
    # Getting the directory saving the output files
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'statistics/{}/'.format(settings.SCENARIO_NAME)
    utils.check_direc_exist(output_direc)

    # Setting the parameters we want, and creating the output dict. We remove some variables from the variable list
    # since these will essentially never be used for this sort of analysis
    dataset = Dataset(file_dict[0][0])
    variable_list = list(dataset.variables.keys())
    for remove_variable in ['time', 'trajectory', 'lat', 'lon', 'beach', 'age', 'reynolds', 'size']:
        if remove_variable in variable_list:
            variable_list.remove(remove_variable)
    base_array = np.zeros((1, 1), dtype=float)
    statistics_dict = dict.fromkeys(['mean', 'max', 'min', 'std', 'count'], deepcopy(base_array))
    beach_label_dict = {'beach': 1, 'afloat': 0, 'seabed': 3, 'total': 2}
    beach_state_dict = dict.fromkeys(beach_label_dict.keys(), deepcopy(statistics_dict))
    output_dict = dict.fromkeys(variable_list, deepcopy(beach_state_dict))

    # Loading in all the beach state data
    for run in progressbar.progressbar(range(settings.RUN_RANGE)):
        for restart in range(settings.SIM_LENGTH):
            # Load parcels file
            parcels_file = file_dict[run][restart]
            dataset = Dataset(parcels_file)
            if restart == 0:
                beach_array = dataset.variables['beach'][:, :-1]
            else:
                beach_array = np.concatenate((beach_array, dataset.variables['beach'][:, :-1]), axis=1)

    # loop through the runs
    for run in progressbar.progressbar(range(settings.RUN_RANGE)):
        for variable in variable_list:
            # Loading all the data for that specific variable and appending them together. This should be ok in terms of
            # RAM if we only keep one variable in memory at a time
            for restart in range(settings.SIM_LENGTH):
                # Load parcels file
                parcels_file = file_dict[run][restart]
                dataset = Dataset(parcels_file)
                if restart == 0:
                    variable_array = dataset.variables[variable][:, :-1]
                else:
                    variable_array = np.concatenate((variable_array, dataset.variables[variable][:, :-1]), axis=1)
            # Now, starting to compute the various parameters
            for beach_state in beach_label_dict.keys():
                state_variable_array = deepcopy(variable_array)
                state_variable_array[beach_array != beach_label_dict[beach_state]] = np.nan
                output_dict[variable][beach_state]['mean'] = np.concatenate((output_dict[variable][beach_state]['mean'],
                                                                             np.nanmean(state_variable_array, axis=1,
                                                                                        keepdims=True)),
                                                                            axis=0)
                output_dict[variable][beach_state]['max'] = np.concatenate((output_dict[variable][beach_state]['max'],
                                                                            np.nanmax(state_variable_array, axis=1,
                                                                                      keepdims=True)),
                                                                           axis=0)
                output_dict[variable][beach_state]['min'] = np.concatenate((output_dict[variable][beach_state]['min'],
                                                                            np.nanmin(state_variable_array, axis=1,
                                                                                      keepdims=True)),
                                                                           axis=0)
                output_dict[variable][beach_state]['std'] = np.concatenate((output_dict[variable][beach_state]['std'],
                                                                            np.nanstd(state_variable_array, axis=1,
                                                                                      keepdims=True)),
                                                                           axis=0)
                output_dict[variable][beach_state]['count'] = np.concatenate((output_dict[variable][beach_state]['count'],
                                                                              np.nansum(
                                                                                  beach_array == beach_label_dict[
                                                                                      beach_state],
                                                                                  axis=1,
                                                                                  keepdims=True)),
                                                                             axis=0)
    # remove the first element of each array, as this was a dummy that was just there to initialize the array
    for variable in variable_list:
        for beach_state in beach_label_dict.keys():
            for statistic in statistics_dict.keys():
                output_dict[variable][beach_state][statistic] = output_dict[variable][beach_state][statistic][1:, :]

    # Saving the output file
    prefix = 'basic_statistics'
    output_name = output_direc + utils._analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
    utils.save_obj(output_name, output_dict)
    os.system('echo "The concentration has been saved"')
