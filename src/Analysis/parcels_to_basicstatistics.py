import settings as settings
import utils
from advection_scenarios import advection_files
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import progressbar
import os
from copy import deepcopy


class parcels_to_basicstatistics:
    """
    This is a general function that we can use to calculate max, min and mean values of particle values. It won't be
    equally meaningful for each variable, but these properties should be relatively cheap to compute
    :param file_dict:
    :return:
    """

    def __init__(self, file_dict: dict):
        self.parallel_step = settings.PARALLEL_STEP
        self.file_dict = file_dict
        self.temp_direc, self.output_direc = get_directories(scenario_name=settings.SCENARIO_NAME)
        self.variable_list = ['z', 'horizontal_distance', 'vertical_distance']
        self.stats_list = ['mean', 'max', 'min', 'std', 'count']

    def run(self):
        if self.parallel_step == 1:
            parcels_dataset, _ = load_parcels_post_output(file_dict=self.file_dict)
            output_dict = create_output_dict(dataset=parcels_dataset, stats_list=self.stats_list)
            # Loading the data
            variable_list = set_variable_list(parcels_dataset)
            data_dict = {}
            for variable in self.variable_list:
                if variable in parcels_dataset.variables.keys():
                    data_dict[variable] = parcels_dataset.variables[variable][:, :-1]


        elif self.parallel_step == 2:

        else:
            ValueError('settings.PARALLEL_STEP can not have a value of {}'.format(self.parallel_step))

def create_output_dict(dataset, stats_list):
    base_array = np.zeros(dataset.dimensions['traj'].size,)

def parcels_to_basicstatistics(file_dict: dict):
    """
    This is a general function that we can use to calculate max, min and mean values of particle values. It won't be
    equally meaningful for each variable, but these properties should be relatively cheap to compute
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
    statistics_dict = dict.fromkeys(['mean', 'max', 'min', 'std', 'count'])
    for key in statistics_dict.keys():
        statistics_dict[key] = deepcopy(base_array)
    beach_label_dict = {'beach': 1, 'afloat': 0, 'seabed': 3, 'total': 2}
    beach_state_dict = dict.fromkeys(beach_label_dict.keys())
    for key in beach_state_dict.keys():
        beach_state_dict[key] = deepcopy(statistics_dict)
    output_dict = dict.fromkeys(variable_list)
    for key in output_dict:
        output_dict[key] = deepcopy(beach_state_dict)

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
                if beach_state == 'total':
                    state_variable_array[beach_array == beach_label_dict[beach_state]] = np.nan
                else:
                    state_variable_array[beach_array != beach_label_dict[beach_state]] = np.nan

                variable_mean = np.nanmean(state_variable_array, axis=1, keepdims=True)
                output_dict[variable][beach_state]['mean'] = np.concatenate((output_dict[variable][beach_state]['mean'],
                                                                             variable_mean), axis=0)

                variable_max = np.nanmax(state_variable_array, axis=1, keepdims=True)
                output_dict[variable][beach_state]['max'] = np.concatenate((output_dict[variable][beach_state]['max'],
                                                                            variable_max), axis=0)

                variable_min = np.nanmin(state_variable_array, axis=1, keepdims=True)
                output_dict[variable][beach_state]['min'] = np.concatenate((output_dict[variable][beach_state]['min'],
                                                                            variable_min), axis=0)

                variable_std = np.nanstd(state_variable_array, axis=1, keepdims=True)
                output_dict[variable][beach_state]['std'] = np.concatenate((output_dict[variable][beach_state]['std'],
                                                                            variable_std), axis=0)

                variable_count = np.nansum(beach_array == beach_label_dict[beach_state], axis=1, keepdims=True)
                output_dict[variable][beach_state]['count'] = np.concatenate((output_dict[variable][beach_state]['count'],
                                                                              variable_count), axis=0)

    # remove the first element of each array, as this was a dummy that was just there to initialize the array
    for variable in variable_list:
        for beach_state in beach_label_dict.keys():
            for statistic in statistics_dict.keys():
                output_dict[variable][beach_state][statistic] = output_dict[variable][beach_state][statistic][1:, :]
                print("{}, {}, {}".format(variable, beach_state, np.nanmax(output_dict[variable][beach_state][statistic])))

    # Saving the output file
    prefix = 'basic_statistics'
    output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
    utils.save_obj(output_name, output_dict)
    utils.print_statement("The basic statistics output has been saved")


########################################################################################################################
"""
These following functions are used across all scenarios
"""


def get_directories(scenario_name):
    temp_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/{}/temporary/'.format(scenario_name)
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/{}/'.format(scenario_name)
    utils.check_direc_exist(temp_direc)
    utils.check_direc_exist(output_direc)
    return temp_direc, output_direc


def load_parcels_post_output(file_dict, year=settings.STARTYEAR, month=settings.STARTMONTH,
                             run=settings.RUN, restart=settings.RESTART):
    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
        parcels_dataset = Dataset(file_dict['parcels'][year][month][run][restart])
        post_dataset = utils.load_obj(file_dict['postprocess'][year][month][run][restart])
    else:
        parcels_dataset = Dataset(file_dict[run][restart])
        post_dataset = None
    return parcels_dataset, post_dataset


def get_file_names(file_dict, directory, final, year=settings.STARTYEAR, month=settings.STARTMONTH,
                   run=settings.RUN, restart=settings.RESTART):
    split = {True: None, False: '.nc'}[final]
    prefix = 'basic_statistics'
    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict['postprocess'][year][month][run][restart],
        prefix=prefix, split=split)
    else:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix, split=split)
    return output_name


def set_variable_list(parcels_dataset):
    variable_list = list(parcels_dataset.variables.keys())
    for remove_variable in ['time', 'trajectory', 'lat', 'lon', 'beach', 'age', 'reynolds', 'size',
                            'to_split', 'size_class', 'size', 'parent']:
        if remove_variable in variable_list:
            variable_list.remove(remove_variable)
    return variable_list