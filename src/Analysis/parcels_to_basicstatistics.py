import settings as settings
import utils
from advection_scenarios import advection_files
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar
import os
from copy import deepcopy


class parcels_to_basicstatistics:
    """
    This is a general function that we can use to calculate max, min and mean values of particle values. It won't be
    equally meaningful for each variable, but these properties should be relatively cheap to compute.

    For ease of calculation (and because it isn't clear yet if these calculations will be necessary for fragmentation
    scenario analysis), we weigh all parcels particle instances equally, without adding corrections for variations
    in the particle number or mass over time.
    :param file_dict:
    :return:
    """

    def __init__(self, file_dict: dict):
        self.parallel_step = settings.PARALLEL_STEP
        self.file_dict = file_dict
        self.temp_direc, self.output_direc = get_directories(scenario_name=settings.SCENARIO_NAME)
        self.stats_list = ['mean', 'max', 'min', 'std']
        self.variable_list = compute_variable_list(file_dict)
        self.beach_label_dict = set_beach_label_dict()
        self.output_dict = create_output_dict(stats_list=self.stats_list, beach_state_list=self.beach_label_dict.keys(),
                                              variable_list=self.variable_list)

    def run(self):
        if self.parallel_step == 1:
            utils.print_statement('Nothing happens for parcels_to_basicstatistics when settings.PARALLEL_STEP == 1',
                                  to_print=True)
        elif self.parallel_step == 2:
            pbar = ProgressBar()
            for ind_year, year in pbar(enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH))):
                for month in range(1, 13):
                    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial'] or (month == 1) & (year == settings.STARTYEAR):
                        for run in range(0, settings.RUN_RANGE):
                            utils.print_statement('{}-{}, run = {}'.format(year, month, run), to_print=True)
                            run_dict = {}
                            # Loading the data over all restart files
                            for restart in range(0, settings.SIM_LENGTH - ind_year):
                                parcels_dataset, _ = load_parcels_post_output(self.file_dict, year, month, run, restart)
                                for variable in utils.flatten_list_of_lists([['beach', 'size_class'], self.variable_list]):
                                    run_dict = concatenate_variable(run_dict, variable, restart, parcels_dataset)
                            # Looping through the beach states
                            # Calculating the statistical values
                            for beach_state in self.beach_label_dict.keys():
                                beach_data_dict = deepcopy(run_dict)
                                if beach_state == 'total':
                                    beach_selection = beach_data_dict['beach'] != self.beach_label_dict[beach_state]
                                else:
                                    beach_selection = beach_data_dict['beach'] == self.beach_label_dict[beach_state]
                                for variable in beach_data_dict.keys():
                                    beach_data_dict[variable][beach_selection == False] = np.nan

                                if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                                    for size_class in range(settings.SIZE_CLASS_NUMBER):
                                        size_data_dict = deepcopy(beach_data_dict)
                                        for variable in self.variable_list:
                                            size_data_dict[variable][size_data_dict['size_class'] != size_class] = np.nan
                                            for statistic in self.stats_list:
                                                result = calculate_statistic(statistic=statistic, data=size_data_dict[variable])
                                                if result is None:
                                                    result = np.ones((beach_data_dict[variable].shape[0], 1)) * 1e20
                                                self.output_dict[variable][beach_state][size_class][statistic] = np.concatenate((self.output_dict[variable][beach_state][size_class][statistic],
                                                                                                                                 result), axis=0)
                                else:
                                    for variable in self.variable_list:
                                        for statistic in self.stats_list:
                                            result = calculate_statistic(statistic=statistic, data=beach_data_dict[variable])
                                            if result is None:
                                                result = np.ones((beach_data_dict[variable].shape[0], 1)) * 1e20
                                            self.output_dict[variable][beach_state][statistic] = np.concatenate(
                                                (self.output_dict[variable][beach_state][statistic], result), axis=0)
            # remove the first element of each array, as this was a dummy that was just there to initialize the array
            for variable in self.variable_list:
                for beach_state in self.beach_label_dict.keys():
                    for statistic in self.stats_list:
                        if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                            for size_class in range(settings.SIZE_CLASS_NUMBER):
                                self.output_dict[variable][beach_state][size_class][statistic] = self.output_dict[variable][beach_state][size_class][statistic][1:, :]
                        else:
                            self.output_dict[variable][beach_state][statistic] = self.output_dict[variable][beach_state][statistic][1:, :]
            # Saving the output
            file_name = get_file_names(file_dict=self.file_dict, directory=self.output_direc, final=True)
            utils.save_obj(filename=file_name, item=self.output_dict)
            utils.print_statement('The statistics have been saved', to_print=True)
        else:
            ValueError('settings.PARALLEL_STEP can not have a value of {}'.format(self.parallel_step))


########################################################################################################################
"""
These following functions are used across all scenarios
"""


def get_directories(scenario_name):
    temp_direc = settings.DATA_OUTPUT_DIREC + 'statistics/{}/temporary/'.format(scenario_name)
    output_direc = settings.DATA_OUTPUT_DIREC + 'statistics/{}/'.format(scenario_name)
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
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict[run][restart], prefix=prefix, split=split)
    return output_name


def set_variable_list(parcels_dataset):
    variable_list = list(parcels_dataset.variables.keys())
    for remove_variable in ['time', 'trajectory', 'lat', 'lon', 'beach', 'age', 'reynolds', 'size',
                            'to_split', 'size_class', 'size', 'parent']:
        if remove_variable in variable_list:
            variable_list.remove(remove_variable)
    return variable_list


def create_output_dict(stats_list, beach_state_list, variable_list):
    output_dict = {}
    for variable in variable_list:
        output_dict[variable] = {}
        for beach_state in beach_state_list:
            output_dict[variable][beach_state] = {}
            if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                for size_class in range(settings.SIZE_CLASS_NUMBER):
                    output_dict[variable][beach_state][size_class] = {}
                    for stats in stats_list:
                        output_dict[variable][beach_state][size_class][stats] = np.zeros((1, 1), dtype=float)
            else:
                for stats in stats_list:
                    output_dict[variable][beach_state][stats] = np.zeros((1, 1), dtype=float)
    return output_dict


def concatenate_variable(run_dict, variable, restart, parcels_dataset):
    if variable in parcels_dataset.variables.keys():
        variable_array = parcels_dataset.variables[variable][:, :-1]
        if restart == 0:
            run_dict[variable] = variable_array
        else:
            # Check if the number of particles matches
            if variable_array.shape[0] == run_dict[variable].shape[0]:
                run_dict[variable] = np.concatenate((run_dict[variable], variable_array), axis=1)
            else:
                additional_row = np.ones(run_dict[variable].shape[1]) * np.nan
                while run_dict[variable].shape[0] < variable_array.shape[0]:
                    run_dict[variable] = np.concatenate((run_dict[variable], [additional_row]), axis=0)
                run_dict[variable] = np.concatenate((run_dict[variable], variable_array), axis=1)
    return run_dict


def set_beach_label_dict():
    return {'beach': 1, 'adrift': 0, 'seabed': 3, 'total': 2}


def compute_variable_list(file_dict):
    parcels_dataset, _ = load_parcels_post_output(file_dict)
    variable_list = set_variable_list(parcels_dataset)
    return variable_list


def calculate_statistic(statistic, data):
    if statistic == 'mean':
        return np.nanmean(data, axis=1, keepdims=True)
    elif statistic == 'max':
        np.nanmax(data, axis=1, keepdims=True)
    elif statistic == 'min':
        np.nanmin(data, axis=1, keepdims=True)
    elif statistic == 'std':
        np.nanstd(data, axis=1, keepdims=True)
    else:
        ValueError('What do you mean by the test {}?'.format(statistic))


