import settings as settings
import utils
from advection_scenarios import advection_files
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar
import os
from copy import deepcopy
from datetime import datetime, timedelta


class parcels_to_vertical_concentration:
    def __init__(self, file_dict: dict):
        self.parallel_step = settings.PARALLEL_STEP
        self.file_dict = file_dict
        self.depth_bins = determine_depth_bins()
        self.temp_direc, self.output_direc = get_directories(scenario_name=settings.SCENARIO_NAME)
        # Some variables specific to FragmentationKaandorpPartial
        self.weight_list = ['particle_mass_sink', 'particle_number_sink']
        self.concentration_list = ['concentration_mass_sink',
                                   'concentration_number_sink']
        self.counts_list = ['counts_mass_sink', 'counts_number_sink']
        # Creating the output_dict
        self.output_dict = create_output_dict(scenario_name=settings.SCENARIO_NAME, depth_bins=self.depth_bins,
                                              concentration_list=self.concentration_list, counts_list=self.counts_list)

    def run(self):
        if self.parallel_step == 1:
            # Setting the time boundaries so that we can get the time points at the end of each month, which we in turn
            # can use to select all floating particles within a given time period
            time_list, days_in_month = determine_month_boundaries()

            # Loading the data
            year, month, run, restart = settings.STARTYEAR, settings.STARTMONTH, settings.RUN, settings.RESTART
            print_statement = 'year {}-{}, run {} restart {}'.format(year, month, run, restart)
            utils.print_statement(print_statement, to_print=True)
            output_name = get_file_names(scenario_name=settings.SCENARIO_NAME, file_dict=self.file_dict,
                                         directory=self.temp_direc, final=False)
            if not utils.check_file_exist(output_name, without_pkl=True):
                parcels_dataset, post_dataset = load_parcels_post_output(scenario_name=settings.SCENARIO_NAME,
                                                                         file_dict=self.file_dict)
                # Load the relevant fields into run_restart_dict and flattening the arrays
                run_restart_dict = set_run_restart_dict(settings.SCENARIO_NAME, parcels_dataset, post_dataset, self.weight_list)

                # Looping through the various months
                for index_time in range(1, time_list.__len__()):
                    month_index = (index_time - 1) % 12
                    key_year = utils.analysis_simulation_year_key((index_time - 1) // 12)
                    # selecting all points within the month in question
                    selection = (run_restart_dict['time'] > time_list[index_time - 1]) & (run_restart_dict['time'] <= time_list[index_time])
                    select_dict = {}
                    for key in run_restart_dict.keys():
                        select_dict[key] = run_restart_dict[key][selection]
                    # Picking out the non-beached particles
                    selection = select_dict['beach'] == 0
                    for key in select_dict.keys():
                        if key != 'beach':
                            select_dict[key] = select_dict[key][selection]
                    # If size_class is a variable, then we break down the calculation to separate size classes, otherwise
                    # we go straight to calculating the distribution
                    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                        for size_class in range(settings.SIZE_CLASS_NUMBER):
                            size_dict = {}
                            for variable in utils.flatten_list_of_lists([['z'], self.weight_list]):
                                size_dict[variable] = select_dict[variable][select_dict['size_class'] == size_class]
                            for index_weight, weight in enumerate(self.weight_list):
                                histogram_counts, _ = np.histogram(size_dict['z'], bins=self.depth_bins, weights=size_dict[weight])
                                # Divide the counts by the number of days in the month
                                histogram_counts /= days_in_month[index_time]
                                self.output_dict[key_year][month_index][size_class][self.concentration_list[index_weight]] += histogram_counts
                                self.output_dict[key_year][month_index][size_class][self.counts_list[index_weight]] += np.nansum(histogram_counts)
                    else:
                        histogram_counts, _ = np.histogram(select_dict['z'], bins=self.depth_bins, weights=select_dict['weights'])
                        # Divide the counts by the number of days in the month
                        histogram_counts /= days_in_month[index_time]
                        self.output_dict[key_year][month_index]['concentration'] += histogram_counts
                        self.output_dict[key_year][month_index]['counts'] += np.nansum(histogram_counts)


                utils.save_obj(output_name, self.output_dict)
                str_format = settings.STARTYEAR, settings.STARTMONTH, settings.RUN, settings.RESTART
                print_statement = 'The vertical concentration for year {}-{}, run {} restart {} has been save'.format(*str_format)
                utils.print_statement(print_statement, to_print=True)

        elif self.parallel_step == 2:
            pbar = ProgressBar()
            for ind_year, year in pbar(enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH))):
                for month in range(1, 13):
                    for run in range(0, settings.RUN_RANGE):
                        for restart in range(0, settings.SIM_LENGTH - ind_year):
                            if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                                file_name = get_file_names(scenario_name=settings.SCENARIO_NAME, file_dict=self.file_dict,
                                                           directory=self.temp_direc, final=False, year=year, month=month,
                                                           run=run, restart=restart)
                                dataset_post = utils.load_obj(filename=file_name)

                                for sim_year in range(settings.SIM_LENGTH):
                                    key_year = utils.analysis_simulation_year_key(sim_year)
                                    for month_index in self.output_dict[key_year].keys():
                                        for size_class in self.output_dict[key_year][month_index].keys():
                                            for index_conc, concentration in enumerate(self.concentration_list):
                                                self.output_dict[key_year][month_index][size_class][concentration] += dataset_post[key_year][month_index][size_class][concentration]
                                                self.output_dict[key_year][month_index][size_class][self.counts_list[index_conc]] += dataset_post[key_year][month_index][size_class][self.counts_list[index_conc]]
                                utils.remove_file(file_name)
                            else:
                                if month == 1 and year == settings.STARTYEAR:
                                    file_name = get_file_names(scenario_name=settings.SCENARIO_NAME, file_dict=self.file_dict,
                                                               directory=self.temp_direc, final=False, year=year, run=run,
                                                               restart=restart)
                                    dataset_post = utils.load_obj(filename=file_name)
                                    for sim_year in range(settings.SIM_LENGTH):
                                        key_year = utils.analysis_simulation_year_key(sim_year)
                                        for month_index in self.output_dict[key_year].keys():
                                            self.output_dict[key_year][month_index]['concentration'] += dataset_post[key_year][month_index]['concentration']
                                            self.output_dict[key_year][month_index]['counts'] += dataset_post[key_year][month_index]['counts']
                                    utils.remove_file(file_name)
            # Saving the output
            output_name = get_file_names(scenario_name=settings.SCENARIO_NAME, file_dict=self.file_dict,
                                         directory=self.output_direc, final=True)
            utils.save_obj(output_name, self.output_dict)
            utils.print_statement("The vertical concentration has been saved")
        else:
            ValueError('settings.PARALLEL_STEP can not have a value of {}'.format(self.parallel_step))


########################################################################################################################
"""
These following functions are used across all scenarios
"""


def determine_depth_bins():
    advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                        advection_scenario=settings.ADVECTION_DATA,
                                                        repeat_dt=None)
    adv_file_dict = advection_scenario.file_names
    depth_min, depth_max = np.nanmin(adv_file_dict['DEPTH']) - 0.1, np.nanmax(adv_file_dict['DEPTH'])
    if depth_min > 0:
        log_min, log_max = np.log10(depth_min), np.log10(depth_max)
    else:
        log_min, log_max = -1, np.log10(depth_max)
    num = 100
    depth_bins = np.logspace(log_min, log_max, num=num)
    return depth_bins


def determine_month_boundaries():
    reference_time = datetime(2010, 1, 1, 12, 0)
    time_list = [-1e6]
    days_in_month = [1]
    for year in range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH + 1):
        for month in range(1, 13):
            last_of_month = datetime(year, month, 1, 0, 0) - timedelta(seconds=1)
            days_in_month.append(last_of_month.day)
            time_list.append((last_of_month - reference_time).total_seconds())
    return time_list, days_in_month


def create_output_dict(scenario_name, depth_bins, concentration_list, counts_list):
    # Creating the output dict containing the depth bins, and creating the base_dict, which is lowest dictionary level
    # that will contain the particle counts and the concentration
    depth_mid = 0.5 * depth_bins[1:] + 0.5 * depth_bins[:-1]
    output_dict = {'depth': depth_mid}
    if scenario_name in ['FragmentationKaandorpPartial']:
        base_dict = {}
        for count in counts_list:
            base_dict[count] = 0.0
        for concentration in concentration_list:
            base_dict[concentration] = np.zeros(depth_bins.__len__() - 1, dtype=np.float32)
        simulation_range = settings.SIM_LENGTH + 1
    else:
        base_dict = {'counts': 0.0, 'concentration': np.zeros(depth_bins.__len__() - 1, dtype=np.float32)}
        simulation_range = settings.SIM_LENGTH + 1
    for simulation_year in range(simulation_range):
        key_year = utils.analysis_simulation_year_key(simulation_year)
        output_dict[key_year] = {}
        for month in range(0, 12):
            if scenario_name in ['FragmentationKaandorpPartial']:
                output_dict[key_year][month] = {}
                for size_class in range(settings.SIZE_CLASS_NUMBER):
                    output_dict[key_year][month][size_class] = deepcopy(base_dict)
            else:
                output_dict[key_year][month] = deepcopy(base_dict)
    return output_dict


def get_directories(scenario_name):
    # temp_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/temporary/'.format(scenario_name)
    temp_direc = settings.SCRATCH_DIR
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format(scenario_name)
    utils.check_direc_exist(temp_direc)
    utils.check_direc_exist(output_direc)
    return temp_direc, output_direc


def get_file_names(scenario_name, file_dict, directory, final, year=settings.STARTYEAR, month=settings.STARTMONTH,
                   run=settings.RUN, restart=settings.RESTART):
    split = {True: None, False: '.nc'}[final]
    prefix = 'vertical_concentration'
    if scenario_name in ['FragmentationKaandorpPartial']:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict['postprocess'][year][month][run][restart],
        prefix=prefix, split=split)
    else:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict[run][restart], prefix=prefix, split=split)
    return output_name


def load_parcels_post_output(scenario_name, file_dict, year=settings.STARTYEAR, month=settings.STARTMONTH,
                             run=settings.RUN, restart=settings.RESTART):
    if scenario_name in ['FragmentationKaandorpPartial']:
        parcels_dataset = Dataset(file_dict['parcels'][year][month][run][restart])
        post_dataset = utils.load_obj(file_dict['postprocess'][year][month][run][restart])
    else:
        parcels_dataset = Dataset(file_dict[run][restart])
        post_dataset = None
    return parcels_dataset, post_dataset


def set_run_restart_dict(scenario_name, parcels_dataset, post_dataset, weight_list):
    run_restart_dict = {}
    for variable in ['z', 'beach', 'time']:
        run_restart_dict[variable] = parcels_dataset.variables[variable][:, :-1].flatten()
    if scenario_name in ['FragmentationKaandorpPartial']:
        run_restart_dict['size_class'] = parcels_dataset.variables['size_class'][:, :-1].flatten()
        for weight in weight_list:
            run_restart_dict[weight] = post_dataset[weight][:, :-1].flatten()
    else:
        run_restart_dict['weights'] = np.ones(run_restart_dict[variable].shape, dtype=float)
    return run_restart_dict
