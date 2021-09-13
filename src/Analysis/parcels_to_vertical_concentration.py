import settings as settings
import utils
from advection_scenarios import advection_files
from netCDF4 import Dataset
import numpy as np
import progressbar
import os
from copy import deepcopy
from datetime import datetime, timedelta

if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
    if settings.PARALLEL_STEP == 1:
        def parcels_to_vertical_concentration(file_dict: dict):
            # Set the bins used to calculate the vertical concentration
            depth_bins = determine_depth_bins()

            # Getting the directory saving the output files
            temp_direc, _ = get_directories(scenario_name=settings.SCENARIO_NAME)

            # Setting the time boundaries so that we can get the time points at the end of each month, which we in turn
            # can use to select all floating particles within a given time period
            time_list, days_in_month = determine_month_boundaries()

            # Create the output dictionary
            output_dict = create_output_dict(scenario_name=settings.SCENARIO_NAME, depth_bins=depth_bins)

            # Looping through all the simulation years and runs
            year, month, run, restart = settings.STARTYEAR, settings.STARTMONTH, settings.RUN, settings.RESTART
            print_statement = 'year {}-{}, run {} restart {}'.format(year, month, run, restart)

            utils.print_statement(print_statement, to_print=True)
            # Load the depth data, with monthly intervals
            parcels_file = file_dict['parcels'][year][month][run][restart]
            post_file = file_dict['postprocess'][year][month][run][restart]
            parcels_dataset = Dataset(parcels_file)
            post_dataset = utils.load_obj(post_file)
            run_restart_dict = {}
            for key in ['z', 'beach', 'size_class']:
                run_restart_dict[key] = parcels_dataset.variables[key][:, :-1].flatten()
            run_restart_dict['particle_number'] = post_dataset['particle_number'][:, :-1].flatten()
            time = parcels_dataset.variables['time'][:, :-1].flatten()

            for index_time in range(1, time_list.__len__()):
                month_index = (index_time - 1) % 12
                year_index = (index_time - 1) // 12
                selection = (time > time_list[index_time - 1]) & (time <= time_list[index_time])
                select_dict = {}
                for key in run_restart_dict.keys():
                    select_dict[key] = run_restart_dict[key][selection]
                # Picking out the non-beached particles
                selection = select_dict['beach'] == 0
                for key in ['z', 'particle_number', 'size_class']:
                    select_dict[key] = select_dict[key][selection]
                for size_class in range(settings.SIZE_CLASS_NUMBER):
                    size_dict = {}
                    for variable in ['z', 'particle_number']:
                        size_dict[variable] = select_dict[variable][select_dict['size_class'] == size_class]
                    histogram_counts, _ = np.histogram(size_dict['z'], bins=depth_bins, weights=size_dict['particle_number'])
                    # Divide the counts by the number of days in the month
                    histogram_counts /= days_in_month[index_time]
                    key_year = utils.analysis_simulation_year_key(year_index)
                    output_dict[key_year][month_index][size_class]['concentration'] += histogram_counts
                    output_dict[key_year][month_index][size_class]['counts'] += np.nansum(histogram_counts)

            # Saving the output
            output_name = get_file_names(scenario_name=settings.SCENARIO_NAME, file_dict=file_dict, directory=temp_direc,
                                         final=False)
            utils.save_obj(output_name, output_dict)
            utils.print_statement("The vertical concentration has been saved")

    if settings.PARALLEL_STEP == 2:
        def parcels_to_vertical_concentration(file_dict: dict):
            # Set the bins used to calculate the vertical concentration
            depth_bins = determine_depth_bins()

            # Getting the directory saving the output files
            temp_direc, output_direc = get_directories(scenario_name=settings.SCENARIO_NAME)

            # Setting the time boundaries so that we can get the time points at the end of each month, which we in turn
            # can use to select all floating particles within a given time period
            time_list, days_in_month = determine_month_boundaries()

            # Create the output dictionary
            output_dict = create_output_dict(scenario_name=settings.SCENARIO_NAME, depth_bins=depth_bins)

            # Looping through all the simulation years and runs
            pbar = progressbar.ProgressBar()
            for ind_year, year in pbar(enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH))):
                for month in range(1, 13):
                    for run in range(0, settings.RUN_RANGE):
                        # Loop through the restart files
                        for restart in range(0, settings.SIM_LENGTH - ind_year):
                            print_statement = 'year {}-{}, run {} restart {}'.format(year, month, run, restart)
                            utils.print_statement(print_statement, to_print=True)
                            # Load the depth data, with monthly intervals
                            parcels_file = file_dict['parcels'][year][month][run][restart]
                            post_file = file_dict['postprocess'][year][month][run][restart]
                            parcels_dataset = Dataset(parcels_file)
                            post_dataset = utils.load_obj(post_file)
                            run_restart_dict = {}
                            for key in ['z', 'beach', 'size_class']:
                                run_restart_dict[key] = parcels_dataset.variables[key][:, :-1].flatten()
                            run_restart_dict['particle_number'] = post_dataset['particle_number'][:, :-1].flatten()
                            time = parcels_dataset.variables['time'][:, :-1].flatten()

                            for index_time in range(1, time_list.__len__()):
                                month_index = (index_time - 1) % 12
                                year_index = (index_time - 1) // 12
                                selection = (time > time_list[index_time - 1]) & (time <= time_list[index_time])
                                select_dict = {}
                                for key in run_restart_dict.keys():
                                    select_dict[key] = run_restart_dict[key][selection]
                                # Picking out the non-beached particles
                                selection = select_dict['beach'] == 0
                                for key in ['z', 'particle_number', 'size_class']:
                                    select_dict[key] = select_dict[key][selection]
                                for size_class in range(settings.SIZE_CLASS_NUMBER):
                                    size_dict = {}
                                    for variable in ['z', 'particle_number']:
                                        size_dict[variable] = select_dict[variable][
                                            select_dict['size_class'] == size_class]
                                    histogram_counts, _ = np.histogram(size_dict['z'], bins=depth_bins,
                                                                       weights=size_dict['particle_number'])
                                    # Divide the counts by the number of days in the month
                                    histogram_counts /= days_in_month[index_time]
                                    key_year = utils.analysis_simulation_year_key(year_index)
                                    output_dict[key_year][month_index][size_class]['concentration'] += histogram_counts
                                    output_dict[key_year][month_index][size_class]['counts'] += np.nansum(
                                        histogram_counts)

            # Saving the output
            output_name = get_file_names(scenario_name=settings.SCENARIO_NAME, file_dict=file_dict,
                                         directory=output_direc, final=True)
            utils.save_obj(output_name, output_dict)
            utils.print_statement("The vertical concentration has been saved")

else:
    def parcels_to_vertical_concentration(file_dict: dict):
        # Set the bins used to calculate the vertical concentration
        depth_bins = determine_depth_bins()

        # Getting the directory saving the output files
        output_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format(settings.SCENARIO_NAME)
        utils.check_direc_exist(output_direc)

        # Setting the time boundaries so that we can get the time points at the end of each month, which we in turn can
        # use to select all floating particles within a given time period
        time_list, days_in_month = determine_month_boundaries()

        # Create the output dictionary
        output_dict = create_output_dict(scenario_name=settings.SCENARIO_NAME, depth_bins=depth_bins)

        # Looping through all the simulation years and runs
        pbar = progressbar.ProgressBar()
        for run in pbar(range(settings.RUN_RANGE)):
            # Loop through the restart files
            for restart in range(settings.SIM_LENGTH):
                key_year = utils.analysis_simulation_year_key(restart)
                # Load the depth data, with monthly intervals
                parcels_file = file_dict[run][restart]
                dataset = Dataset(parcels_file)
                depth = dataset.variables['z'][:]
                beach = dataset.variables['beach'][:]
                time_interval = depth.shape[1] // 365 * 31
                depth = depth[:, ::time_interval]
                beach = beach[:, ::time_interval]
                # Start looping through the different time intervals
                for month in range(depth.shape[1]):
                    # Checking to make sure we don't have any nan values, or particles with beach==2 (indicating a particle
                    # that was previously removed)
                    non_nan_values = (~np.isnan(depth[:, month])) & (beach[:, month] != 2)
                    depth_month = depth[:, month][non_nan_values]
                    # Now, calculating the vertical histogram
                    histogram_counts, _ = np.histogram(depth_month, bins=depth_bins)
                    output_dict[key_year][month] += histogram_counts

        # Saving the output
        output_name = get_file_names(scenario_name=settings.SCENARIO_NAME, file_dict=file_dict, directory=output_direc,
                                     final=True)
        utils.save_obj(output_name, output_dict)
        utils.print_statement("The vertical concentration has been saved")


########################################################################################################################
"""
These following functions are used across all scenarios
"""


def determine_depth_bins():
    advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                        advection_scenario=settings.ADVECTION_DATA,
                                                        repeat_dt=None)
    adv_file_dict = advection_scenario.file_names
    depth_min, depth_max = np.nanmin(adv_file_dict['DEPTH']), np.nanmax(adv_file_dict['DEPTH'])
    step = 0.3
    depth_bins = np.arange(depth_min, depth_max + step, step)
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


def create_output_dict(scenario_name, depth_bins):
    # Creating the output dict containing the depth bins, and creating the base_dict, which is lowest dictionary level
    # that will contain the particle counts and the concentration
    output_dict = {'depth': depth_bins}
    base_dict = {'counts': 0.0, 'concentration': np.zeros(depth_bins.__len__() - 1, dtype=np.float32)}
    for simulation_year in range(settings.SIM_LENGTH + 1):
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
    temp_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/temporary/'.format(scenario_name)
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
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix, split=split)
    return output_name

