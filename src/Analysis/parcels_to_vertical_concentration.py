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
    def parcels_to_vertical_concentration(file_dict: dict):
        # Get the depth range for the region in question, and using that to set the bins
        advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                            advection_scenario=settings.ADVECTION_DATA,
                                                            repeat_dt=None)
        adv_file_dict = advection_scenario.file_names
        depth_min, depth_max = np.nanmin(adv_file_dict['DEPTH']), np.nanmax(adv_file_dict['DEPTH'])
        step = 0.1
        depth_bins = np.arange(depth_min, depth_max + step, step)

        # Getting the directory saving the output files
        output_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format(
            settings.SCENARIO_NAME)
        utils.check_direc_exist(output_direc)

        # Determining the time selections so that we always compute the vertical concentration on the first day of the
        # month
        reference_time = datetime(2010, 1, 1, 12, 0)
        time_list = [-1e6]
        for year in range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH):
            for month in range(1, 13):
                last_of_month = datetime(year, month, 1, 0, 0) - timedelta(seconds=1)
                time_list.append((last_of_month - reference_time).total_seconds())

        # Create the output dictionary, and a dictionary to keep track of particle counts for the normalization
        output_dict = {'depth': depth_bins}
        for simulation_year in range(settings.SIM_LENGTH):
            key_year = utils.analysis_simulation_year_key(simulation_year)
            output_dict[key_year] = {}
            for month in range(0, 12):
                output_dict[key_year][month] = {}
                for size_class in range(settings.SIZE_CLASS_NUMBER):
                    output_dict[key_year][month][size_class] = {'counts': 0.0,
                                                                'concentration': np.zeros(depth_bins.__len__() - 1,
                                                                                          dtype=np.float32)}

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
                            month_index = index_time % 12
                            year_index = index_time // 12
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
                                key_year = utils.analysis_simulation_year_key(year_index)
                                output_dict[key_year][month_index][size_class]['concentration'] += histogram_counts
                                output_dict[key_year][month_index][size_class]['counts'] += np.nansum(histogram_counts)

        # Saving the output
        prefix = 'vertical_concentration'
        output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict['postprocess'][settings.STARTYEAR][1][0][restart],
                                                                   prefix=prefix)
        utils.save_obj(output_name, output_dict)
        utils.print_statement("The vertical concentration has been saved")

else:
    def parcels_to_vertical_concentration(file_dict: dict):
        # Get the depth range for the region in question, and using that to set the bins
        advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                            advection_scenario=settings.ADVECTION_DATA,
                                                            repeat_dt=None)
        adv_file_dict = advection_scenario.file_names
        depth_min, depth_max = np.nanmin(adv_file_dict['DEPTH']), np.nanmax(adv_file_dict['DEPTH'])
        step = 0.1
        depth_bins = np.arange(depth_min, depth_max + step, step)

        # Getting the directory saving the output files
        output_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format(settings.SCENARIO_NAME)
        utils.check_direc_exist(output_direc)

        # Create the output dictionary, and a dictionary to keep track of particle counts for the normalization
        output_dict = {'depth': depth_bins}
        counts_dict = {}
        for simulation_year in range(settings.SIM_LENGTH):
            key_year = utils.analysis_simulation_year_key(simulation_year)
            output_dict[key_year] = {}
            counts_dict[key_year] = {}
            for month in range(12):
                output_dict[key_year][month] = np.zeros(depth_bins.__len__() - 1, dtype=np.float32)
                counts_dict[key_year][month] = 0.0

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
                    # Updating the number of particles used within a certain profile so that we can normalize later
                    counts_dict[key_year][month] += depth_month.size

        # Now, we loop through all the profiles and normalize them by the number of particles within the profile
        for key_year in counts_dict.keys():
            for key_month in counts_dict[key_year].keys():
                output_dict[key_year][key_month] /= counts_dict[key_year][key_month]

        # Saving the output
        prefix = 'vertical_concentration'
        output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
        utils.save_obj(output_name, output_dict)
        utils.print_statement("The vertical concentration has been saved")
