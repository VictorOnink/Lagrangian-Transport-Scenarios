import settings as settings
import utils
from advection_scenarios import advection_files
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar
from copy import deepcopy


if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
    def parcels_to_sizespectrum(file_dict: dict):
        """
        For the FragmentationKaandorp scenario, getting a histogram of how many particles we have in the various size
        categories, where we distinguish between different beach states
        :param file_dict:
        :return:
        """
        # Setting the size bins
        bin_number = 20
        size_bins = np.logspace(start=-5, stop=-2, num=20)

        output_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/{}/'.format(
            settings.SCENARIO_NAME)
        utils.check_direc_exist(output_direc)

        # Loading the time axis
        time_list = np.array([], dtype=np.int)
        for restart in range(settings.SIM_LENGTH):
            parcels_file = file_dict['parcels'][settings.STARTYEAR][1][0][restart]
            parcels_dataset = Dataset(parcels_file)
            time = parcels_dataset.variables['time'][:, :-1]
            for time_case in np.unique(time):
                if np.nansum(time_case == time) > 10 and time_case != np.nan:
                    time_list = np.append(time_list, time_case)

        var_list = ['time', 'size', 'beach', 'z', 'distance2coast']
        for key in var_list:
            assert key in parcels_dataset.variables.keys(), '{} is not in the output file!!!'.format(key)

        # Creating the output dict
        beach_label = {'beach': 1, 'afloat': 0, 'seabed': 3, 'removed': 2}
        output_dict = {'size_bins': size_bins, 'beach': {}, 'adrift': {}, 'adrift_5m': {}, 'adrift_2m': {},
                       'adrift_10km': {}, 'adrift_20km': {}, 'seabed': {}, 'total': {}}
        time_step = 60
        for key in output_dict.keys():
            for index_time in range(0, len(time_list), time_step):
                output_dict[key][index_time] = np.zeros(bin_number - 1, dtype=float)

        pbar = ProgressBar()
        for ind_year, year in pbar(enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH))):
            for month in range(1, 13):
                for run in range(0, settings.RUN_RANGE):
                    # Loop through the restart files
                    for restart in range(0, settings.SIM_LENGTH - ind_year):
                        # Load the lon, lat, time, beach and weight data
                        parcels_file = file_dict['postprocess'][year][month][run][restart]
                        post_file = file_dict['postprocess'][year][month][run][restart]
                        parcels_dataset = Dataset(parcels_file)
                        dataset_post = utils.load_obj(post_file)
                        run_restart_dict = {}
                        for key in var_list:
                            run_restart_dict[key] = parcels_dataset.variables[key][:, :-1]
                        run_restart_dict['particle_number'] = dataset_post['particle_number']
                        # Calculating the spectrum every 30 days (note, one output step is 12 hours)
                        for index_time in range(0, len(time_list), time_step):
                            time_selection = time_list[index_time] == time
                            if time_selection.size > 0:
                                time_selection_dict = {}
                                for key in run_restart_dict.keys():
                                    time_selection_dict[key] = run_restart_dict[key][time_selection]
                                # All particles
                                size_counts, _ = np.histogram(time_selection_dict['size'], bins=size_bins,
                                                              weights=time_selection_dict['particle_number'])
                                output_dict['total'][index_time] += size_counts
                                # beached particles
                                selection = time_selection_dict['beach'] == beach_label['beach']
                                size_counts, _ = np.histogram(time_selection_dict['size'][selection], bins=size_bins,
                                                              weights=time_selection_dict['particle_number'][selection])
                                output_dict['beach'][index_time] += size_counts
                                # seabed particles
                                selection = time_selection_dict['beach'] == beach_label['seabed']
                                size_counts, _ = np.histogram(time_selection_dict['size'][selection], bins=size_bins,
                                                              weights=time_selection_dict['particle_number'][selection])
                                output_dict['seabed'][index_time] += size_counts
                                # floating particles
                                selection = time_selection_dict['beach'] == beach_label['adrift']
                                size_counts, _ = np.histogram(time_selection_dict['size'][selection], bins=size_bins,
                                                              weights=time_selection_dict['particle_number'][selection])
                                output_dict['adrift'][index_time] += size_counts
                                # floating particles within 5m of surface
                                selection = (time_selection_dict['beach'] == beach_label['adrift']) & (time_selection_dict['z'] < 5)
                                size_counts, _ = np.histogram(time_selection_dict['size'][selection], bins=size_bins,
                                                              weights=time_selection_dict['particle_number'][selection])
                                output_dict['adrift_5m'][index_time] += size_counts
                                # floating particles within 2m of surface
                                selection = (time_selection_dict['beach'] == beach_label['adrift']) & (time_selection_dict['z'] < 2)
                                size_counts, _ = np.histogram(time_selection_dict['size'][selection], bins=size_bins,
                                                              weights=time_selection_dict['particle_number'][selection])
                                output_dict['adrift_2m'][index_time] += size_counts
                                # Floating within 10 km of the model coastline
                                selection = (time_selection_dict['beach'] == beach_label['adrift']) & (time_selection_dict['distance2coast'] < 10)
                                size_counts, _ = np.histogram(time_selection_dict['size'][selection], bins=size_bins,
                                                              weights=time_selection_dict['particle_number'][selection])
                                output_dict['adrift_10km'][index_time] += size_counts
                                # Floating within 20 km of the model coastline
                                selection = (time_selection_dict['beach'] == beach_label['adrift']) & (time_selection_dict['distance2coast'] < 20)
                                size_counts, _ = np.histogram(time_selection_dict['size'][selection], bins=size_bins,
                                                              weights=time_selection_dict['particle_number'][selection])
                                output_dict['adrift_20km'][index_time] += size_counts

        # Adding the index of the final timestep for ease later on
        output_dict['final_index'] = index_time

        # Saving the output
        prefix = 'size_distribution'
        output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
        utils.save_obj(output_name, output_dict)
        utils.print_statement("The size distribution has been saved")

else:
    def parcels_to_sizespectrum(file_dict: dict):
        """
        For the FragmentationKaandorp scenario, getting a histogram of how many particles we have in the various size
        categories
        :param file_dict:
        :return:
        """
        # Setting the size bins
        size_bins = np.logspace(start=-5, stop=-2, num=20)

        output_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/{}/'.format(
            settings.SCENARIO_NAME)
        utils.check_direc_exist(output_direc)

        # Loading the time axis
        time_list = np.array([], dtype=np.int)
        for restart in range(settings.SIM_LENGTH):
            parcels_file = file_dict[0][restart]
            parcels_dataset = Dataset(parcels_file)
            time = parcels_dataset.variables['time'][:, :-1]
            for time_case in np.unique(time):
                if np.nansum(time_case == time) > 10 and time_case != np.nan:
                    time_list = np.append(time_list, time_case)

        # Creating the output dict
        output_dict = {'size_bins': size_bins}

        # loop through the runs
        for run in range(settings.RUN_RANGE):
            # Loop through the restart files
            for restart in range(settings.SIM_LENGTH):
                # Load the lon, lat, time, beach and weight data
                parcels_file = file_dict[run][restart]
                parcels_dataset = Dataset(parcels_file)
                time = parcels_dataset.variables['time'][:, :-1]
                size = parcels_dataset.variables['size'][:, :-1]
                if 'particle_number' in parcels_dataset.variables.keys():
                    particle_number = parcels_dataset.variables['particle_number'][:, :-1]

                # Calculating the spectrum every 30 days (note, one output step is 12 hours)
                for index_time in range(0, len(time_list), 60):
                    size_selection = size[time_list[index_time] == time]
                    utils.print_statement(size_selection.size, to_print=True)
                    if 'particle_number' in parcels_dataset.variables.keys():
                        particle_number_selection = particle_number[time_list[index_time] == time]
                        size_counts, _ = np.histogram(size_selection, bins=size_bins, weights=particle_number_selection)
                    else:
                        size_counts, _ = np.histogram(size_selection, bins=size_bins)
                    output_dict[index_time] = size_counts

        # Saving the output
        prefix = 'size_distribution'
        output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
        utils.save_obj(output_name, output_dict)
        utils.print_statement("The size distribution has been saved")
