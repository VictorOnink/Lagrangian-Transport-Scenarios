import settings as settings
import utils
from advection_scenarios import advection_files
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar
from copy import deepcopy


if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
    def parcels_to_sizespectrum(file_dict: dict, scenario):
        """
        For the FragmentationKaandorp scenario, getting a histogram of how many particles we have in the various size
        categories, where we distinguish between different beach states
        :param file_dict:
        :return:
        """
        # Setting the size bins
        bin_number = settings.SIZE_CLASS_NUMBER

        output_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/{}/'.format(
            settings.SCENARIO_NAME)
        utils.check_direc_exist(output_direc)

        # Getting the minimum depth of the advection scenario
        min_depth = np.nanmin(scenario.file_dict['DEPTH'])

        # Loading the time axis
        time_list = np.array([], dtype=np.int)
        for restart in range(settings.SIM_LENGTH):
            parcels_file = file_dict['parcels'][settings.STARTYEAR][1][0][restart]
            parcels_dataset = Dataset(parcels_file)
            time = parcels_dataset.variables['time'][:, :-1]
            for time_case in np.unique(time):
                if np.nansum(time_case == time) > 10 and time_case != np.nan:
                    time_list = np.append(time_list, time_case)

        var_list = ['time', 'size_class', 'beach', 'z', 'distance2coast']
        for key in var_list:
            assert key in parcels_dataset.variables.keys(), '{} is not in the output file!!!'.format(key)

        # Creating the output dict
        beach_label = {'beach': 1, 'adrift': 0, 'seabed': 3, 'removed': 2}
        output_dict = {'size_bins': range(bin_number), 'beach': {}, 'adrift': {}, 'adrift_5m': {}, 'adrift_2m': {},
                       'adrift_10km': {}, 'adrift_10km_surf': {}, 'seabed': {}, 'total': {}, 'adrift_open': {},
                       'adrift_open_surf': {}}
        time_step = 60
        for key in output_dict.keys():
            if key not in ['size_bins']:
                for index_time in range(0, len(time_list), time_step):
                    output_dict[key][index_time] = np.zeros(shape=bin_number, dtype=float)

        pbar = ProgressBar()
        for ind_year, year in pbar(enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH))):
            for month in range(1, 13):
                for run in range(0, settings.RUN_RANGE):
                    # Loop through the restart files
                    for restart in range(0, settings.SIM_LENGTH - ind_year):
                        print_statement = 'year {}-{}, run {} restart {}'.format(year, month, run, restart)
                        utils.print_statement(print_statement, to_print=True)
                        # Load the lon, lat, time, beach and weight data
                        parcels_file = file_dict['parcels'][year][month][run][restart]
                        post_file = file_dict['postprocess'][year][month][run][restart]
                        parcels_dataset = Dataset(parcels_file)
                        dataset_post = utils.load_obj(post_file)
                        run_restart_dict = {}
                        for key in var_list:
                            run_restart_dict[key] = parcels_dataset.variables[key][:, :-1]
                        run_restart_dict['particle_number'] = dataset_post['particle_number'][:, :-1]
                        # Calculating the spectrum every 30 days (note, one output step is 12 hours)
                        for index_time in range(0, len(time_list), time_step):
                            time_selection = time_list[index_time] == run_restart_dict['time']
                            if time_selection.size > 0:
                                time_sel = {}
                                for key in run_restart_dict.keys():
                                    time_sel[key] = run_restart_dict[key][time_selection]
                                # All particles
                                output_dict['total'][index_time] += number_per_size_class(time_sel['size_class'], time_sel['particle_number'], bin_number)
                                # beached particles
                                selection = time_sel['beach'] == beach_label['beach']
                                output_dict['beach'][index_time] += number_per_size_class(time_sel['size_class'], time_sel['particle_number'], bin_number, selection=selection)
                                # seabed particles
                                selection = time_sel['beach'] == beach_label['seabed']
                                output_dict['seabed'][index_time] += number_per_size_class(time_sel['size_class'], time_sel['particle_number'], bin_number, selection=selection)
                                # floating particles
                                selection = time_sel['beach'] == beach_label['adrift']
                                output_dict['adrift'][index_time] += number_per_size_class(time_sel['size_class'], time_sel['particle_number'], bin_number, selection=selection)
                                # floating particles within 5m of surface
                                selection = (time_sel['beach'] == beach_label['adrift']) & (time_sel['z'] < (min_depth + 5))
                                output_dict['adrift_5m'][index_time] += number_per_size_class(time_sel['size_class'], time_sel['particle_number'], bin_number, selection=selection)
                                # floating particles within 2m of surface
                                selection = (time_sel['beach'] == beach_label['adrift']) & (time_sel['z'] < (min_depth + 2))
                                output_dict['adrift_2m'][index_time] += number_per_size_class(time_sel['size_class'], time_sel['particle_number'], bin_number, selection=selection)
                                # Floating within 10 km of the model coastline
                                selection = (time_sel['beach'] == beach_label['adrift']) & (time_sel['distance2coast'] < 10)
                                output_dict['adrift_10km'][index_time] += number_per_size_class(time_sel['size_class'], time_sel['particle_number'], bin_number, selection=selection)
                                # Floating within 20 km of the model coastline
                                selection = (time_sel['beach'] == beach_label['adrift']) & (time_sel['distance2coast'] < 10 & (time_sel['z'] < (min_depth + 0.26)))
                                output_dict['adrift_10km_surf'][index_time] += number_per_size_class(time_sel['size_class'], time_sel['particle_number'], bin_number, selection=selection)
                                # Floating beyond 10 km of the model coastline
                                selection = (time_sel['beach'] == beach_label['adrift']) & (time_sel['distance2coast'] > 10)
                                output_dict['adrift_open'][index_time] += number_per_size_class(time_sel['size_class'], time_sel['particle_number'], bin_number, selection=selection)
                                # Floating beyond 10 km of the model coastline and within 26cm of the surface
                                selection = (time_sel['beach'] == beach_label['adrift']) & (time_sel['distance2coast'] > 10) & (time_sel['z'] < (min_depth + 0.26))
                                output_dict['adrift_open_surf'][index_time] += number_per_size_class(time_sel['size_class'], time_sel['particle_number'], bin_number, selection=selection)

        # Adding the index of the final timestep for ease later on
        output_dict['final_index'] = index_time

        # Saving the output
        prefix = 'size_distribution'
        output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict['postprocess'][settings.STARTYEAR][1][0][restart], prefix=prefix)
        utils.save_obj(output_name, output_dict)
        utils.print_statement("The size distribution has been saved")

else:
    def parcels_to_sizespectrum(file_dict: dict, scenario):
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


def number_per_size_class(size_class_array, particle_number_array, bin_number, selection=None):
    assert size_class_array.size == particle_number_array.size, 'The size and number arrays must have the same size'
    output_array = np.zeros(shape=bin_number, dtype=float)
    if selection is not None:
        assert selection.size == size_class_array.size, 'The selection and data arrays must have the same size'
        for size_class in range(bin_number):
            in_size_class = size_class_array[selection] == size_class
            if np.nansum(in_size_class) > 0:
                output_array[size_class] += np.nansum(particle_number_array[selection][in_size_class])
            else:
                output_array[size_class] += 0
    else:
        for size_class in range(bin_number):
            output_array[size_class] += np.nansum(particle_number_array[size_class_array == size_class])
    return output_array