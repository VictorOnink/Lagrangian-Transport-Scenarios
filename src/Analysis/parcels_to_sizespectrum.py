import settings as settings
import utils
from advection_scenarios import advection_files
from netCDF4 import Dataset
import numpy as np
import progressbar
from copy import deepcopy


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


def parcels_to_sizespectrum_beachstate(file_dict: dict):
    """
    For the FragmentationKaandorp scenario, getting a histogram of how many particles we have in the various size
    categories, where we distinguish between different beach states
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

    var_list = ['time', 'size', 'beach', 'z', 'particle_number']
    for key in var_list:
        assert key in parcels_dataset.variables.keys(), '{} is not in the output file!!!'.format(key)

    # Creating the output dict
    beach_state_dict = {'beach': np.zeros(time_list.shape, dtype=float),
                        'afloat': np.zeros(time_list.shape, dtype=float),
                        'seabed': np.zeros(time_list.shape, dtype=float),
                        'removed': np.zeros(time_list.shape, dtype=float),
                        'total': np.zeros(time_list.shape, dtype=float),
                        'time': time_list}
    beach_label = {'beach': 1, 'afloat': 0, 'seabed': 3, 'removed': 2}

    output_dict = {'size_bins': size_bins, 'beach': {}, 'afloat': {}, 'afloat_5m': {}, 'seabed': {},
                   'total': {}}

    # loop through the runs
    for run in range(settings.RUN_RANGE):
        # Loop through the restart files
        for restart in range(settings.SIM_LENGTH):
            # Load the lon, lat, time, beach and weight data
            parcels_file = file_dict[run][restart]
            parcels_dataset = Dataset(parcels_file)
            run_restart_dict = {}
            for key in var_list:
                run_restart_dict[key] = parcels_dataset.variables[key][:, :-1]

            # Calculating the spectrum every 30 days (note, one output step is 12 hours)
            for index_time in range(0, len(time_list), 60):
                time_selection = time_list[index_time] == time
                time_selection_dict = {}
                for key in var_list:
                    time_selection_dict[key] = run_restart_dict[key][time_selection]
                # All particles
                size_counts, _ = np.histogram(time_selection_dict['size'], bins=size_bins,
                                              weights=time_selection_dict['particle_number'])
                output_dict['total'][index_time] = size_counts
                # beached particles
                beach_selection = time_selection_dict['beach'] == beach_label['beach']
                size_counts, _ = np.histogram(time_selection_dict['size'][beach_selection], bins=size_bins,
                                              weights=time_selection_dict['particle_number'][beach_selection])
                output_dict['beach'][index_time] = size_counts
                # seabed particles
                seabed_selection = time_selection_dict['beach'] == beach_label['seabed']
                size_counts, _ = np.histogram(time_selection_dict['size'][seabed_selection], bins=size_bins,
                                              weights=time_selection_dict['particle_number'][seabed_selection])
                output_dict['seabed'][index_time] = size_counts
                # floating particles
                afloat_selection = time_selection_dict['beach'] == beach_label['afloat']
                size_counts, _ = np.histogram(time_selection_dict['size'][afloat_selection], bins=size_bins,
                                              weights=time_selection_dict['particle_number'][afloat_selection])
                output_dict['afloat'][index_time] = size_counts
                # floating particles within 5m of surface
                afloat5m_selection = (time_selection_dict['beach'] == beach_label['afloat']) & (time_selection_dict['z'] < 5)
                size_counts, _ = np.histogram(time_selection_dict['size'][afloat5m_selection], bins=size_bins,
                                              weights=time_selection_dict['particle_number'][afloat5m_selection])
                output_dict['afloat_5m'][index_time] = size_counts

    # Adding the index of the final timestep for ease later on
    output_dict['final_index'] = index_time

    # Saving the output
    prefix = 'size_distribution'
    output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
    utils.save_obj(output_name, output_dict)
    utils.print_statement("The size distribution has been saved")