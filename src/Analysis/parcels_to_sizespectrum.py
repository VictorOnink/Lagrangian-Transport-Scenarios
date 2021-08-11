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
    size_bins = np.logspace(start=-6, stop=-2, num=50)

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

            # Calculating the spectrum every 30 days
            for index_time, time_select in enumerate(range(0, len(time_list), 60)):
                size_selection = size[time_select == time]
                size_counts, _ = np.histogram(size_selection, bins=size_bins)
                output_dict[index_time] = size_counts

    # Saving the output
    prefix = 'size_distribution'
    output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
    utils.save_obj(output_name, output_dict)
    utils.print_statement("The timeseries has been saved")
