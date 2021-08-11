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

    # Loading the time axis
    time_list = np.array([], dtype=np.int)
    for restart in range(settings.SIM_LENGTH):
        parcels_file = file_dict[0][restart]
        parcels_dataset = Dataset(parcels_file)
        time = parcels_dataset.variables['time'][:, :-1]
        for time_case in np.unique(time):
            if np.nansum(time_case == time) > 10 and time_case != np.nan:
                time_list = np.append(time_list, time_case)
    utils.print_statement(len(time_list), to_print=True)