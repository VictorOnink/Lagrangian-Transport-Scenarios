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
        time_list = np.append(time_list, np.unique(parcels_dataset.variables['time'][:, :-1]))
    utils.print_statement('we have {} time points'.format(time_list.size))