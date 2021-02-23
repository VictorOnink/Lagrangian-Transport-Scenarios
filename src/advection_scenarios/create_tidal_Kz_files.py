import numpy as np
from numpy import array
import xarray
import progressbar
from copy import deepcopy
import os
from netCDF4 import Dataset
import settings
import utils

def create_tidal_Kz_files(LON: array, LAT: array, DEPTH: array, BATH_filenames: list, BATH_variables: dict):
    """
    Determining Kz from the estimates of production of turbulent kinetic energy, as outlined in
    de Lavergne et al. (2020) and based on code shared by Clement Vic
    """
    # Loading the global data from de Lavergne et al. (2020)
    TIDAL_filename = settings.DATA_INPUT_DIR_SERVERS + 'global_tidal_energy_dissipation.nc'
    TIDAL_data = {}
    for key in Dataset(TIDAL_filename).variables.keys:
        TIDAL_data[key] = Dataset(TIDAL_filename).variables[key][:]

    TIDAL_data[10000]