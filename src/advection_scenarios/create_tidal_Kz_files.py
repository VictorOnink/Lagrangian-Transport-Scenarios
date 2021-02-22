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
    de Lavergne et al. (2020)
    """
    # Loading the global data from de Lavergne et al. (2020)
    TIDAL_filename = settings.DATA_INPUT_DIR_SERVERS + 'global_tidal_energy_dissipation.nc'
    TIDAL_data = Dataset(TIDAL_filename)
    TIDAL_fields = {'LON_T': TIDAL_data.variables['lon'][:], 'LAT_T': TIDAL_data.variables['lat'][:], }