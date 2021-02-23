import numpy as np
from numpy import array
import scipy.interpolate as interpolate
import xarray
import progressbar
from copy import deepcopy
import os
from netCDF4 import Dataset
import settings
import utils

def create_tidal_Kz_files(LON: array, LAT: array, DEPTH: array, BATH_filenames: str, BATH_variables: dict):
    """
    Determining Kz from the estimates of production of turbulent kinetic energy, as outlined in
    de Lavergne et al. (2020) and based on code shared by Clement Vic
    """
    # Loading the global data from de Lavergne et al. (2020)
    TIDAL_filename = utils._get_input_directory(server=settings.SERVER) + 'global_tidal_energy_dissipation.nc'
    TIDAL_data = {}
    for key in Dataset(TIDAL_filename).variables.keys():
        TIDAL_data[key] = Dataset(TIDAL_filename).variables[key][:]

    # Load the bathymetry data
    bathymetry = Dataset(BATH_filenames).variables[BATH_variables['DEPTH']][:]
    bathymetry[bathymetry.mask] = 0

    # Computing Kz on the TIDAL_data grid according to Kv = gamma * epsilon / N^2
    gamma = 0.2  # Mixing efficiency
    TIDAL_Kz = np.divide(gamma * TIDAL_data['epsilon_tid'], TIDAL_data['buoyancy_frequency_squared'])

    # The TIDAL_data gridding isn't regular in the z-direction. We will first interpolate the TIDAL_Kz fields onto the
    # DEPTH levels
    TIDAL_Kz_inter = interpolate_to_DEPTH(TIDAL_Kz, TIDAL_data, DEPTH)

    TIDAL_data[10000]


def interpolate_to_DEPTH(TIDAL_Kz: array, TIDAL_data: dict, DEPTH: array):
    TIDAL_depth = TIDAL_data['depth_midpoint']
    TIDAL_depth[TIDAL_depth.mask] = np.nanmax(TIDAL_depth)
    TIDAL_Kz[TIDAL_Kz.mask] = 0
    TIDAL_Kz_inter = np.zeros((DEPTH.shape[0], TIDAL_Kz.shape[1], TIDAL_Kz.shape[2]))
    for lat in range(TIDAL_Kz.shape[1]):
        for lon in range(TIDAL_Kz.shape[2]):
            inter_f = interpolate.interp1d(TIDAL_depth[:, lat, lon], TIDAL_Kz[:, lat, lon], bounds_error=False)
            TIDAL_Kz_inter[:, lat, lon] = inter_f(np.array(DEPTH))
    return TIDAL_Kz_inter
