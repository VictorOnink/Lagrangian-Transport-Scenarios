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

    # A number of physical parameters
    gamma = 0.2  # Mixing efficiency

    # The TIDAL_data gridding isn't regular in the z-direction. We will first interpolate the epsilon_tid and
    # buoyancy_frequency_squared fields onto the DEPTH levels
    TIDAL_inter = interpolate_to_DEPTH(TIDAL_data, DEPTH)

    TIDAL_data[10000]


def interpolate_to_DEPTH(TIDAL_data, DEPTH):
    depth_midpoint = TIDAL_data['depth_midpoint']
    TIDAL_inter = {}
    for key in ['buoyancy_frequency_squared', 'epsilon_tid']:
        field = TIDAL_data[key]
        field[field.mask] = 0
        field_inter = np.zeros(field.shape)
        for lat in range(field.shape[1]):
            for lon in range(field.shape[2]):
                inter_f = interpolate.interp1d(depth_midpoint[:, lat, lon], field[:, lat, lon], bounds_error=False)
                field_inter[:, lat, lon] = inter_f(DEPTH)
        TIDAL_inter[key] = field_inter
    return TIDAL_inter
