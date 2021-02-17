import numpy as np
from numpy import array
import xarray
import progressbar
from copy import deepcopy
import os
from netCDF4 import Dataset
import settings


def create_MLD_files(UV_filenames: list, UV_variables: dict, TEMP_filenames: list, TEMP_variables: dict,
                        SALINITY_filenames: list, SALINITY_variables: dict, lon: array, lat: array, depth: array):
    """
    Computing the MLD based on the critical Richardson number approach adapted from van Roekel et al. (2018)
    """
    for step in range(len(UV_filenames)):
        # Loading the relevant UV, temperature and salinity fields
        UV_file, TEMP_file, SAL_file = UV_filenames[step], TEMP_filenames[step], SALINITY_filenames[step]
        UV_var, TEMP_var, SAL_var = [*UV_variables.keys()], [*TEMP_variables.keys()][0], [*SALINITY_variables.keys()][0]
        U, V = Dataset(UV_file).variables[UV_variables[UV_var[0]]][:], Dataset(UV_file).variables[UV_variables[UV_var[1]]][:]
        TEMP = Dataset(TEMP_file).variables[TEMP_variables[TEMP_var]][:]
        SAL = Dataset(SAL_file).variables[SALINITY_variables[SAL_var]][:]

        # Computing the buoyancy fields
        BUO = buoyancy_field(TEMP, SAL)

        # Computing the velocity shear
        SHEAR = velocity_shear(U, V)
        print(SHEAR)


def buoyancy_field(TEMP, SAL):
    """
    Computing the buoyancy field, using the buoyancy definition according to eq. 28 from van Roekel et al. (2018).
    The buoyancy is calculated based on b = -g (1 - a (T - T_r) + b (S - S_r))
    a = 2e-4 K^-1
    b = 8e-4 ppt^-1 (PSU^-1)
    T_r = 298.15 K = 25 C
    S_r = 35 ppt (PSU)
    g = 9.81 m s^-2
    """
    g = settings.G
    a, b = settings.A_T, settings.B_S
    T_r, S_r = settings.T_R, settings.S_R
    return -g * (1 - a * (TEMP - T_r) + b * (SAL - S_r))


def velocity_shear(U, V):
    magnitude = np.sqrt(np.square(U) + np.square(V))
    shear = magnitude - magnitude[:, 0, :, :]
    return shear
