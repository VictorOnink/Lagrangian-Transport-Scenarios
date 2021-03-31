import numpy as np
from numpy import array
import xarray
import progressbar
from copy import deepcopy
import os
import sys
from netCDF4 import Dataset
import settings
import utils


def create_MLD_files(UV_filenames: list, UV_variables: dict, TEMP_filenames: list, TEMP_variables: dict,
                     SALINITY_filenames: list, SALINITY_variables: dict, LON: array, LAT: array, DEPTH: array):
    """
    Computing the MLD based on the critical Richardson number approach adapted from van Roekel et al. (2018)
    """
    # The critical Richardson Number that delineates the MLD
    Ri_c = 0.3
    # Expanding the depth so that the array size is the same as the salinity fields
    DEPTH = np.tile(DEPTH[np.newaxis, :, np.newaxis, np.newaxis], (1, 1, LAT.shape[0], LON.shape[0]))
    # Creating a list to which we save all MLD filenames
    MLD_filenames = []
    if settings.ADVECTION_DATA is 'CMEMS_MEDITERRANEAN' and settings.START_YEAR is not 2010:
        UV_filenames = UV_filenames[1:]

    for step in range(len(UV_filenames)):
        # Loading the relevant UV, temperature and salinity fields
        UV_file, TEMP_file, SAL_file = UV_filenames[step], TEMP_filenames[step], SALINITY_filenames[step]
        print(UV_file)
        print(TEMP_file)
        print(SAL_file)
        print(' ')
        MLD_file = TEMP_file.replace('TEMP', 'MLD')
        MLD_filenames.append(MLD_file)

        if not utils._check_file_exist(MLD_file):
            UV_var, TEMP_var, SAL_var = [*UV_variables.keys()], [*TEMP_variables.keys()][0], [*SALINITY_variables.keys()][0]
            U, V = Dataset(UV_file).variables[UV_variables[UV_var[0]]][:], Dataset(UV_file).variables[
                                                                               UV_variables[UV_var[1]]][:]
            try:
                TEMP = Dataset(TEMP_file).variables[TEMP_variables[TEMP_var]][:]
            except:
                print('{} {}'.format(TEMP_file, TEMP_variables[TEMP_var]))
            try:
                SAL = Dataset(SAL_file).variables[SALINITY_variables[SAL_var]][:]
            except:
                print('{} {}'.format(SAL_file, SALINITY_variables[SAL_var]))
                sys.exit()
            TIME = Dataset(SAL_file).variables['time'][:]

            # Computing the buoyancy fields
            BUO = buoyancy_field(TEMP, SAL)

            # Computing the velocity shear
            SHEAR = velocity_shear(U, V)

            # Getting the Richardson Number
            Ri = richardson_number(BUO, SHEAR, DEPTH)

            # Determining which cells have Ri < Ri_c
            criteria = Ri > Ri_c

            # Getting the last depth at which Ri < Ri_c, which would be the MLD
            MLD = deepcopy(DEPTH)
            MLD[criteria] = np.nan
            MLD = np.nanmax(MLD, axis=(0, 1), keepdims=True)

            # Creating a NETCDF4 file containing the MLD field
            to_netcdf = MLD, TIME, LON, LAT, np.nanmin(DEPTH)
            create_netcdf(filename=MLD_file, to_netcdf=to_netcdf)
    return MLD_filenames


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
    # Velocity shear relative to surface
    magnitude = np.sqrt(np.square(U) + np.square(V))
    shear = magnitude - magnitude[:, 0, :, :]
    # We return the squared shear relative to the ocean surface, and add an added to term to prevent division by 0 at a
    # later point
    return np.square(shear) + 0.0001


def richardson_number(BUO, SHEAR, DEPTH):
    # Difference in buoyancy
    buo_diff = BUO[0, 0, :, :] - BUO
    # Ratio of buoyancy and shear
    ratio = np.divide(buo_diff, SHEAR)
    # Calculating the final richardson number
    Ri = np.multiply(ratio, DEPTH)
    return Ri


def create_netcdf(filename: str, to_netcdf: tuple):
    MLD, TIME, LON, LAT, DEPTH = to_netcdf
    root_grp = Dataset(filename, 'w', format='NETCDF4')
    # Create Dimensions
    root_grp.createDimension('time', TIME.size)
    root_grp.createDimension('lat', LAT.size)
    root_grp.createDimension('lon', LON.size)
    root_grp.createDimension('depth', 1)
    # variables
    time = root_grp.createVariable('time', 'f8', ('time',))
    time.units = 'seconds since 1970-01-01 00:00:00'
    time.calender = 'standard'
    lat = root_grp.createVariable('lat', 'f8', ('lat',))
    lon = root_grp.createVariable('lon', 'f8', ('lon',))
    depth = root_grp.createVariable('depth', 'f8', ('depth',))
    mld = root_grp.createVariable('MLD', 'f8', ('time', 'depth', 'lat', 'lon'))
    # Now assigning values to the variables
    lat[:], lon[:], time[:], mld[:], depth[:] = LAT, LON, TIME, MLD, DEPTH
    # Finally, close the new field
    root_grp.close()

