from parcels import FieldSet
import numpy as np
from numpy import array
import xarray
import progressbar
from copy import deepcopy
import os

def create_border_current(output_name: str, filenames: list, variables: dict, dimensions: dict, grid: array):
    # Setting up the fieldset containing the UV fields
    fieldset = set_fieldset(filenames, variables, dimensions)
    fieldset.computeTimeChunk(time=0, dt=-1)
    nx = fieldset.U.lon.size
    ny = fieldset.U.lat.size

    # Creating the arrays that we will use for the first round of getting the border currents
    u_vel = np.zeros(fieldset.U.data.shape)
    os.system('echo '+str(u_vel.shape))
    v_vel = np.zeros(fieldset.V.data.shape)

    # Looping through all the cells, checking if they are border currents or not
    for i in progressbar.progressbar(range(0, nx)):
        for j in range(1, ny - 1):
            if is_ocean(fieldset.U.data[0, 0, j, i], fieldset.V.data[0, 0, j, i]):
                mask = land_borders(fieldset.U.data[0, 0, :, :], fieldset.V.data[0, 0, :, :], j, i, nx)
                if not mask.all():
                    u_vel[0, j, i] = sum(mask[:, 2]) - sum(mask[:, 0])
                    v_vel[0, j, i] = sum(mask[2, :]) - sum(mask[0, :])
    # Get the shore and coastal arrays
    shore = get_shore_cells(grid=grid)
    coastal = get_coastal_cells(grid=grid)
    # Assure that all shore and coastal cells have border current values
    u_vel_all = deepcopy(u_vel)
    v_vel_all = deepcopy(v_vel)
    for i in progressbar.progressbar(range(1, nx - 1)):
        for j in range(1, ny - 1):
            if shore[j, i] == 1 or coastal[j, i] == 1:
                k = [-1, 1]
                u_vel_all[0, j, i] += (u_vel[0, j, (i + k[0]) % nx] + u_vel[0, j, (i + k[1]) % nx])
                v_vel_all[0, j, i] += (v_vel[0, (j + k[0]) % ny, i] + v_vel[0, (j + k[1]) % ny, i])
    # Carry out normalisation
    u_vel_all, v_vel_all = normalisation(u_vel_all, v_vel_all)

    # We don't care about the depth dimensions, we only want to keep lon and lat parts
    u_vel_all, v_vel_all = u_vel_all[0, :, :], v_vel_all[0, :, :]

    # Creating the netCDF file
    coords = [('lat', fieldset.U.lat), ('lon', fieldset.U.lon)]
    u_vel_xarray = xarray.DataArray(u_vel_all, coords=coords)
    v_vel_xarray = xarray.DataArray(v_vel_all, coords=coords)
    coord_dict = {'lon': fieldset.U.lon, 'lat': fieldset.U.lat}
    dset = xarray.Dataset({'border_u': u_vel_xarray, 'border_v': v_vel_xarray}, coords=coord_dict)
    dset.to_netcdf(output_name)


def set_fieldset(filenames: list, variables: dict, dimensions: dict):
    filenames = {'U': filenames[0],
                 'V': filenames[0]}
    return FieldSet.from_netcdf(filenames, variables, dimensions)


def is_ocean(u: float, v: float):
    return u == 0 and v == 0


def land_borders(u: array, v: array, m: int, k: int, nx: int):
    mask = np.ones((3, 3), dtype=bool)
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            mask[j + 1, i + 1] = is_ocean(u[j + m, (i + k) % nx], v[j + m, (i + k) % nx])
    return mask


def get_coastal_cells(grid: array):
    # Going through ocean cells to see which are next to land
    mask = grid.mask
    coastal = np.zeros(grid.shape)
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            if not mask[i, j]:
                for k in [-1, 0, 1]:
                    for m in [-1, 0, 1]:
                        if abs(k) != abs(m):  # find if there is land directly adjacent
                            if mask[(i + k) % mask.shape[0], (j + m) % mask.shape[1]]:
                                coastal[i, j] = 1
    return coastal


def get_shore_cells(grid: array):
    # Going through land cells to see which are next to the ocean
    mask = grid.mask
    shore = np.zeros(mask.shape)
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            if mask[i, j]:
                for k in [-1, 0, 1]:
                    for m in [-1, 0, 1]:
                        if abs(k) != abs(m):  # find if there is ocean directly adjacent
                            if not mask[(i + k) % mask.shape[0], (j + m) % mask.shape[1]]:
                                shore[i, j] = 1
    return shore


def normalisation(u: array, v: array):
    magnitude = np.sqrt(np.square(u) + np.square(v))
    non_zero = magnitude != 0
    u[non_zero] = np.divide(u[non_zero], magnitude[non_zero])
    v[non_zero] = np.divide(v[non_zero], magnitude[non_zero])
    return u, v
