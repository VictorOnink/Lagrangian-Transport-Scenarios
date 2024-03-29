from parcels import FieldSet
import numpy as np
from numpy import array
import xarray
import progressbar
from copy import deepcopy
import utils


def create_border_current(output_name: str, filenames: list, variables: dict, dimensions: dict, grid: array):
    # Setting up the fieldset containing the UV fields
    fieldset = set_fieldset(filenames, variables, dimensions)
    fieldset.computeTimeChunk(time=0, dt=-1)
    nx = fieldset.U.lon.size
    ny = fieldset.U.lat.size

    # U and V data
    u_data = reduce_array(fieldset.U.data)
    v_data = reduce_array(fieldset.V.data)

    # Creating the arrays that we will use for the first round of getting the border currents
    u_vel = np.zeros(u_data.shape)
    v_vel = np.zeros(v_data.shape)

    # Get the shore and coastal arrays
    shore = get_shore_cells(grid=grid)
    coastal = get_coastal_cells(grid=grid)

    # Looping through all the cells, checking if they are border currents or not
    for i in progressbar.progressbar(range(0, nx)):
        for j in range(1, ny - 1):
            if coastal[j, i] == 1:
                j_steps, i_steps = np.array([j - 1, j, j + 1]) % nx, np.array([i - 1, i, i + 1]) % ny
                grid_select = np.ix_(j_steps, i_steps)
                mask = land_borders(u_data[grid_select], v_data[grid_select])
                if not mask.all():
                    u_vel[j, i] = sum(mask[:, 2]) - sum(mask[:, 0])

    # Assure that all shore and coastal cells have border current values
    u_vel_all = deepcopy(u_vel)
    v_vel_all = deepcopy(v_vel)
    for i in progressbar.progressbar(range(1, nx - 1)):
        for j in range(1, ny - 1):
            if shore[j, i] == 1 or coastal[j, i] == 1:
                k = [-1, 1]
                u_vel_all[j, i] += (u_vel[j, (i + k[0]) % nx] + u_vel[j, (i + k[1]) % nx])
                v_vel_all[j, i] += (v_vel[(j + k[0]) % ny, i] + v_vel[(j + k[1]) % ny, i])
    # Carry out normalisation
    u_vel_all, v_vel_all = normalisation(u_vel_all, v_vel_all)

    # Creating the netCDF file
    coords = [('time', np.array([0])), ('lat', fieldset.U.lat), ('lon', fieldset.U.lon)]
    u_vel_xarray = xarray.DataArray(u_vel_all[np.newaxis, :, :], coords=coords)
    v_vel_xarray = xarray.DataArray(v_vel_all[np.newaxis, :, :], coords=coords)
    coord_dict = {'time': np.array([0]), 'lon': fieldset.U.lon, 'lat': fieldset.U.lat}
    dset = xarray.Dataset({'border_u': u_vel_xarray, 'border_v': v_vel_xarray}, coords=coord_dict)
    dset.to_netcdf(output_name)

    # Just to check some basics to see if it did what I want
    magnitude = np.sqrt(np.square(u_vel_all) + np.square(v_vel_all))
    assert np.max(magnitude, axis=(0, 1)) != 1, utils.print_statement("WARNING: The maximum magnitude is too high, namely {} m/s".format(np.max(magnitude, axis=(0, 1)) != 0))
    assert np.max(np.abs(u_vel_all), axis=(0, 1)) <= 1, utils.print_statement("WARNING: The maximum u component is too high, namely {} m/s".format(np.max(np.abs(u_vel_all), axis=(0, 1))))
    assert np.max(np.abs(v_vel_all), axis=(0, 1)) <= 1, utils.print_statement("WARNING: The maximum v component is too high, namely {} m/s".format(np.max(np.abs(v_vel_all), axis=(0, 1))))


def set_fieldset(filenames: list, variables: dict, dimensions: dict):
    filenames = {'U': filenames[0],
                 'V': filenames[0]}
    return FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=True)


def reduce_array(data: array):
    if len(data.shape) == 4:
        data = data[0, 0, :, :]
    elif len(data.shape) == 3:
        data = data[0, :, :]
    else:
        utils.print_statement("What weird data are you working with? It has shape {}".format(data.shape))
    return np.array(data)


def is_ocean(u: float, v: float):
    return u == 0 and v == 0


def land_borders(u: array, v: array):
    mask = np.ones((3, 3), dtype=bool)
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            mask[i, j] = is_ocean(u[i, j], v[i, j])
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
                        if k != 0 and m != 0:  # find if there is land directly adjacent
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
                        if k != 0 and m != 0:  # find if there is ocean directly adjacent
                            if not mask[(i + k) % mask.shape[0], (j + m) % mask.shape[1]]:
                                shore[i, j] = 1
    return shore


def normalisation(u: array, v: array):
    magnitude = np.sqrt(np.square(u) + np.square(v))
    non_zero = magnitude != 0
    u[non_zero] = np.divide(u[non_zero], magnitude[non_zero])
    v[non_zero] = np.divide(v[non_zero], magnitude[non_zero])
    return u, v
