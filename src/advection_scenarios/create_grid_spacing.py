import settings
import os
import numpy as np
import xarray
import progressbar


def create_grid_spacing(output_name: str, grid: np.array, lon: np.array, lat: np.array):
    grid_lon_spacing = np.zeros(grid.shape)
    grid_lat_spacing = np.zeros(grid.shape)

    for lat_step in range(grid.shape[0] - 1):
        for lon_step in range(grid.shape[1] - 1):
            grid_lon_spacing[lat_step, lon_step] = np.abs(lon[lon_step] - lon[lon_step + 1])
            grid_lat_spacing[lat_step, lon_step] = np.abs(lat[lat_step] - lat[lat_step + 1])

    grid_lon_spacing = fill_last(grid_lon_spacing)
    grid_lat_spacing = fill_last(grid_lat_spacing)

    # Saving the entire grid spacing fields
    coords = [('lat', lat), ('lon', lon)]
    lon_space = xarray.DataArray(grid_lon_spacing, coords=coords)
    lat_space = xarray.DataArray(grid_lat_spacing, coords=coords)
    dcoo = {'lat': lat, 'lon': lon}
    dset = xarray.Dataset({'lon_spacing': lon_space, 'lat_spacing': lat_space}, coords=dcoo)
    dset.to_netcdf(output_name)

    # Checks to see if the grid spacing calculation works as expected
    os.system('echo "The maximum lon spacing is {}, and the minimum is {}"'.format(grid_lon_spacing.max(),
                                                                                   grid_lon_spacing.min()))
    os.system('echo "The maximum lat spacing is {}, and the minimum is {}"'.format(grid_lat_spacing.max(),
                                                                                   grid_lat_spacing.min()))


def fill_last(array: np.array):
    array[:, -1] = array[:, -2]
    array[-1, :] = array[-2, :]
    return array
