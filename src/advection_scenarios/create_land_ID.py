import numpy as np
import xarray
import progressbar


def create_land_ID(output_name: str, grid: np.array, lon: np.array, lat: np.array):
    land_ID = np.zeros(grid.shape)
    mask = grid.mask

    for lat_index in progressbar.progressbar(range(grid.shape[0])):
        for lon_index in range(grid.shape[1]):
            if not mask[lat_index, lon_index]:
                for lat_step in [-1, 0]:
                    for lon_step in [-1, 0]:
                        if abs(lat_step) != abs(lon_step):
                            if mask[(lat_index + lat_step) % grid.shape[0], (lon_index + lon_step) % grid.shape[1]]:
                                land_ID[lat_index, lon_index] = 1
    # Saving the entire land ID field
    coords = [('lat', lat), ('lon', lon)]
    id = xarray.DataArray(land_ID, coords=coords)
    dcoo = {'lat': lat, 'lon': lon}
    dset = xarray.Dataset({'land_ID': id}, coords=dcoo)
    dset.to_netcdf(output_name)
