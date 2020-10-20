import numpy as np
import xarray
import progressbar
import geopy.distance
import settings


def create_distance_to_shore(output_name: str, grid: np.array, lon: np.array, lat: np.array):
    # Getting the dimensions of the model grid
    mask = grid.mask
    Lon, Lat = np.meshgrid(lon, lat)
    n_lon, n_lat = len(lon), len(lat)

    # Initializing the distance array
    distance = np.zeros(mask.shape)

    # Creating a memory variable that keeps track of how many cells away land is. That way, when we move from one
    # cell to it's neighbor, then we don't need to check everything again...
    memory_var = 1

    # The actual distance computation loop
    for lat_index in progressbar.progressbar(range(mask.shape[0])):
        for lon_index in range(mask.shape[1]):
            if not mask[lat_index, lon_index]:
                if memory_var > 2:
                    cells = memory_var - 2
                else:
                    cells = 1
                land_lon, land_lat, dis = [], [], []
                while len(land_lon) == 0:
                    for lat_step in [-cells, cells]:
                        for lon_step in range(-cells, cells + 1):
                            if boundary_conditions(n_lat, n_lon, lat_index, lat_step, lon_index, lon_step):
                                if mask[(lat_index + lat_step) % n_lat, (lon_index + lon_step) % n_lon]:
                                    land_lat.append((lat_index + lat_step) % n_lat)
                                    land_lon.append((lon_index + lon_step) % n_lon)
                    for lat_step in range(-cells, cells + 1):
                        for lon_step in [-cells, cells]:
                            if boundary_conditions(n_lat, n_lon, lat_index, lat_step, lon_index, lon_step):
                                if mask[(lat_index + lat_step) % n_lat, (lon_index + lon_step) % n_lon]:
                                    land_lat.append((lat_index + lat_step) % n_lat)
                                    land_lon.append((lon_index + lon_step) % n_lon)
                    cells += 1
                memory_var = cells

                for points in range(len(land_lon)):
                    dis.append(geopy.distance.distance((lat[lat_index], lon[lon_index]),
                                                       (lat[land_lat[points]], lon[land_lon[points]])).km)
                distance[lat_index, lon_index] += np.min(dis)
            else:
                memory_var = 1
    # Saving the entire distance field
    coords = [('lat', lat), ('lon', lon)]
    dist = xarray.DataArray(distance, coords=coords)
    dcoo = {'lat': lat, 'lon': lon}
    dset = xarray.Dataset({'distance': dist}, coords=dcoo)
    dset.to_netcdf(output_name)

def create_distance_to_shore_land(output_name: str, grid: np.array, lon: np.array, lat: np.array):
    # Getting the dimensions of the model grid
    mask = grid.mask
    Lon, Lat = np.meshgrid(lon, lat)
    n_lon, n_lat = len(lon), len(lat)

    # Initializing the distance array
    distance = np.zeros(mask.shape)

    # Creating a memory variable that keeps track of how many cells away land is. That way, when we move from one
    # cell to it's neighbor, then we don't need to check everything again...
    memory_var = 1

    # The actual distance computation loop
    for lat_index in progressbar.progressbar(range(mask.shape[0])):
        for lon_index in range(mask.shape[1]):
            if mask[lat_index, lon_index]:
                if memory_var > 2:
                    cells = memory_var - 2
                else:
                    cells = 1
                land_lon, land_lat, dis = [], [], []
                while len(land_lon) == 0:
                    for lat_step in [-cells, cells]:
                        for lon_step in range(-cells, cells + 1):
                            if boundary_conditions(n_lat, n_lon, lat_index, lat_step, lon_index, lon_step):
                                if not mask[(lat_index + lat_step) % n_lat, (lon_index + lon_step) % n_lon]:
                                    land_lat.append((lat_index + lat_step) % n_lat)
                                    land_lon.append((lon_index + lon_step) % n_lon)
                    for lat_step in range(-cells, cells + 1):
                        for lon_step in [-cells, cells]:
                            if boundary_conditions(n_lat, n_lon, lat_index, lat_step, lon_index, lon_step):
                                if not mask[(lat_index + lat_step) % n_lat, (lon_index + lon_step) % n_lon]:
                                    land_lat.append((lat_index + lat_step) % n_lat)
                                    land_lon.append((lon_index + lon_step) % n_lon)
                    cells += 1
                memory_var = cells

                for points in range(len(land_lon)):
                    dis.append(geopy.distance.distance((lat[lat_index], lon[lon_index]),
                                                       (lat[land_lat[points]], lon[land_lon[points]])).km)
                distance[lat_index, lon_index] += np.min(dis)
            else:
                memory_var = 1
    # Saving the entire distance field
    coords = [('lat', lat), ('lon', lon)]
    dist = xarray.DataArray(distance, coords=coords)
    dcoo = {'lat': lat, 'lon': lon}
    dset = xarray.Dataset({'distance': dist}, coords=dcoo)
    dset.to_netcdf(output_name)

def boundary_conditions(n_lat: int, n_lon: int, lat_index: int, k: int, lon_index: int, m: int):
    if settings.ADVECTION_DATA == 'HYCOM_GLOBAL':
        return (lat_index + k) < n_lat
    if settings.ADVECTION_DATA == 'HYCOM_CARIBBEAN':
        return ((lat_index + k) < n_lat) & ((lat_index + k) >= 0) & ((lon_index + m) < n_lon) & ((lon_index + m) >= 0)
