import numpy as np
import xarray
import geopy.distance
import settings
import utils


def create_distance_to_shore(output_name: str, grid: np.array, lon: np.array, lat: np.array):
    # Getting the dimensions of the model grid
    land = grid.mask
    n_lon, n_lat = len(lon), len(lat)

    # Initializing the distance array
    distance = np.zeros(land.shape)

    # Creating a search_memory variable that keeps track of how many cells away land is. That way, when we move from one
    # cell to it's neighbor, then we don't need to check everything again...
    search_memory = 1

    # Looping through all the cells
    for lat_index in range(land.shape[0]):
        for lon_index in range(land.shape[1]):
            # If the land mask is false, find the nearest land cell
            if not land[lat_index, lon_index]:
                # Reduce the search_memory from the previous cell by two. This saves computational effort because if for
                # example you are in the middle of the Pacific and you know that for the neighboring cell the nearest
                # land cell was 500 cells away, you don't need to check all those cells again for the current cell since
                # you've only moved one cell.
                if search_memory > 2:
                    box_size = search_memory - 2
                else:
                    box_size = 1

                # These lists keep track of all land cells that are encountered in the search
                land_lon, land_lat, land_distance = [], [], []

                # Keep looking for land cells until you encounter at least one.
                while land_lon.__len__() == 0:
                    # First check the top and bottom rows of cells in the search box
                    for lat_step in [-box_size, box_size]:
                        for lon_step in range(-box_size, box_size + 1):
                            # Check if any boundary conditions are being broken
                            if boundary_conditions(n_lat, n_lon, lat_index, lat_step, lon_index, lon_step):
                                # If the cell being checked is land, save the lon and lat indices
                                if land[(lat_index + lat_step) % n_lat, (lon_index + lon_step) % n_lon]:
                                    land_lat.append((lat_index + lat_step) % n_lat)
                                    land_lon.append((lon_index + lon_step) % n_lon)
                    # Then check the left and right columns of cells in the search box.
                    for lat_step in range(-box_size, box_size + 1):
                        for lon_step in [-box_size, box_size]:
                            # Check if any boundary conditions are being broken
                            if boundary_conditions(n_lat, n_lon, lat_index, lat_step, lon_index, lon_step):
                                # If the cell being checked is land, save the lon and lat indices
                                if land[(lat_index + lat_step) % n_lat, (lon_index + lon_step) % n_lon]:
                                    land_lat.append((lat_index + lat_step) % n_lat)
                                    land_lon.append((lon_index + lon_step) % n_lon)
                    # If we don't encounter land cells, we increase the size of the search box by 1.
                    box_size += 1

                # Once we have found land, save the size of the search box in the memory
                search_memory = box_size

                # For all the encountered land cells, determine the distance to the ocean cell in kilometers
                for points in range(land_lon.__len__()):
                    land_distance.append(geopy.distance.distance((lat[lat_index], lon[lon_index]),
                                                       (lat[land_lat[points]], lon[land_lon[points]])).km)

                # Update the distance array with the shortest distance to land in the land_distance list
                distance[lat_index, lon_index] = np.min(land_distance)
            else:
                # If the cell we are looking at is land, set the search_memory to 1
                search_memory = 1

    # Saving the entire distance field
    dset = xarray.Dataset({"distance": xarray.DataArray(distance, coords=[("lat", lat), ("lon", lon)])},
                          coords={"lat": lat, "lon": lon})
    dset.to_netcdf(output_name)


def create_distance_to_shore_land(output_name: str, grid: np.array, lon: np.array, lat: np.array):
    # Getting the dimensions of the model grid
    # utils.print_statement("grid shape {} {}".format(grid.shape[0], grid.shape[1]))
    mask = grid.mask
    n_lat, n_lon = mask.shape
    # Initializing the distance array
    distance = np.zeros(mask.shape)
    print("distance shape {}".format(distance.shape[0], distance.shape[1]))
    # Creating a memory variable that keeps track of how many cells away land is. That way, when we move from one
    # cell to it's neighbor, then we don't need to check everything again...
    memory_var = 1
    # The actual distance computation loop
    for lat_index in range(n_lat):
        for lon_index in range(n_lon):
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

    distance[mask == False] = 0

    # Saving the entire distance field
    coords = [("time", np.array([0])), ("lat", lat), ("lon", lon)]
    dist = xarray.DataArray(distance[np.newaxis, :, :], coords=coords)
    dcoo = {"time": np.array([0]), "lat": lat, "lon": lon}
    dset = xarray.Dataset({"distance": dist}, coords=dcoo)
    dset.to_netcdf(output_name)


def boundary_conditions(n_lat: int, n_lon: int, lat_index: int, k: int, lon_index: int, m: int):
    if settings.ADVECTION_DATA == "HYCOM_GLOBAL":
        return (lat_index + k) < n_lat & (lat_index + k) > 0
    if settings.ADVECTION_DATA == "HYCOM_CARIBBEAN":
        return ((lat_index + k) < n_lat) & ((lat_index + k) >= 0) & ((lon_index + m) < n_lon) & ((lon_index + m) >= 0)
    if settings.ADVECTION_DATA == "CMEMS_MEDITERRANEAN":
        return ((lat_index + k) < n_lat) & ((lat_index + k) >= 0) & ((lon_index + m) < n_lon) & ((lon_index + m) >= 0)
