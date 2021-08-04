import settings
import utils
from advection_scenarios.create_distance_to_shore import create_distance_to_shore_land

import numpy as np
import glob
import os
from netCDF4 import Dataset
import pandas as pd
import geopy.distance
from datetime import timedelta
import math
import fiona
from shapely import vectorized
from shapely.geometry import shape
import xarray
import progressbar
from parcels import Field


def create_input_files(prefix: str, grid: np.array, lon: np.array, lat: np.array, repeat_dt):
    # Create the output prefix and then check if any files with such prefixes exist
    if settings.INPUT in ['Jambeck', 'Lebreton']:
        output_prefix = settings.INPUT_DIREC + settings.INPUT + '_{}_{}_'.format(prefix, settings.START_YEAR)
    elif settings.INPUT == 'Point_Release':
        str_format = (prefix, settings.START_YEAR, settings.INPUT_LAT, settings.INPUT_LON)
        output_prefix = settings.INPUT_DIREC + settings.INPUT + '_{}_{}_{}_{}_'.format(*str_format)
    elif settings.INPUT == 'Uniform':
        str_format = (prefix, settings.START_YEAR)
        output_prefix = settings.INPUT_DIREC + settings.INPUT + '_{}_{}_'.format(*str_format)
    else:
        utils.print_statement("Perhaps take another look at what input you are using?")

    # Check existence
    if len(glob.glob(output_prefix + '*')) > 0:
        utils.print_statement("The input files {} are already present".format(output_prefix))
    else:
        utils.print_statement("We need to create the input files {}".format(output_prefix))
        # Calculating the number of particle releases per year
        releases = number_of_releases(repeat_dt)
        if settings.INPUT == 'Jambeck':
            # Get the population data
            dataset = Dataset(settings.INPUT_DIREC + 'gpw_v4_population_count_adjusted_rev11_2pt5_min.nc')
            var_name = 'UN WPP-Adjusted Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'
            population = np.array(dataset.variables[var_name][2, :, :])
            lon_population = dataset.variables['longitude'][:]
            lat_population = dataset.variables['latitude'][:]
            Lon_population, Lat_population = np.meshgrid(lon_population, lat_population)
            # Get the mismanaged plastic fraction
            mismanaged = get_mismanaged_fraction_Jambeck(dataset=dataset)
            # Get the distance from land to shore
            distance_file = settings.INPUT_DIREC + prefix + '_distance_to_coast_land.nc'
            distance = get_distance_to_shore(filename=distance_file, grid=grid, lon=lon, lat=lat)
            # The yearly mismanaged plastic
            transport_to_ocean = 0.15  # percentage of mismanaged plastic that reaches the ocean
            kg_to_tons = 1000
            mismanaged_total = np.multiply(mismanaged, population) * 365 * transport_to_ocean / kg_to_tons
            # Get everything in column arrays to load
            lon_inputs = Lon_population[mismanaged_total > 0].flatten()
            lat_inputs = Lat_population[mismanaged_total > 0].flatten()
            plastic_inputs = mismanaged_total[mismanaged_total > 0].flatten()
        elif settings.INPUT == 'Lebreton':
            utils.print_statement("Loading the lebreton data")
            lebData = pd.read_csv(settings.INPUT_DIREC + 'PlasticRiverInputs.csv')
            lon_inputs, lat_inputs = np.array(lebData['X']), np.array(lebData['Y'])
            plastic_inputs = np.array(lebData['i_low'])
        elif settings.INPUT == 'Point_Release':
            lon_inputs, lat_inputs = np.ones((1, 1)) * settings.INPUT_LON, np.ones((1, 1)) * settings.INPUT_LAT
            plastic_inputs = np.ones((1, 1)) * settings.INPUT_MAX
        elif settings.INPUT == 'Uniform':
            lon_min, lon_max, lat_min, lat_max = np.min(lon), np.max(lon), np.min(lat), np.max(lat)
            release_grid = np.mgrid[lon_min:lon_max:settings.RELEASE_GRID, lat_min:lat_max:settings.RELEASE_GRID]
            n = release_grid[0].size
            lon_inputs, lat_inputs = np.reshape(release_grid[0], n), np.reshape(release_grid[1], n)
            # Remove cells on land
            Land = land_mark_checker(grid, lon, lat)
            [lon_inputs, lat_inputs] = [
                np.array([lo for lo, la in zip(lon_inputs, lat_inputs) if Land[0, 0, la, lo] == 0.0]),
                np.array([la for lo, la in zip(lon_inputs, lat_inputs) if Land[0, 0, la, lo] == 0.0])]
            # If we are dealing with a uniform release, we are going to assume that weight doesn't really matter, and
            # we only care about the release positions.
            utils.print_statement('The number of particles released per year is {} particles'.format(len(lon_inputs) * releases))
            split_to_runs(particle_lat=lat_inputs, particle_lon=lon_inputs, particle_weight=None,
                          output_prefix=output_prefix)
            return output_prefix
        # Only keep the particles that are within the domain
        lon_inputs, lat_inputs, plastic_inputs = within_domain(lon=lon, lat=lat, lon_inputs=lon_inputs,
                                                               lat_inputs=lat_inputs, plastic_inputs=plastic_inputs)

        # Get the inputs onto the grid of the advection data
        inputs_grid = utils.histogram(lon_data=lon_inputs, lat_data=lat_inputs, bins_Lon=lon, bins_Lat=lat,
                                      weight_data=plastic_inputs, area_correc=False)
        if settings.INPUT == 'Jambeck':
            inputs_grid[distance > 50] = 0
        # Get the ocean cells adjacent to coastal ocean cells, as these are the ones in which the particles will
        # be placed. For brevity, we will refer to these cells as coastal in the code here
        coastal = get_coastal_extended(grid=grid, coastal=get_coastal_cells(grid=grid))
        # Getting the inputs onto the coastal cells, and then computing the inputs per release_dt
        inputs_coastal_grid = input_to_nearest_coastal(coastal=coastal, inputs_grid=inputs_grid, lon=lon, lat=lat)
        inputs_coastal_grid /= releases
        # How many particles will be released per cell, and what are the assigned weights
        particle_number, particle_weight, particle_remain_num, particle_remain = number_weights_releases(
            input_grid=inputs_coastal_grid)
        # Get arrays of all particles release points and associated weights
        particle_lat, particle_lon, particle_weight = particle_grid_to_list(particle_number=particle_number,
                                                                            particle_weight=particle_weight,
                                                                            particle_remain_num=particle_remain_num,
                                                                            particle_remain=particle_remain,
                                                                            lon=lon, lat=lat)
        # How much of the plastic is being accounted for
        missing_percent = np.divide(np.sum(inputs_coastal_grid) - np.sum(particle_weight),
                                    np.sum(inputs_coastal_grid)) * 100
        utils.print_statement("The particles account for {}% of the total inputs".format(100 - missing_percent))
        str_format = len(particle_lat), len(particle_lat) * releases, releases
        utils.print_statement("We release {} particles per release step, so {} per year over {} steps".format(*str_format))
        # Dividing the particles into runs
        split_to_runs(particle_lat=particle_lat, particle_lon=particle_lon, particle_weight=particle_weight,
                      output_prefix=output_prefix)
        utils.print_statement("The input files have been created")
    # Returning the output prefix
    return output_prefix


def get_mismanaged_fraction_Jambeck(dataset: Dataset):
    mismanaged_file = settings.INPUT_DIREC + 'Jambeck_mismanaged_grid.nc'
    if utils.check_file_exist(mismanaged_file):
        utils.print_statement("The mismanaged grid already exists")
        return Dataset(mismanaged_file).variables['mismanaged_plastic'][:]
    else:
        utils.print_statement("We need to generate the mismanaged grid")
        # Load the grid of the population data
        lon_pop, lat_pop = dataset.variables['longitude'][:], dataset.variables['latitude'][:]
        Lon, Lat = np.meshgrid(lon_pop, lat_pop)
        # Load the Jambeck estimates of mismanaged plastic per capita
        jambeck_excel = pd.read_excel(settings.INPUT_DIREC + 'Jambeck2010data.xlsx')
        jambeck_country = list(jambeck_excel['Country'])
        jambeck_data = jambeck_excel['Mismanaged plastic waste [kg/person/day]7']
        # Initialize grid for mismanaged plastic estimates
        mismanaged_grid = np.zeros(Lon.shape)
        # Getting the country shapefiles
        countries = fiona.open(settings.INPUT_DIREC +
                               'country_shapefile/gpw_v4_national_identifier_grid_rev11_30_sec.shp')
        # Looping through all the countries, and marking which of the grid cells fall within which
        for country_index in progressbar.progressbar(range(len(countries))):
            if countries[country_index]['properties']['NAME0'] in jambeck_country:
                country_geometry = shape(countries[country_index]['geometry'])
                country_mask = vectorized.contains(country_geometry, Lon, Lat)
                if country_index in [145, 146]:
                    # 122 is the Netherlands Antilles value
                    mismanaged_grid[country_mask] = jambeck_data[122]
                else:
                    mismanaged_grid[country_mask] = jambeck_data[
                        jambeck_country.index(countries[country_index]['properties']['NAME0'])]
        # Saving the entire distance field
        utils.print_statement("Starting to save the mismanaged grid")
        coords = [('lat', lat_pop), ('lon', lon_pop)]
        misman = xarray.DataArray(mismanaged_grid, coords=coords)
        dcoo = {'lat': lat_pop, 'lon': lon_pop}
        dset = xarray.Dataset({'mismanaged_plastic': misman}, coords=dcoo)
        dset.to_netcdf(mismanaged_file)
        utils.print_statement("The mismanaged grid has now been created and saved for future use")
        return mismanaged_grid


def get_distance_to_shore(filename: str, grid: np.array, lon: np.array, lat: np.array):
    if utils.check_file_exist(filename):
        utils.print_statement("The mismanaged grid already exists")
    else:
        create_distance_to_shore_land(output_name=filename, grid=grid, lon=lon, lat=lat)
    return slicing_correction(Dataset(filename).variables['distance'][:, :])


def slicing_correction(array: np.array):
    if len(array.shape) == 2:
        return array
    elif len(array.shape) == 3:
        return array[0, :, :]


def within_domain(lon: np.array, lat: np.array, lon_inputs: np.array, lat_inputs: np.array, plastic_inputs: np.array):
    lon_max, lon_min = np.max(lon), np.min(lon)
    lat_max, lat_min = np.max(lat), np.min(lat)
    # Check which cells are within the domain and non-zero inputs
    domain = (lon_inputs <= lon_max) & (lon_inputs >= lon_min) & (lat_inputs >= lat_min) & \
             (lat_inputs <= lat_max) & (plastic_inputs > 0)
    str_format = (len(lon_inputs), settings.INPUT, np.sum(domain * 1))
    utils.print_statement("Of the original {} input sites in the {} scenario, {} are within the domain".format(*str_format))
    return lon_inputs[domain], lat_inputs[domain], plastic_inputs[domain]


def get_coastal_cells(grid: np.array):
    mask = grid.mask
    coastal = np.zeros(mask.shape, dtype=bool)
    for lat_index in range(mask.shape[0]):
        for lon_index in range(mask.shape[1]):
            if not mask[lat_index, lon_index]:
                for lat_step in [-1, 0, 1]:
                    for lon_step in [-1, 0, 1]:
                        if mask[(lat_index + lat_step) % mask.shape[0], (lon_index + lon_step) % mask.shape[1]]:
                            coastal[lat_index, lon_index] = True
    return coastal


def get_coastal_extended(grid: np.array, coastal: np.array):
    mask = grid.mask
    # We want the next layer of ocean cells after the coastal cells
    coastal_extended = np.zeros(mask.shape, dtype=bool)
    for lat_index in range(mask.shape[0]):
        for lon_index in range(mask.shape[1]):
            if not mask[lat_index, lon_index] and coastal[lat_index, lon_index]:
                for lat_step in [-1, 0, 1]:
                    for lon_step in [-1, 0, 1]:
                        if coastal[(lat_index + lat_step) % mask.shape[0], (lon_index + lon_step) % mask.shape[1]]:
                            coastal_extended[lat_index, lon_index] = True
    return coastal_extended


def input_to_nearest_coastal(coastal: np.array, inputs_grid: np.array, lon: np.array, lat: np.array):
    # Find the positions with non-zero inputs
    non_zero_inputs = np.where(inputs_grid > 0)
    inputs_coastal_grid = np.zeros(coastal.shape)
    # Looping through the non-zero points to find the nearest one with
    for point in range(len(non_zero_inputs[0])):
        lat_point, lon_point = non_zero_inputs[0][point], non_zero_inputs[1][point]
        # If the input is already in the coastal cell
        if coastal[lat_point, lon_point]:
            inputs_coastal_grid[lat_point, lon_point] += inputs_grid[lat_point, lon_point]
        # Else find the nearest coastal cell
        else:
            step = 1
            coastal_lon, coastal_lat, coastal_distance = [], [], []
            while len(coastal_lon) < 1:
                for lat_step in [-step, step]:
                    for lon_step in range(-step, step + 1):
                        if coastal[
                            (lat_point + lat_step) % coastal.shape[0], (lon_point + lon_step) % coastal.shape[1]]:
                            coastal_lat.append((lat_point + lat_step) % coastal.shape[0])
                            coastal_lon.append((lon_point + lon_step) % coastal.shape[1])
                step += 1
            # Now find the coastal cell that is geographically the closest to the input
            nearest_dist = geopy.distance.distance((lat[lat_point], lon[lon_point]),
                                                   (lat[coastal_lat[0]], lon[coastal_lon[0]]))
            nearest_lat, nearest_lon = coastal_lat[0], coastal_lon[0]
            if len(coastal_lon) > 0:
                for candidate in range(len(coastal_lon)):
                    distance = geopy.distance.distance((lat[lat_point], lon[lon_point]),
                                                       (lat[coastal_lat[candidate]], lon[coastal_lon[candidate]]))
                    if distance < nearest_dist:
                        nearest_dist = distance
                        nearest_lat, nearest_lon = coastal_lat[candidate], coastal_lon[candidate]
            inputs_coastal_grid[nearest_lat, nearest_lon] += inputs_grid[lat_point, lon_point]
    return inputs_coastal_grid


def number_weights_releases(input_grid: np.array):
    particle_number = np.zeros(input_grid.shape)
    particle_weight = np.zeros(input_grid.shape)
    particle_remain = np.zeros(input_grid.shape)
    particle_remain_num = np.zeros(input_grid.shape)
    # First all the inputs lower than the cutoff and greater than minimum
    selection = (input_grid < settings.INPUT_MAX) & (input_grid >= settings.INPUT_MIN)
    particle_number[selection] += 1
    particle_weight[selection] += input_grid[selection]
    # Next we consider the cells where we have inputs greater than the cutoff
    selection = input_grid >= settings.INPUT_MAX
    particle_number[selection] += np.floor(input_grid[selection] / settings.INPUT_MAX)
    particle_weight[selection] += settings.INPUT_MAX
    # Finally, we have the remainder for when inputs are not multiples of the cutoff, and so we have a remainder
    selection = (input_grid - np.multiply(particle_number, particle_weight)) > settings.INPUT_MIN
    particle_remain_num[selection] += 1
    particle_remain[selection] += (input_grid[selection] - np.multiply(particle_number, particle_weight)[selection])
    return particle_number.astype('int'), particle_weight, particle_remain_num.astype('int'), particle_remain


def particle_grid_to_list(particle_number: np.array, particle_weight: np.array, particle_remain_num: np.array,
                          particle_remain: np.array, lon: np.array, lat: np.array):
    particle_lat, particle_lon, particle_mass = [], [], []
    for lat_index in range(particle_number.shape[0]):
        for lon_index in range(particle_number.shape[1]):
            if particle_number[lat_index, lon_index] > 0:
                for reps in range(particle_number[lat_index, lon_index]):
                    particle_lat.append(lat[lat_index])
                    particle_lon.append(lon[lon_index])
                    particle_mass.append(particle_weight[lat_index, lon_index])
            if particle_remain_num[lat_index, lon_index] > 0:
                particle_lat.append(lat[lat_index])
                particle_lon.append(lon[lon_index])
                particle_mass.append(particle_remain[lat_index, lon_index])
    particle_lat, particle_lon, particle_mass = np.array(particle_lat), np.array(particle_lon), np.array(particle_mass)
    non_zero = particle_mass > 0
    return particle_lat[non_zero], particle_lon[non_zero], particle_mass[non_zero]


def split_to_runs(particle_lat: np.array, particle_lon: np.array, particle_weight,
                  output_prefix: str, sub_division: int = settings.INPUT_DIV):
    run_number = len(particle_lat) // sub_division + 1
    if particle_weight is not None:
        var_dict = {'lon': particle_lon, 'lat': particle_lat, 'weight': particle_weight}
    else:
        var_dict = {'lon': particle_lon, 'lat': particle_lat}
    for run in range(run_number):
        for variables in var_dict:
            var_run = var_dict[variables][run * sub_division:(run + 1) * sub_division]
            np.save(output_prefix + '{}_run={}.npy'.format(variables, run), var_run)


def land_mark_checker(grid: np.array, lon: np.array, lat: np.array):
    mask = np.ma.getmask(grid)
    land = Field('Land', mask, lon=lon, lat=lat, transpose=False, mesh='spherical')
    return land


def number_of_releases(repeat_dt):
    if repeat_dt is None:
        releases = 1
    else:
        releases = math.floor(timedelta(days=365) / repeat_dt) + 1
    return releases
