import settings
import utils

import numpy as np
import glob
import os
import xarray
from netCDF4 import Dataset
import pandas as pd
import fiona
from collections import Iterable
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cpf
from copy import deepcopy
from geopy import distance
from progressbar import ProgressBar
import math

def create_input_files(prefix: str, grid: np.array, lon: np.array, lat: np.array):
    # Create the output prefix and then check if any files with such prefixes exist
    if settings.INPUT in ['Jambeck', 'Lebreton']:
        output_prefix = settings.INPUT_DIREC + settings.INPUT + '_' + prefix+'_2010_'
    else:
        os.system('echo "Perhaps take another look at what input you are using?"')

    #Check existence
    if len(glob.glob(output_prefix+'*') > 0):
        os.system('echo "The input files {} are already present"'.format(output_prefix))
    else:
        os.system('echo "We need to create the input files {}"'.format(output_prefix))
        if settings.INPUT == 'Jambeck':
            # Get the population data
            dataset = Dataset(settings.INPUT_DIREC + 'gpw_v4_population_count_adjusted_rev11_2pt5_min.nc')
            population = np.array(dataset.variables['UN WPP-Adjusted Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][2, :, :])
            lon_population = dataset.variables['longitude'][:]
            lat_population = dataset.variables['latitude'][:]
            # Get the mismanaged plastic fraction
            mismanaged = get_mismanaged_fraction_Jambeck(grid=grid, lon=lon, lat=lat, prefix=prefix)
            # Get the distance from land to shore
        elif settings.INPUT == 'Lebreton':
            lebData = pd.read_csv(settings.INPUT_DIREC + 'PlasticRiverInputs.csv')
            lon_inputs, lat_inputs = np.array(lebData['X']), np.array(lebData['Y'])
            plastic_inputs = np.array(lebData['i_low'])
        # Only keep the particles that are within the domain
        lon_inputs, lat_inputs, plastic_inputs = within_domain(lon=lon, lat=lat, lon_inputs=lon_inputs,
                                                               lat_inputs=lat_inputs, plastic_inputs=plastic_inputs)
        # Get the inputs onto the grid of the advection data
        inputs_grid = utils.histogram(lon_data=lon_inputs, lat_data=lat_inputs, bins_Lon=lon, bins_Lat=lat,
                                      weight_data=plastic_inputs, area_correc=False)
        # Get the ocean cells adjacent to coastal ocean cells, as these are the ones in which the particles will
        # be placed. For brevity, we will refer to these cells as coastal in the code here
        coastal = get_coastal_extended(grid=grid, coastal=get_coastal_cells(grid=grid))
        os.system('echo "Working up to getting the coastal fraction"')

def get_mismanaged_fraction_Jambeck(grid: np.array, lon: np.array, lat: np.array, prefix: str):
    mismanaged_file = settings.INPUT_DIREC + 'Jambeck_mismanaged_fraction_2010.nc'
    if utils._check_file_exist(mismanaged_file):
        dataset = Dataset(mismanaged_file)
        return dataset.variables['mismanaged fraction'][:]
    else:
        os.system('echo "I can not find the original code that was used to generate this so I will fill this in '
                  'later..."')


def within_domain(lon: np.array, lat: np.array, lon_inputs: np.array, lat_inputs: np.array, plastic_inputs: np.array):
    lon_max, lon_min = np.max(lon), np.min(lon)
    lat_max, lat_min = np.max(lat), np.min(lon)
    # Check which cells are within the domain
    domain = (lon_inputs <= lon_max) & (lon_inputs >= lon_min) & (lat_inputs >= lat_min) & (lat_inputs <= lat_max)
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