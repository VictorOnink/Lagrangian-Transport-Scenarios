import settings as settings
import utils

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from scipy import io
import progressbar

def parcels_to_concentration(file_dict: dict):
    dataset = Dataset(settings.DATA_DIR_SERVERS[settings.SERVER] + 'HYCOM/HYCOM_Surface_3h_2000-01-01.nc')
    uo = dataset.variables['water_u'][0, 0, :, :]
    # Load HYCOM latitude, longitude and set a zero array with the same grid size as
    lat_hyc = dataset.variables['lat'][:]
    lon_hyc = dataset.variables['lon'][:]
    overall_concentration = np.zeros(uo.shape)
    # Generate the grid for the hexbin operation
    bin_number = 1000
    lon_min, lon_max = -180, 180
    lat_min, lat_max = -90, 90
    hex_grid = Hexagonal2DGrid((bin_number, bin_number), [lon_min, lon_max, lat_min, lat_max])

    output_direc = utils._get_output_directory(server=settings.SERVER) + 'concentrations/'
    # Counter for the number of files
    counter = 0
    # loop through the runs
    for run in progressbar.progressbar(range(settings.RUN_RANGE)):
        # Loop through the restart files
        for restart in range(settings.SIM_LENGTH):
            # Load the lon, lat, weight data
            parcels_file = file_dict[run][restart]
            dataset = Dataset(parcels_file)
            lon, lat = dataset.variables['lon'][:, :-1], dataset.variables['lat'][:, :-1]
            weight = dataset.variables['weights'][:, :-1] * settings.BUOYANT
            # Flatten the arrays
            lon_f, lat_f, weight_f = lon.flatten(), lat.flatten(), weight.flatten()
            # Carry out the hex bin operation
            hexagon_cumulative_sum, hexagon_coord = hexbin(lon_f, lat_f, weight_f, hex_grid)
            # Get the average per day, so divide by the number of days in the year
            hexagon_cumulative_sum /= lon.shape[1]
            # Get the concentration onto the HYCOM grid
            overall_concentration += utils.histogram(lon_data=hexagon_coord[:, 0], lat_data=hexagon_coord[:, 1],
                                                     bins_Lon=lon_hyc, bins_Lat=lat_hyc,
                                                     weight_data=hexagon_cumulative_sum
                                                     )
            counter += 1
    # Divide by the number of files we are averaging over
    overall_concentration /= counter
    # Saving the computed concentration
    prefix = 'concentration'
    output_name = output_direc + utils._analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
    io.savemat(output_name, {'concentration': overall_concentration})


class Hexagonal2DGrid(object):
    """
    Used to fix the gridsize, extent, and bins
    """

    def __init__(self, gridsize, extent, bins=None):
        self.gridsize = gridsize
        self.extent = extent
        self.bins = bins


def hexbin(x, y, c, hexgrid):
    """
    To hexagonally bin the data in 2 dimensions
    """
    # Basically you fix the gridsize, extent, and bins to keep them the same
    # then the resulting count array is the same
    hexbin = plt.hexbin(x, y, C=c, mincnt=1,
                        gridsize=hexgrid.gridsize,
                        extent=hexgrid.extent,
                        bins=hexgrid.bins,
                        reduce_C_function=np.sum)
    # you could close the figure if you don't want it
    plt.close('all')
    counts = hexbin.get_array().copy()
    center = hexbin.get_offsets()
    return counts, center
