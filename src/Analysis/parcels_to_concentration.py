import settings as settings
import utils
from advection_scenarios import advection_files
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from scipy import io
import progressbar
import os


def parcels_to_concentration(file_dict: dict):
    # Generate the grid for the hexbin operation
    advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                        advection_scenario=settings.ADVECTION_DATA,
                                                        repeat_dt=None)
    adv_file_dict = advection_scenario.file_names
    bin_number = 1000
    LON, LAT, GRID = adv_file_dict['LON'], adv_file_dict['LAT'], adv_file_dict['GRID']
    lon_min, lon_max = np.nanmin(LON), np.nanmax(LON)
    lat_min, lat_max = np.nanmin(LAT), np.nanmax(LAT)
    hex_grid = Hexagonal2DGrid((bin_number, bin_number), [lon_min, lon_max, lat_min, lat_max])
    # Getting the grid for the final concentrations

    output_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format(settings.SCENARIO_NAME)
    utils.check_direc_exist(output_direc)
    # Counter for the number of files
    counter = 0
    # Create the output dictionary
    output_dict = {'overall_concentration': np.zeros(GRID.shape), 'lon': LON, 'lat': LAT}
    for simulation_years in range(settings.SIM_LENGTH):
        output_dict[simulation_years] = np.zeros(GRID.shape)

    # loop through the runs
    for run in progressbar.progressbar(range(settings.RUN_RANGE)):
        # Loop through the restart files
        for restart in range(settings.SIM_LENGTH):
            # Load the lon, lat, weight data
            parcels_file = file_dict[run][restart]
            dataset = Dataset(parcels_file)
            lon, lat = dataset.variables['lon'][:, :-1], dataset.variables['lat'][:, :-1]
            if 'weights' in dataset.variables.keys():
                weight = dataset.variables['weights'][:, :-1] * settings.BUOYANT
            else:
                weight = np.ones(lon.shape,dtype=np.float32)
            # Flatten the arrays
            lon_f, lat_f, weight_f = lon.flatten(), lat.flatten(), weight.flatten()
            # Carry out the hex bin operation
            hexagon_cumulative_sum, hexagon_coord = hexbin(lon_f, lat_f, weight_f, hex_grid)
            # Get the average per day, so divide by the number of days in the year
            hexagon_cumulative_sum /= lon.shape[1]
            # Get the concentration onto the advection grid
            output_dict['overall_concentration'] += utils.histogram(lon_data=hexagon_coord[:, 0],
                                                                    lat_data=hexagon_coord[:, 1],
                                                                    bins_Lon=LON, bins_Lat=LAT,
                                                                    weight_data=hexagon_cumulative_sum
                                                                    )
            counter += 1
            # Getting the concentration for the last particle positions with a simulation year
            hexagon_last, hexagon_coord = hexbin(lon[:, -1], lat[:, -1], weight[:, -1], hex_grid)
            output_dict[restart] += utils.histogram(lon_data=hexagon_coord[:, 0], lat_data=hexagon_coord[:, 1],
                                                    bins_Lon=LON, bins_Lat=LAT, weight_data=hexagon_last)

    # Divide by the number of files we are averaging over
    output_dict['overall_concentration'] /= counter

    # Dividing the end of year concentrations by the number of runs
    for restart in range(settings.SIM_LENGTH):
        output_dict[restart] /= settings.RUN_RANGE

    # Saving the computed concentration
    prefix = 'concentration'
    output_name = output_direc + utils._analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
    utils.save_obj(output_name, output_dict)
    os.system('echo "The concentration has been saved"')


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

