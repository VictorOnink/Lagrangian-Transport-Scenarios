import settings as settings
import utils
from advection_scenarios import advection_files
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import progressbar
import os
from copy import deepcopy


def parcels_to_concentration(file_dict: dict):
    """
    We calculate the horizontal concentrations on the grid of the circulation data. The concentration is calculated
    either as an average over the entire simulation, or over each given year. We also distinguish between whether
    particles are beached, afloat, or stuck at the sea floor
    :param file_dict:
    :return:
    """
    # Generate the grid for the hexbin operation
    advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                        advection_scenario=settings.ADVECTION_DATA,
                                                        repeat_dt=None)
    adv_file_dict = advection_scenario.file_names
    bin_number = 5000
    LON, LAT, GRID = adv_file_dict['LON'], adv_file_dict['LAT'], adv_file_dict['GRID']
    lon_min, lon_max = np.nanmin(LON), np.nanmax(LON)
    lat_min, lat_max = np.nanmin(LAT), np.nanmax(LAT)
    hex_grid = Hexagonal2DGrid((bin_number, bin_number), [lon_min, lon_max, lat_min, lat_max])

    # Getting the directory saving the output files
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format(settings.SCENARIO_NAME)
    utils.check_direc_exist(output_direc)

    # Counter for the number of files
    counter = 0

    # Create the output dictionarynp.zeros(GRID.shape)
    beach_state_dict = {'beach': np.zeros(GRID.shape), 'afloat': np.zeros(GRID.shape), 'seabed': np.zeros(GRID.shape)}
    beach_label_dict = {'beach': 1, 'afloat': 0, 'seabed': 3}
    output_dict = {'overall_concentration': deepcopy(beach_state_dict), 'lon': LON, 'lat': LAT}
    for simulation_years in range(settings.SIM_LENGTH):
        output_dict[utils.analysis_simulation_year_key(simulation_years)] = deepcopy(beach_state_dict)

    # loop through the runs
    for run in progressbar.progressbar(range(settings.RUN_RANGE)):
        # Loop through the restart files
        for restart in range(settings.SIM_LENGTH):
            # Load the lon, lat, weight data
            parcels_file = file_dict[run][restart]
            dataset = Dataset(parcels_file)
            full_data_dict = {}
            for variable in ['lon', 'lat', 'beach']:
                full_data_dict[variable] = dataset.variables[variable][:, :-1]
            if 'weights' in dataset.variables.keys():
                full_data_dict['weights'] = dataset.variables['weights'][:, :-1] * settings.BUOYANT
            else:
                full_data_dict['weights'] = np.ones(full_data_dict['lon'].shape, dtype=np.float32)
            # Saving the number of time steps saved in the data file
            time_steps = full_data_dict['lon'].shape[1]
            # Flatten the arrays
            for variable in full_data_dict.keys():
                full_data_dict[variable] = full_data_dict[variable].flatten()
            # Now, we sort through the various particle states
            for beach_state in beach_state_dict.keys():
                state_data = {}
                for variable in ['lon', 'lat', 'weights']:
                    state_data[variable] = full_data_dict[variable][full_data_dict['beach'] == beach_label_dict[beach_state]]
                # Carry out the hex bin operation
                hexagon_cumulative_sum, hexagon_coord = hexbin(state_data['lon'], state_data['lat'],
                                                               state_data['weights'], hex_grid)
                # Get the average per day, so divide by the number of days in the year
                hexagon_cumulative_sum /= time_steps
                # Get the concentration onto the advection grid
                key_year = utils.analysis_simulation_year_key(restart)
                output_dict[key_year][beach_state] += utils.histogram(lon_data=hexagon_coord[:, 0],
                                                                      lat_data=hexagon_coord[:, 1],
                                                                      bins_Lon=LON, bins_Lat=LAT,
                                                                      weight_data=hexagon_cumulative_sum
                                                                      )

    # Dividing the end of year concentrations by the number of runs
    for simulation_years in range(settings.SIM_LENGTH):
        key_year = utils.analysis_simulation_year_key(restart)
        for beach_state in output_dict[key_year].keys():
            output_dict[key_year][beach_state] /= settings.RUN_RANGE

    # Calculating the average concentrations over the entire length of the simulation from the individual years
    for beach_state in beach_state_dict.keys():
        for simulation_years in range(settings.SIM_LENGTH):
            key_year = utils.analysis_simulation_year_key(restart)
            output_dict['overall_concentration'][beach_state] += output_dict[key_year][beach_state]
    for beach_state in beach_state_dict.keys():
        output_dict['overall_concentration'][beach_state] /= settings.SIM_LENGTH

    # Saving the computed concentration
    prefix = 'horizontal_concentration'
    output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
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


