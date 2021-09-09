import settings as settings
import utils
from advection_scenarios import advection_files
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar
from copy import deepcopy


if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
    def parcels_to_concentration(file_dict: dict):
        """
        We calculate the horizontal concentrations on the grid of the circulation data. The concentration is calculated
        either as an average over the entire simulation, or over each given year. We also distinguish between whether
        particles are beached, adrift, or stuck at the sea floor
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

        # Create the output dictionary
        lat_dim, lon_dim = GRID.shape
        base_grid = np.zeros(shape=(lat_dim - 1, lon_dim - 1), dtype=float)
        size_dict = dict.fromkeys(range(settings.SIZE_CLASS_NUMBER))
        for size in size_dict.keys():
            size_dict[size] = deepcopy(base_grid)
        beach_state_dict = {'beach': deepcopy(size_dict), 'adrift': deepcopy(size_dict)}
        beach_label_dict = {'beach': 1, 'adrift': 0}
        output_dict = {'overall_concentration': deepcopy(beach_state_dict)}
        for simulation_years in range(settings.SIM_LENGTH):
            output_dict[utils.analysis_simulation_year_key(simulation_years)] = deepcopy(beach_state_dict)

        # loop through the runs
        pbar = ProgressBar()
        for ind_year, year in pbar(enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH))):
            for month in range(1, 13):
                for run in range(0, settings.RUN_RANGE):
                    # Loop through the restart files
                    for restart in range(0, settings.SIM_LENGTH - ind_year):
                        print_statement = 'year {}-{}, run {} restart {}'.format(year, month, run, restart)
                        utils.print_statement(print_statement, to_print=True)
                        # Load the lon, lat, weight data
                        parcels_file = file_dict['parcels'][year][month][run][restart]
                        post_file = file_dict['postprocess'][year][month][run][restart]
                        dataset_parcels = Dataset(parcels_file)
                        dataset_post = utils.load_obj(post_file)
                        full_data_dict = {}
                        for variable in ['lon', 'lat', 'beach', 'size_class']:
                            full_data_dict[variable] = dataset_parcels.variables[variable][:, :-1]
                        # Adding the particle number data, which we use to weigh the particle concentration
                        full_data_dict['weights'] = dataset_post['particle_number'][:, :-1]
                        # Saving the number of time steps saved in the data file
                        time_steps = full_data_dict['lon'].shape[1]
                        # Flatten the arrays
                        for variable in full_data_dict.keys():
                            full_data_dict[variable] = full_data_dict[variable].flatten()
                        # Now, we sort through the various particle states
                        for beach_state in beach_state_dict.keys():
                            state_data = {}
                            for variable in ['lon', 'lat', 'weights', 'size_class']:
                                state_data[variable] = full_data_dict[variable][full_data_dict['beach'] == beach_label_dict[beach_state]]
                            for size_class in range(settings.SIZE_CLASS_NUMBER):
                                size_class_data = {}
                                for variable in ['lon', 'lat', 'weights', ]:
                                    size_class_data[variable] = state_data[variable][state_data['size_class'] == size_class]
                                    # Carry out the hex bin operation
                                    hexagon_cumulative_sum, hexagon_coord = hexbin(state_data['lon'], state_data['lat'],
                                                                                   state_data['weights'], hex_grid)
                                    # Get the average per day, so divide by the number of days in the year
                                    hexagon_cumulative_sum /= time_steps
                                    # Get the concentration onto the advection grid
                                    key_year = utils.analysis_simulation_year_key(restart)
                                    concentration, lat_mid, lon_mid = utils.histogram(lon_data=hexagon_coord[:, 0],
                                                                                      lat_data=hexagon_coord[:, 1],
                                                                                      bins_Lon=LON, bins_Lat=LAT,
                                                                                      weight_data=hexagon_cumulative_sum
                                                                                      )
                                    output_dict[key_year][beach_state][size_class] += concentration

        # Adding the lon/lat arrays
        output_dict['lon'], output_dict['lat'] = lon_mid, lat_mid

        # Dividing the end of year concentrations by the number of runs
        for simulation_years in range(settings.SIM_LENGTH):
            key_year = utils.analysis_simulation_year_key(simulation_years)
            for beach_state in output_dict[key_year].keys():
                for size_class in output_dict[key_year][beach_state].keys():
                    output_dict[key_year][beach_state][size_class] /= settings.RUN_RANGE

        # Calculating the average concentrations over the entire length of the simulation from the individual years
        for simulation_years in range(settings.SIM_LENGTH):
            key_year = utils.analysis_simulation_year_key(simulation_years)
            for beach_state in beach_state_dict.keys():
                for size_class in output_dict[key_year][beach_state].keys():
                    output_dict['overall_concentration'][beach_state][size_class] += output_dict[key_year][beach_state][size_class]

        for beach_state in beach_state_dict.keys():
            for size_class in output_dict[key_year][beach_state].keys():
                output_dict['overall_concentration'][beach_state][size_class] /= settings.SIM_LENGTH

        # Saving the computed concentration
        prefix = 'horizontal_concentration'
        output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict['postprocess'][settings.STARTYEAR][1][0][restart], prefix=prefix)
        utils.save_obj(output_name, output_dict)
        utils.print_statement("The concentration has been saved")

else:
    def parcels_to_concentration(file_dict: dict):
        """
        We calculate the horizontal concentrations on the grid of the circulation data. The concentration is calculated
        either as an average over the entire simulation, or over each given year. We also distinguish between whether
        particles are beached, adrift, or stuck at the sea floor
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
        output_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format(
            settings.SCENARIO_NAME)
        utils.check_direc_exist(output_direc)

        # Create the output dictionarynp.zeros(GRID.shape)
        beach_state_dict = {'beach': np.zeros(GRID.shape), 'adrift': np.zeros(GRID.shape),
                            'seabed': np.zeros(GRID.shape)}
        beach_label_dict = {'beach': 1, 'adrift': 0, 'seabed': 3}
        output_dict = {'overall_concentration': deepcopy(beach_state_dict)}
        for simulation_years in range(settings.SIM_LENGTH):
            output_dict[utils.analysis_simulation_year_key(simulation_years)] = deepcopy(beach_state_dict)

        # loop through the runs
        pbar = ProgressBar()
        for run in pbar(range(settings.RUN_RANGE)):
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
                        state_data[variable] = full_data_dict[variable][
                            full_data_dict['beach'] == beach_label_dict[beach_state]]
                    # Carry out the hex bin operation
                    hexagon_cumulative_sum, hexagon_coord = hexbin(state_data['lon'], state_data['lat'],
                                                                   state_data['weights'], hex_grid)
                    # Get the average per day, so divide by the number of days in the year
                    hexagon_cumulative_sum /= time_steps
                    # Get the concentration onto the advection grid
                    key_year = utils.analysis_simulation_year_key(restart)
                    concentration, lat_mid, lon_mid = utils.histogram(lon_data=hexagon_coord[:, 0],
                                                                      lat_data=hexagon_coord[:, 1],
                                                                      bins_Lon=LON, bins_Lat=LAT,
                                                                      weight_data=hexagon_cumulative_sum
                                                                      )

                    output_dict[key_year][beach_state] += concentration

        # Adding the lon/lat arrays
        output_dict['lon'], output_dict['lat'] = lon_mid, lat_mid

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
        utils.print_statement("The concentration has been saved")


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


