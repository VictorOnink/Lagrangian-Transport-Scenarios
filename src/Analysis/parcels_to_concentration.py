import settings as settings
import utils
from advection_scenarios import advection_files
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar
from copy import deepcopy


class parcels_to_concentration():
    def __init__(self, file_dict: dict):
        self.parallel_step = settings.PARALLEL_STEP
        self.file_dict = file_dict
        self.LON, self.LAT, self.GRID, self.hexgrid = create_hex_grid()
        self.beach_label_dict = set_beach_label_dict(scenario_name=settings.SCENARIO_NAME)
        self.temp_direc, self.output_direc = get_directories(scenario_name=settings.SCENARIO_NAME)
        self.output_file_dict = create_output_file_dict(scenario_name=settings.SCENARIO_NAME, grid=self.GRID,
                                                        beach_states=self.beach_label_dict.keys(), lon=self.LON,
                                                        lat=self.LAT)

    def run(self):
        if self.parallel_step == 1:
            print_statement = 'year {}-{}, run {} restart {}'.format(settings.STARTYEAR, settings.STARTMONTH,
                                                                     settings.RUN, settings.RESTART)
            utils.print_statement(print_statement, to_print=True)
            # Loading the data
            parcels_dataset, post_dataset = load_parcels_post_output(scenario_name=settings.SCENARIO_NAME,
                                                                     file_dict=self.file_dict)
            # Loading, flattening and removing nan values for necessary data arrays
            full_data_dict = {}
            for variable in ['lon', 'lat', 'beach']:
                full_data_dict[variable] = parcels_dataset.variables[variable][:, :-1]
            full_data_dict, time_steps = complete_full_data_dict(scenario_name=settings.SCENARIO_NAME,
                                                                 full_data_dict=full_data_dict,
                                                                 parcels_dataset=parcels_dataset,
                                                                 post_dataset=post_dataset)
            # Looping through the beach states
            for beach_state in self.beach_label_dict.keys():
                state_data = {}
                for variable in ['lon', 'lat', 'weights', 'size_class']:
                    if variable in full_data_dict.keys():
                        beach_selection = full_data_dict['beach'] == self.beach_label_dict[beach_state]
                        state_data[variable] = full_data_dict[variable][beach_selection]
                    print(state_data.keys())
                    print(full_data_dict.keys())
                    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                        for size_class in range(settings.SIZE_CLASS_NUMBER):
                            size_class_data = {}
                            for variable in ['lon', 'lat', 'weights']:
                                size_selection = state_data['size_class'] == size_class
                                size_class_data[variable] = state_data[variable][size_selection]
                            key_year = utils.analysis_simulation_year_key(settings.RESTART)
                            self.output_dict[key_year][beach_state][size_class] = calculate_concentration(lon=size_class_data['lon'],
                                                                                                          lat=size_class_data['lat'],
                                                                                                          weights=size_class_data['weights'],
                                                                                                          hex_grid=self.hexgrid, time_steps=time_steps,
                                                                                                          lon_bin=self.LON, lat_bin=self.LAT)
                    else:
                        self.output_dict[key_year][beach_state] = calculate_concentration(lon=state_data['lon'],
                                                                                          lat=state_data['lat'],
                                                                                          weights=state_data['weights'],
                                                                                          hex_grid=self.hexgrid,
                                                                                          time_steps=time_steps,
                                                                                          lon_bin=self.LON,
                                                                                          lat_bin=self.LAT)
            output_name = get_file_names(scenario_name=settings.SCENARIO_NAME, file_dict=self.file_dict,
                                         directory=self.temp_direc, final=False)
            utils.save_obj(output_name, self.output_dict)
            str_format = settings.STARTYEAR, settings.STARTMONTH, settings.RUN, settings.RESTART
            print_statement = 'The concentration for year {}-{}, run {} restart {} has been save'.format(*str_format)
            utils.print_statement(print_statement, to_print=True)

        elif self.parallel_step == 2:
            pbar = ProgressBar()
            for ind_year, year in pbar(enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH))):
                for run in range(0, settings.RUN_RANGE):
                    for restart in range(0, settings.SIM_LENGTH - ind_year):
                        if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                            for month in range(1, 13):
                                file_name = get_file_names(scenario_name=settings.SCENARIO_NAME, file_dict=self.file_dict,
                                                           directory=self.temp_direc, final=False, year=year, month=month,
                                                           run=run, restart=restart)
                                dataset_post = utils.load_obj(filename=file_name)
                                for beach_state in self.beach_label_dict.keys():
                                    for size_class in range(settings.SIZE_CLASS_NUMBER):
                                        key_year = utils.analysis_simulation_year_key(restart)
                                        self.output_dict[key_year][beach_state][size_class] += dataset_post[key_year][beach_state][size_class]
                                utils.remove_file(file_name)

                        else:
                            file_name = get_file_names(scenario_name=settings.SCENARIO_NAME, file_dict=self.file_dict,
                                                       directory=self.temp_direc, final=False, year=year, run=run,
                                                       restart=restart)
                            dataset_post = utils.load_obj(filename=file_name)
                            for beach_state in self.beach_label_dict.keys():
                                key_year = utils.analysis_simulation_year_key(restart)
                                self.output_dict[key_year][beach_state] += dataset_post[key_year][beach_state]
                            utils.remove_file(file_name)
            # Calculating the average concentrations over the entire length of the simulation from the individual years
            for simulation_years in range(settings.SIM_LENGTH):
                key_year = utils.analysis_simulation_year_key(simulation_years)
                for beach_state in self.beach_label_dict.keys():
                    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                        for size_class in self.output_dict[key_year][beach_state].keys():
                            self.output_dict['overall_concentration'][beach_state][size_class] += self.output_dict[key_year][beach_state][size_class]
                    else:
                        self.output_dict['overall_concentration'][beach_state] += self.output_dict[key_year][beach_state]

            for beach_state in self.beach_label_dict.keys():
                if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                    for size_class in self.output_dict[key_year][beach_state].keys():
                        self.output_dict['overall_concentration'][beach_state][size_class] /= settings.SIM_LENGTH
                else:
                    self.output_dict['overall_concentration'][beach_state][size_class] /= settings.SIM_LENGTH

            # Saving the computed concentration
            output_name = get_file_names(scenario_name=settings.SCENARIO_NAME, file_dict=self.file_dict,
                                         directory=self.output_direc, final=True)
            utils.save_obj(output_name, self.output_dict)
            utils.print_statement("The concentration has been saved")

        else:
            ValueError('settings.PARALLEL_STEP can not have a value of {}'.format(self.parallel_step))


########################################################################################################################
"""
These following functions are general ones that are used by all scenarios to compute concentrations
"""


def create_output_file_dict(scenario_name, grid, beach_states, lon, lat):
    # Creating the base grid that has the dimensions of the output grid of the histogram
    lat_dim, lon_dim = grid.shape
    base_grid = np.zeros(shape=(lat_dim - 1, lon_dim - 1), dtype=float)
    beach_state_dict = dict.fromkeys(beach_states)
    if scenario_name in ['FragmentationKaandorpPartial']:
        size_dict = dict.fromkeys(range(settings.SIZE_CLASS_NUMBER))
        for size in size_dict.keys():
            size_dict[size] = deepcopy(base_grid)

    # Calculating the cell midpoints for lon and lat
    bin_mid_lon = 0.5 * lon[1:] + 0.5 * lon[:-1]
    bin_mid_lat = 0.5 * lat[1:] + 0.5 * lat[:-1]

    # Creating the beach state dictionary
    for state in beach_state_dict.keys():
        if scenario_name in ['FragmentationKaandorpPartial']:
            beach_state_dict[state] = deepcopy(size_dict)
        else:
            beach_state_dict[state] = deepcopy(base_grid)

    # Creating the final output dictionary
    output_dict = {'overall_concentration': deepcopy(beach_state_dict), 'lon': bin_mid_lon, 'lat': bin_mid_lat}
    for simulation_years in range(settings.SIM_LENGTH):
        output_dict[utils.analysis_simulation_year_key(simulation_years)] = deepcopy(beach_state_dict)

    return output_dict


def flatten_and_nan_removal(full_data_dict):
    # Flatten the arrays
    for variable in full_data_dict.keys():
        full_data_dict[variable] = full_data_dict[variable].flatten()
    # Remove all nan values
    is_nan = ~np.isnan(deepcopy(full_data_dict['lon']))
    for variable in full_data_dict.keys():
        full_data_dict[variable] = full_data_dict[variable][is_nan]
    return full_data_dict


def create_hex_grid():
    advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                        advection_scenario=settings.ADVECTION_DATA,
                                                        repeat_dt=None)
    adv_file_dict = advection_scenario.file_names
    LON, LAT, GRID = adv_file_dict['LON'], adv_file_dict['LAT'], adv_file_dict['GRID']
    bin_number = max(GRID.shape)
    lon_min, lon_max = np.nanmin(LON), np.nanmax(LON)
    lat_min, lat_max = np.nanmin(LAT), np.nanmax(LAT)
    hex_grid = Hexagonal2DGrid((bin_number, bin_number), [lon_min, lon_max, lat_min, lat_max])
    return LON, LAT, GRID, hex_grid


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


def get_directories(scenario_name):
    temp_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/temporary/'.format(scenario_name)
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format(scenario_name)
    utils.check_direc_exist(temp_direc)
    utils.check_direc_exist(output_direc)
    return temp_direc, output_direc


def get_file_names(scenario_name, file_dict, directory, final, year=settings.STARTYEAR, month=settings.STARTMONTH,
                   run=settings.RUN, restart=settings.RESTART):
    split = {True: None, False: '.nc'}[final]
    prefix = 'horizontal_concentration'
    if scenario_name in ['FragmentationKaandorpPartial']:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict['postprocess'][year][month][run][restart],
        prefix=prefix, split=split)
    else:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix, split=split)
    return output_name


def set_beach_label_dict(scenario_name):
    if scenario_name in ['FragmentationKaandorpPartial']:
        return {'beach': 1, 'adrift': 0}
    else:
        return {'beach': 1, 'adrift': 0, 'seabed': 3}


def load_parcels_post_output(scenario_name, file_dict, year=settings.STARTYEAR, month=settings.STARTMONTH,
                             run=settings.RUN, restart=settings.RESTART):
    if scenario_name in ['FragmentationKaandorpPartial']:
        parcels_dataset = Dataset(file_dict['parcels'][year][month][run][restart])
        post_dataset = utils.load_obj(file_dict['postprocess'][year][month][run][restart])
    else:
        parcels_dataset = Dataset(file_dict[run][restart])
        post_dataset = None
    return parcels_dataset, post_dataset


def complete_full_data_dict(scenario_name, full_data_dict, parcels_dataset, post_dataset):
    if scenario_name in ['FragmentationKaandorpPartial']:
        full_data_dict['weights'] = post_dataset['particle_number'][:, :-1]
        full_data_dict['size_class'] = parcels_dataset.variables['size_class'][:, :-1]
    else:
        if 'weights' in parcels_dataset.variables.keys():
            full_data_dict['weights'] = parcels_dataset.variables['weights'][:, :-1] * settings.BUOYANT
        else:
            full_data_dict['weights'] = np.ones(full_data_dict['lon'].shape, dtype=np.float32)
    time_steps = full_data_dict['lon'].shape[1]
    # Flatten the arrays and remove nana values
    full_data_dict = flatten_and_nan_removal(full_data_dict)
    return full_data_dict, time_steps


def calculate_concentration(lon, lat, weights, hex_grid, time_steps, lon_bin, lat_bin):
    hexagon_cumulative_sum, hexagon_coord = hexbin(lon, lat,
                                                   weights, hex_grid)
    # Get the average per day, so divide by the number of days in the year
    hexagon_cumulative_sum /= time_steps
    # Get the concentration onto the advection grid
    key_year = utils.analysis_simulation_year_key(settings.RESTART)
    concentration, _, _ = utils.histogram(lon_data=hexagon_coord[:, 0],
                                                      lat_data=hexagon_coord[:, 1],
                                                      bins_Lon=lon_bin, bins_Lat=lat_bin,
                                                      weight_data=hexagon_cumulative_sum
                                                      )
    return concentration