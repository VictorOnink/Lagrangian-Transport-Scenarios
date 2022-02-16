import settings as settings
import utils
from advection_scenarios import advection_files
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar
from copy import deepcopy
from datetime import datetime, timedelta


class parcels_to_concentration_monthly:
    def __init__(self, file_dict: dict):
        self.parallel_step = settings.PARALLEL_STEP
        self.file_dict = file_dict
        self.LON, self.LAT, self.GRID, self.MIN_DEPTH = self.get_lonlat_grid_depth()
        self.beach_label_dict = self.set_beach_label_dict()
        self.temp_direc, self.output_direc = self.get_directories()
        self.weight_list = ['particle_mass_sink', 'particle_number_sink']
        self.depth_level = ['surface_1m', 'surface_5m', 'column']
        self.data_variable_list = self.set_data_variable_list()
        self.output_dict = self.create_output_file_dict()
        self.parcels_dataset, self.post_dataset = self.load_parcels_post_output()
        self.month_indices = self.get_month_indices()

    def run(self):
        if self.parallel_step == 1:
            print_statement = 'year {}-{}, run {} restart {}'.format(settings.STARTYEAR, settings.STARTMONTH,
                                                                     settings.RUN, settings.RESTART)
            utils.print_statement(print_statement, to_print=True)
            output_name = self.get_file_names(directory=self.temp_direc, final=False)
            if not utils.check_file_exist(output_name, without_pkl=True):
                for month in range(1, 13):
                    # Getting the indices for the month
                    start_ind, end_ind = self.month_indices[month][0], self.month_indices[month][1]
                    # Loading, flattening and removing nan values for necessary data arrays
                    full_data_dict = {}
                    for variable in self.data_variable_list:
                        full_data_dict[variable] = self.parcels_dataset.variables[variable][:, start_ind:end_ind]
                    full_data_dict, time_steps = self.complete_full_data_dict(full_data_dict=full_data_dict,
                                                                              parcels_dataset=self.parcels_dataset,
                                                                              post_dataset=self.post_dataset,
                                                                              start_ind=start_ind, end_ind=end_ind)
                    # Looping through the beach states
                    for beach_state in self.beach_label_dict.keys():
                        state_data = {}
                        for variable in utils.flatten_list_of_lists([['lon', 'lat', 'weights', 'size_class', 'z'], self.weight_list]):
                            if variable in full_data_dict.keys():
                                beach_selection = full_data_dict['beach'] == self.beach_label_dict[beach_state]
                                state_data[variable] = full_data_dict[variable][beach_selection]
                        if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                            for size_class in range(settings.SIZE_CLASS_NUMBER):
                                size_class_data = {}
                                for variable in utils.flatten_list_of_lists([['lon', 'lat'], self.weight_list]):
                                    size_selection = state_data['size_class'] == size_class
                                    size_class_data[variable] = state_data[variable][size_selection]
                                key_year = utils.analysis_simulation_year_key(settings.RESTART)
                                for weight in self.weight_list:
                                    # First, the total water column
                                    output = self.calculate_concentration(lon=size_class_data['lon'],
                                                                          lat=size_class_data['lat'],
                                                                          weights=size_class_data[weight],
                                                                          time_steps=time_steps)
                                    self.output_dict[key_year][month][beach_state][weight][size_class]['column'] = output
                                    # Next, the concentration within 5m of the ocean surface (defined at self.MIN_DEPTH)
                                    selec = size_class_data['z'] < (self.MIN_DEPTH + 5)
                                    if np.sum(selec) > 0:
                                        output = self.calculate_concentration(lon=size_class_data['lon'][selec],
                                                                              lat=size_class_data['lat'][selec],
                                                                              weights=size_class_data[weight][selec],
                                                                              time_steps=time_steps)
                                        self.output_dict[key_year][month][beach_state][weight][size_class]['surface_5m'] = output
                                    # Finally, the concentration within 1m of the ocean surface
                                    selec = size_class_data['z'] < (self.MIN_DEPTH + 1)
                                    if np.sum(selec) > 0:
                                        output = self.calculate_concentration(lon=size_class_data['lon'][selec],
                                                                              lat=size_class_data['lat'][selec],
                                                                              weights=size_class_data[weight][selec],
                                                                              time_steps=time_steps)
                                        self.output_dict[key_year][month][beach_state][weight][size_class]['surface_1m'] = output
                        elif settings.SCENARIO_NAME in ['SizeTransport']:
                            key_year = utils.analysis_simulation_year_key(settings.RESTART)
                            # First, the total water column
                            output = self.calculate_concentration(lon=state_data['lon'], lat=state_data['lat'],
                                                                  weights=state_data['weights'], time_steps=time_steps)
                            self.output_dict[key_year][month][beach_state]['column'] = output
                            # Next, the concentration within 5m of the ocean surface (defined at self.MIN_DEPTH)
                            selec = (state_data['z'] < 11) & (state_data['z'] > 10)
                            if np.sum(selec) > 0:
                                output = self.calculate_concentration(lon=state_data['lon'][selec],
                                                                      lat=state_data['lat'][selec],
                                                                      weights=state_data['weights'],
                                                                      time_steps=time_steps)
                                self.output_dict[key_year][month][beach_state]['surface_5m'] = output
                            # Finally, the concentration within 1m of the ocean surface
                            selec = state_data['z'] < (self.MIN_DEPTH + 1)
                            if np.sum(selec) > 0:
                                output = self.calculate_concentration(lon=state_data['lon'][selec],
                                                                      lat=state_data['lat'][selec],
                                                                      weights=state_data['weights'],
                                                                      time_steps=time_steps)
                                self.output_dict[key_year][month][beach_state]['surface_1m'] = output

                        else:
                            key_year = utils.analysis_simulation_year_key(settings.RESTART)
                            self.output_dict[key_year][month][beach_state] = self.calculate_concentration(lon=state_data['lon'],
                                                                                                   lat=state_data['lat'],
                                                                                                   weights=state_data['weights'],
                                                                                                   time_steps=time_steps)
                utils.save_obj(output_name, self.output_dict)
                str_format = settings.STARTYEAR, settings.STARTMONTH, settings.RUN, settings.RESTART
                print_statement = 'The concentration for year {}-{}, run {} restart {} has been save'.format(*str_format)
                utils.print_statement(print_statement, to_print=True)

        elif self.parallel_step == 2:
            pbar = ProgressBar()
            for ind_year, year in pbar(enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH))):
                for month_ind in range(1, 13):
                    for run in range(0, settings.RUN_RANGE):
                        for restart in range(0, settings.SIM_LENGTH - ind_year):
                            file_name = self.get_file_names(directory=self.temp_direc, final=False, year=year,
                                                            month=month_ind, run=run, restart=restart)
                            if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                                dataset_post = utils.load_obj(filename=file_name)
                                for month in range(1, 13):
                                    for beach_state in self.beach_label_dict.keys():
                                        for weight in self.weight_list:
                                            for size_class in range(settings.SIZE_CLASS_NUMBER):
                                                for depth in self.depth_level:
                                                    key_year = utils.analysis_simulation_year_key(restart + ind_year)
                                                    self.output_dict[key_year][month][beach_state][weight][size_class][depth] += dataset_post[key_year][month][beach_state][weight][size_class][depth]
                                utils.remove_file(file_name + '.pkl')
                            else:
                                if month_ind == 1 and year == settings.STARTYEAR:
                                    dataset_post = utils.load_obj(filename=file_name)
                                    for month in range(1, 13):
                                        for beach_state in self.beach_label_dict.keys():
                                            key_year = utils.analysis_simulation_year_key(restart + ind_year)
                                            if settings.SCENARIO_NAME in ['SizeTransport']:
                                                for depth in self.depth_level:
                                                    self.output_dict[key_year][month][beach_state][depth] += dataset_post[key_year][month][beach_state][depth]
                                            else:
                                                self.output_dict[key_year][month][beach_state] += dataset_post[key_year][month][beach_state]
                                    utils.remove_file(file_name + '.pkl')
            # Saving the computed concentration
            output_name = self.get_file_names(directory=self.output_direc, final=True)
            utils.save_obj(output_name, self.output_dict)
            utils.print_statement("The concentration has been saved")

        else:
            ValueError('settings.PARALLEL_STEP can not have a value of {}'.format(self.parallel_step))

    def create_output_file_dict(self):
        # Creating the base grid that has the dimensions of the output grid of the histogram
        lat_dim, lon_dim = self.GRID.shape
        base_grid = np.zeros(shape=(lat_dim - 1, lon_dim - 1), dtype=float)
        beach_state_dict = dict.fromkeys(self.beach_label_dict.keys())
        if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
            size_dict = dict.fromkeys(range(settings.SIZE_CLASS_NUMBER))
            for size in size_dict.keys():
                size_dict[size] = deepcopy(base_grid)

        # Calculating the cell midpoints for lon and lat
        bin_mid_lon = 0.5 * self.LON[1:] + 0.5 * self.LON[:-1]
        bin_mid_lat = 0.5 * self.LAT[1:] + 0.5 * self.LAT[:-1]

        # Creating the beach state dictionary
        for state in beach_state_dict.keys():
            if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                beach_state_dict[state] = {}
                for weight in self.weight_list:
                    beach_state_dict[state][weight] = {}
                    for depth in self.depth_level:
                        beach_state_dict[state][weight][depth] = deepcopy(size_dict)
            elif settings.SCENARIO_NAME in ['SizeTransport']:
                beach_state_dict[state] = {}
                for depth in self.depth_level:
                    beach_state_dict[state][depth] = deepcopy(base_grid)
            else:
                beach_state_dict[state] = deepcopy(base_grid)

        # Creating the final output dictionary
        output_dict = {'lon': bin_mid_lon, 'lat': bin_mid_lat}
        for simulation_years in range(settings.SIM_LENGTH):
            year_key = utils.analysis_simulation_year_key(simulation_years)
            output_dict[year_key] = {}
            for month in range(1, 13):
                output_dict[year_key][month] = deepcopy(beach_state_dict)

        return output_dict

    def complete_full_data_dict(self, full_data_dict, parcels_dataset, post_dataset, start_ind, end_ind):
        if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
            for weight in self.weight_list:
                full_data_dict[weight] = post_dataset[weight][:, start_ind:end_ind]
            full_data_dict['size_class'] = parcels_dataset.variables['size_class'][:, start_ind:end_ind]
        else:
            if 'weights' in parcels_dataset.variables.keys():
                full_data_dict['weights'] = parcels_dataset.variables['weights'][:, start_ind:end_ind] * settings.BUOYANT
            else:
                full_data_dict['weights'] = np.ones(full_data_dict['lon'].shape, dtype=np.float32)
        time_steps = full_data_dict['lon'].shape[1]
        # Flatten the arrays and remove nana values
        full_data_dict = self.flatten_and_nan_removal(full_data_dict)
        return full_data_dict, time_steps

    def get_file_names(self, directory, final, year=settings.STARTYEAR, month=settings.STARTMONTH,
                       run=settings.RUN, restart=settings.RESTART):
        split = {True: None, False: '.nc'}[final]
        prefix = 'horizontal_concentration'
        if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
            output_name = directory + utils.analysis_save_file_name(
                input_file=self.file_dict['postprocess'][year][month][run][restart],
                prefix=prefix, split=split)
        else:
            output_name = directory + utils.analysis_save_file_name(input_file=self.file_dict[run][restart],
                                                                    prefix=prefix, split=split)
        return output_name

    def load_parcels_post_output(self, year=settings.STARTYEAR, month=settings.STARTMONTH, run=settings.RUN,
                                 restart=settings.RESTART):
        if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
            parcels_dataset = Dataset(self.file_dict['parcels'][year][month][run][restart])
            post_dataset = utils.load_obj(self.file_dict['postprocess'][year][month][run][restart])
        else:
            parcels_dataset = Dataset(self.file_dict[run][restart])
            post_dataset = None
        return parcels_dataset, post_dataset

    def calculate_concentration(self, lon, lat, time_steps, weights):
        # Get the concentration onto the advection grid
        concentration, _, _ = utils.histogram(lon_data=lon, lat_data=lat, bins_Lon=self.LON, bins_Lat=self.LAT,
                                              weight_data=weights)
        # Divide the concentration by the number of timesteps within the concentration
        concentration /= time_steps
        return concentration

    @staticmethod
    def get_month_indices():
        # Setting the start time of the simulation, and defining the output dt and time index
        current_time = datetime(settings.STARTYEAR + settings.RESTART, settings.STARTMONTH, settings.STARTDAY)
        end_time = datetime(settings.STARTYEAR + settings.RESTART + 1, 1, 1)
        output_dt = settings.OUTPUT_TIME_STEP
        index_dt = 0
        # Creating the dictionary containing the start and end indices for each month
        month_indices = {}
        while current_time < end_time:
            if current_time.month not in month_indices.keys():
                month_indices[current_time.month] = [index_dt, index_dt]
            else:
                month_indices[current_time.month][1] = index_dt
            current_time += output_dt
            index_dt += 1
        return month_indices

    @staticmethod
    def get_lonlat_grid_depth():
        advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                            advection_scenario=settings.ADVECTION_DATA,
                                                            repeat_dt=None)
        adv_file_dict = advection_scenario.file_names
        LON, LAT, GRID = adv_file_dict['LON'], adv_file_dict['LAT'], adv_file_dict['GRID']
        MIN_DEPTH = np.nanmin(adv_file_dict['DEPTH'])
        return LON, LAT, GRID, MIN_DEPTH

    @staticmethod
    def flatten_and_nan_removal(full_data_dict):
        # Flatten the arrays
        for variable in full_data_dict.keys():
            full_data_dict[variable] = full_data_dict[variable].flatten()
        # Remove all nan values
        is_nan = ~np.isnan(deepcopy(full_data_dict['lon']))
        for variable in full_data_dict.keys():
            full_data_dict[variable] = full_data_dict[variable][is_nan]
        return full_data_dict

    @staticmethod
    def get_directories():
        temp_direc = settings.SCRATCH_DIR
        output_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format(settings.SCENARIO_NAME)
        utils.check_direc_exist(temp_direc)
        utils.check_direc_exist(output_direc)
        return temp_direc, output_direc

    @staticmethod
    def set_beach_label_dict():
        if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
            return {'beach': 1, 'adrift': 0}
        else:
            return {'beach': 1, 'adrift': 0, 'seabed': 3}

    @staticmethod
    def set_data_variable_list():
        if settings.SCENARIO_NAME in ["SizeTransport", 'FragmentationKaandorpPartial']:
            return ['lon', 'lat', 'beach', 'z']
        else:
            return ['lon', 'lat', 'beach']
