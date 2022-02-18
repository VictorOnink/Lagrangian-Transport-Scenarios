import settings as settings
import utils
from advection_scenarios import advection_files
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar
from copy import deepcopy
from datetime import datetime


class parcels_to_spatial_vertical_profiles:
    def __init__(self, file_dict: dict):
        assert settings.SCENARIO_NAME in ['SizeTransport'], 'The spatial vertical profiles is currently only set ' \
                                                            'up for SizeTransport'
        self.file_dict = file_dict
        self.parallel_step = settings.PARALLEL_STEP
        self.file_dict = file_dict
        self.LON, self.LAT, self.GRID, self.MIN_DEPTH = self.get_lonlat_grid_depth()
        self.LON_bin, self.LAT_bin, self.LON_mid, self.LAT_mid = self.get_lonlat_bins()
        self.depth_bins = self.determine_depth_bins()
        self.temp_direc, self.output_direc = self.get_directories()
        self.season_indices = self.get_season_indices()
        self.output_dict = self.create_output_file_dict()

    def run(self):
        if self.parallel_step == 1:
            # Loading the data
            year, month, run, restart = settings.STARTYEAR, settings.STARTMONTH, settings.RUN, settings.RESTART
            key_year = utils.analysis_simulation_year_key(restart)
            print_statement = 'year {}-{}, run {} restart {}'.format(year, month, run, restart)
            utils.print_statement(print_statement, to_print=True)
            output_name = self.get_file_names(directory=self.temp_direc, final=False)
            if not utils.check_file_exist(output_name, without_pkl=True):
                parcels_dataset = self.load_parcels_output()

                # Looping through the various seasons
                pbar = ProgressBar()
                for season in pbar(self.season_indices.keys()):
                    # Selecting the necessary variables for a particular season
                    season_dict = {}
                    ind_min, ind_max = self.season_indices[season][0], self.season_indices[season][1]
                    for variable in self.set_data_variable_list():
                        season_dict[variable] = parcels_dataset.variables[variable][:, ind_min:ind_max]

                    # Looping through the various lonlat bins and
                    for lo_i in range(self.LON_bin.size - 1):
                        for la_i in range(self.LAT_bin.size - 1):
                            selection = (season_dict['lon'] > self.LON_bin[lo_i]) & \
                                        (season_dict['lon'] < self.LON_bin[lo_i + 1]) & \
                                        (season_dict['lat'] > self.LAT_bin[la_i]) & \
                                        (season_dict['lat'] > self.LAT_bin[la_i + 1]) & \
                                        (season_dict['beach'] == 0)
                            if selection.max() > 0:
                                depth_selection = season_dict['z'][selection]
                                depth_counts, _ = np.histogram(depth_selection, bins=self.depth_bins)
                                self.output_dict[key_year][season][(self.LON_mid[lo_i], self.LAT_mid[la_i])] += depth_counts
                            # If there are no particles within this region, remove this key from the output dict to
                            # reduce storage requirements
                            else:
                                self.output_dict[key_year][season].pop((self.LON_mid[lo_i], self.LAT_mid[la_i]))

            utils.save_obj(output_name, self.output_dict)
            str_format = settings.STARTYEAR, settings.STARTMONTH, settings.RUN, settings.RESTART
            print_statement = 'The concentration for year {}-{}, run {} restart {} has been save'.format(*str_format)
            utils.print_statement(print_statement, to_print=True)

        elif self.parallel_step == 2:
            for run in range(0, settings.RUN_RANGE):
                for restart in range(0, settings.SIM_LENGTH):
                    file_name = self.get_file_names(final=False, directory=self.temp_direc, run=run, restart=restart)
                    dataset_post = utils.load_obj(filename=file_name)
                    for year in range(settings.SIM_LENGTH):
                        key_year = utils.analysis_simulation_year_key(year)
                        for season in self.season_indices.keys():
                            for locations in dataset_post[key_year][season].keys():
                                if type(locations) == tuple:
                                    self.output_dict[key_year][season][locations] += dataset_post[key_year][season][locations]
                    utils.remove_file(file_name + '.pkl')

            utils.save_obj(self.get_file_names(final=True, directory=self.output_direc), self.output_dict)
            utils.print_statement("The spatial vertical concentration has been saved", to_print=True)
        else:
            ValueError('settings.PARALLEL_STEP can not have a value of {}'.format(self.parallel_step))

    def get_file_names(self, directory, final, run=settings.RUN, restart=settings.RESTART):
        split = {True: None, False: '.nc'}[final]
        prefix = 'spatial_vertical_profile'
        output_name = directory + utils.analysis_save_file_name(input_file=self.file_dict[run][restart], prefix=prefix,
                                                                split=split)
        return output_name

    def get_lonlat_bins(self):
        lon_bin = np.arange(np.round(self.LON.min()), np.round(self.LON.max()) + 1)
        lat_bin = np.arange(np.round(self.LAT.min()), np.round(self.LAT.max()) + 1)
        lon_mid, lat_mid = (lon_bin[1:] + lon_bin[:-1]) / 2, (lat_bin[1:] + lat_bin[:-1]) / 2
        return lon_bin, lat_bin, lon_mid, lat_mid

    def load_parcels_output(self, run=settings.RUN, restart=settings.RESTART):
        parcels_dataset = Dataset(self.file_dict[run][restart])
        return parcels_dataset

    def create_output_file_dict(self):
        depth_mid = 0.5 * self.depth_bins[1:] + 0.5 * self.depth_bins[:-1]
        output_dict = {'depth_bins': depth_mid}
        # Looping through the years and seasons:
        for year in range(settings.SIM_LENGTH):
            key_year = utils.analysis_simulation_year_key(year)
            output_dict[key_year] = {}
            for season in self.season_indices.keys():
                output_dict[key_year][season] = {}
                for lon_mid in self.LON_mid:
                    for lat_mid in self.LAT_mid:
                        output_dict[key_year][season][(lon_mid, lat_mid)] = np.zeros(self.depth_bins.__len__() - 1,
                                                                                     dtype=np.float32)
        return output_dict

    @staticmethod
    def determine_depth_bins():
        advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                            advection_scenario=settings.ADVECTION_DATA,
                                                            repeat_dt=None)
        adv_file_dict = advection_scenario.file_names
        depth_min, depth_max = np.nanmin(adv_file_dict['DEPTH']) - 0.1, np.nanmax(adv_file_dict['DEPTH'])
        if depth_min > 0:
            log_min, log_max = np.log10(depth_min), np.log10(depth_max)
        else:
            log_min, log_max = -1, np.log10(depth_max)
        num = 100
        depth_bins = np.logspace(log_min, log_max, num=num)
        return depth_bins

    @staticmethod
    def get_season_indices():
        # Setting the start time of the simulation, and defining the output dt and time index
        current_time = datetime(settings.STARTYEAR + settings.RESTART, settings.STARTMONTH, settings.STARTDAY)
        end_time = datetime(settings.STARTYEAR + settings.RESTART + 1, 1, 1)
        output_dt, index_dt = settings.OUTPUT_TIME_STEP, 0
        # Creating the dictionary containing the start and end indices for each month
        month_indices = {}
        while current_time < end_time:
            if current_time.month not in month_indices.keys():
                month_indices[current_time.month] = [index_dt, index_dt]
            else:
                month_indices[current_time.month][1] = index_dt
            current_time += output_dt
            index_dt += 1
        # And now then the indices for the seasons, which we do in 3-monthly blocks
        season_indices = {0: (month_indices[1][0], month_indices[3][1]), 1: (month_indices[4][0], month_indices[6][1]),
                          2: (month_indices[7][0], month_indices[9][1]), 3: (month_indices[10][0], month_indices[12][1])
                          }
        return season_indices

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
        return {'beach': 1, 'adrift': 0, 'seabed': 3}

    @staticmethod
    def set_data_variable_list():
        return ['lon', 'lat', 'beach', 'z']
