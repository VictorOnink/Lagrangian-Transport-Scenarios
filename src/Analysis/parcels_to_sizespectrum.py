import settings as settings
import utils
from advection_scenarios import advection_files
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar
from copy import deepcopy
from datetime import datetime, timedelta


class parcels_to_sizespectrum:
    def __init__(self, file_dict, scenario):
        assert settings.SCENARIO_NAME == 'FragmentationKaandorpPartial', 'This analysis only works for FragmentationKaandorpPartial'
        self.file_dict = file_dict
        self.parallel_step = settings.PARALLEL_STEP
        self.scenario = scenario
        self.temp_direc, self.output_direc = get_directories()
        self.min_depth = np.nanmin(self.scenario.file_dict['DEPTH'])
        self.surface_depth = self.min_depth + 0.26
        self.time_list = create_time_list()
        self.time_analysis_step = 60
        self.var_list = set_var_list()
        self.beach_label_dict = set_beach_label_dict()
        self.reservoirs = set_reservoirs()
        self.output_dict = create_output_dict(time_list=self.time_list, time_analysis_step=self.time_analysis_step,
                                              reservoir_list=self.reservoirs)

    def run(self):
        if self.parallel_step == 1:
            year, month, run, restart = settings.STARTYEAR, settings.STARTMONTH, settings.RUN, settings.RESTART
            print_statement = 'year {}-{}, run {} restart {}'.format(year, month, run, restart)
            utils.print_statement(print_statement, to_print=True)
            output_name = get_file_names(file_dict=self.file_dict, directory=self.temp_direc, final=False)
            if not utils.check_file_exist(output_name, without_pkl=True):
                # Loading the data
                parcels_dataset, post_dataset = load_parcels_post_output(scenario_name=settings.SCENARIO_NAME,
                                                                         file_dict=self.file_dict)
                full_data_dict = set_full_data_dict(parcels_dataset, post_dataset, self.var_list)
                # Looping through the time
                for index_time in range(0, self.time_list.__len__(), self.time_analysis_step):
                    time_selection = self.time_list[index_time] == full_data_dict['time']
                    if time_selection.max():
                        time_sel_dict = {}
                        for key in full_data_dict.keys():
                            time_sel_dict[key] = full_data_dict[key][time_selection]
                        for reservoir in self.reservoirs:
                            self.output_dict[reservoir][index_time] += reservoir_calculation(reservoir, time_sel_dict,
                                                                                             self.beach_label_dict,
                                                                                             self.surface_depth)
                # Adding the index of the final timestep for ease later on
                self.output_dict['final_index'] = index_time

                # Saving everything
                utils.save_obj(output_name, self.output_dict)
                str_format = settings.STARTYEAR, settings.STARTMONTH, settings.RUN, settings.RESTART
                print_statement = 'The size distribution for year {}-{}, run {} restart {} has been save'.format(*str_format)
                utils.print_statement(print_statement, to_print=True)

        elif self.parallel_step == 2:
            pbar = ProgressBar()
            for ind_year, year in pbar(enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH))):
                for month in range(1, 13):
                    for run in range(0, settings.RUN_RANGE):
                        for restart in range(0, settings.SIM_LENGTH - ind_year):
                            file_name = get_file_names(file_dict=self.file_dict, directory=self.temp_direc, final=False,
                                                       year=year, month=month, run=run, restart=restart)
                            if utils.check_file_exist(file_name, without_pkl=True):
                                dataset_post = utils.load_obj(filename=file_name)
                                for reservoir in self.reservoirs:
                                    for index_time in range(0, self.time_list.__len__(), self.time_analysis_step):
                                        self.output_dict[reservoir][index_time] += dataset_post[reservoir][index_time]
            # Adding the index of the final timestep for ease later on
            self.output_dict['final_index'] = index_time
            # Saving everything
            output_name = get_file_names(file_dict=self.file_dict, directory=self.temp_direc, final=True)
            utils.save_obj(output_name, self.output_dict)
            print_statement = 'The size distribution has been save'
            utils.print_statement(print_statement, to_print=True)




        else:
            ValueError('settings.PARALLEL_STEP can not have a value of {}'.format(self.parallel_step))

########################################################################################################################
"""
These following functions are used across all scenarios
"""


def get_file_names(file_dict, directory, final, year=settings.STARTYEAR, month=settings.STARTMONTH,
                   run=settings.RUN, restart=settings.RESTART):
    split = {True: None, False: '.nc'}[final]
    prefix = 'size_distribution'
    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict['postprocess'][year][month][run][restart],
        prefix=prefix, split=split)
    else:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict[run][restart], prefix=prefix, split=split)
    return output_name


def load_parcels_post_output(scenario_name, file_dict, year=settings.STARTYEAR, month=settings.STARTMONTH,
                             run=settings.RUN, restart=settings.RESTART):
    if scenario_name in ['FragmentationKaandorpPartial']:
        parcels_dataset = Dataset(file_dict['parcels'][year][month][run][restart])
        post_dataset = utils.load_obj(file_dict['postprocess'][year][month][run][restart])
    else:
        parcels_dataset = Dataset(file_dict[run][restart])
        post_dataset = None
    return parcels_dataset, post_dataset


def get_directories():
    temp_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/{}/temporary/'.format(settings.SCENARIO_NAME)
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/{}/'.format(settings.SCENARIO_NAME)
    utils.check_direc_exist(temp_direc)
    utils.check_direc_exist(output_direc)
    return temp_direc, output_direc


def create_time_list():
    reference_time = datetime(2010, 1, 1, 12, 0)
    current_time, end_time = datetime(2010, 1, 1, 0), datetime(2010 + settings.SIM_LENGTH + 1, 1, 1, 0)
    time_step, time_list = timedelta(hours=12), []
    while current_time < end_time:
        time_list.append((current_time - reference_time).total_seconds())
        current_time += time_step
    return time_list


def set_var_list():
    return ['time', 'size_class', 'beach', 'z', 'distance2coast']


def set_beach_label_dict():
    return {'beach': 1, 'adrift': 0}


def set_reservoirs():
    return ['beach', 'adrift', 'adrift_10km', 'adrift_10km_surf', 'adrift_open', 'adrift_open_surf']


def create_output_dict(time_list, time_analysis_step, reservoir_list):
    output_dict = {'size_bins': np.arange(settings.SIZE_CLASS_NUMBER)}
    for reservoir in reservoir_list:
        output_dict[reservoir] = {}
        for time_index in range(0, time_list.__len__(), time_analysis_step):
            output_dict[reservoir][time_index] = np.zeros(shape=settings.SIZE_CLASS_NUMBER, dtype=float)
    return output_dict


def set_full_data_dict(parcels_dataset, post_dataset, variable_list):
    full_data_dict = {}
    for variable in variable_list:
        full_data_dict[variable] = parcels_dataset.variables[variable][:, :-1].flatten()
    full_data_dict['particle_number'] = post_dataset['particle_number'][:, :-1].flatten()
    # Only return the non-nan values
    is_not_nan = ~np.isnan(full_data_dict['beach'])
    for variable in full_data_dict.keys():
        full_data_dict[variable] = full_data_dict[variable][is_not_nan]
    return full_data_dict


def reservoir_calculation(reservoir, time_sel_dict, beach_label_dict, surface_depth):
    if reservoir == 'beach':
        selection = time_sel_dict['beach'] == beach_label_dict['beach']
    elif reservoir == 'adrift':
        selection = time_sel_dict['beach'] == beach_label_dict['adrift']
    elif reservoir == 'adrift_10km':
        selection = (time_sel_dict['beach'] == beach_label_dict['adrift']) & (time_sel_dict['distance2coast'] < 10)
    elif reservoir == 'adrift_10km_surf':
        selection = (time_sel_dict['beach'] == beach_label_dict['adrift']) & (time_sel_dict['distance2coast'] < 10) & \
                    (time_sel_dict['distance2coast'] < surface_depth)
    elif reservoir == 'adrift_open':
        selection = (time_sel_dict['beach'] == beach_label_dict['adrift']) & (time_sel_dict['distance2coast'] > 10)
    elif reservoir == 'adrift_open_surf':
        selection = (time_sel_dict['beach'] == beach_label_dict['adrift']) & (time_sel_dict['distance2coast'] > 10) & (
                    time_sel_dict['distance2coast'] < surface_depth)
    else:
        ValueError('What do you mean by {}'.format(reservoir))
    return number_per_size_class(time_sel_dict['size_class'], time_sel_dict['particle_number'],
                                 settings.SIZE_CLASS_NUMBER, selection=selection)


def number_per_size_class(size_class_array, particle_number_array, bin_number, selection=None):
    assert size_class_array.size == particle_number_array.size, 'The size and number arrays must have the same size'
    output_array = np.zeros(shape=bin_number, dtype=float)
    if selection is not None:
        assert selection.size == size_class_array.size, 'The selection and data arrays must have the same size'
        for size_class in range(bin_number):
            in_size_class = size_class_array[selection] == size_class
            if np.nansum(in_size_class) > 0:
                output_array[size_class] += np.nansum(particle_number_array[selection][in_size_class])
            else:
                output_array[size_class] += 0
    else:
        for size_class in range(bin_number):
            output_array[size_class] += np.nansum(particle_number_array[size_class_array == size_class])
    return output_array
