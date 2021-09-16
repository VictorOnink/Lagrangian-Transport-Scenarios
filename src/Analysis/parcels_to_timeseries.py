import settings as settings
import utils
from netCDF4 import Dataset
import numpy as np
from copy import deepcopy
from progressbar import ProgressBar
from datetime import datetime, timedelta


class parcels_to_timeseries:
    def __init__(self, file_dict: dict):
        self.parallel_step = settings.PARALLEL_STEP
        self.file_dict = file_dict
        self.domain = (-180, 180, -90, 90) # lon_min, lon_max, lat_min, lat_max
        self.temp_direc, self.output_direc = get_directories(scenario_name=settings.SCENARIO_NAME)
        self.beach_label_dict = set_beach_label_dict()
        self.output_dict, self.time_list = create_output_dict_time_list(beach_label_dict=self.beach_label_dict)

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
                full_data_dict = set_full_data_dict(parcels_dataset, post_dataset)
                # Just get the particles within the domain, which we do by setting all values not within the
                # domain to nan. These will therefore not be taken into account in the calculations of total
                # counts/weights
                within_domain = utils.particles_in_domain(domain=self.domain, lon=full_data_dict['lon'],
                                                          lat=full_data_dict['lat'])
                for variable in full_data_dict.keys():
                    full_data_dict[variable][within_domain is False] = np.nan
                # Going through the timesteps
                for time_index, time_value in enumerate(self.time_list):
                    time_selection = full_data_dict['time'] == time_value
                    if np.nansum(time_selection) > 0:
                        time_dict = {}
                        for variable in ['beach', 'weights', 'size_class']:
                            if variable in full_data_dict.keys():
                                time_dict[variable] = full_data_dict[variable][time_selection]
                        for beach_state in self.beach_label_dict.keys():
                            beach = time_dict['beach'] == self.beach_label_dict[beach_state]
                            if 'size_class' in full_data_dict.keys():
                                for size_class in range(settings.SIZE_CLASS_NUMBER):
                                    size = time_dict['size_class'] == size_class
                                    self.output_dict[beach_state][size_class][time_index] += np.nansum(time_dict['weights'][size & beach])
                            else:
                                self.output_dict[beach_state][time_index] += np.nansum(time_dict['weights'][beach])
                # Calculating the total over all beach states
                for beach_state in self.beach_label_dict.keys():
                    if 'size_class' in full_data_dict.keys():
                        for size_class in range(settings.SIZE_CLASS_NUMBER):
                            self.output_dict['total'][size_class] += self.output_dict[beach_state][size_class]
                    else:
                        self.output_dict['total'] += self.output_dict[beach_state]

                # Saving the output
                utils.save_obj(filename=output_name, item=self.output_dict)
                str_format = settings.STARTYEAR, settings.STARTMONTH, settings.RUN, settings.RESTART
                print_statement = 'The timeseries for year {}-{}, run {} restart {} has been save'.format(*str_format)
                utils.print_statement(print_statement, to_print=True)

        elif self.parallel_step == 2:
            pbar = ProgressBar()
            for ind_year, year in pbar(enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH))):
                for month in range(1, 13):
                    for run in range(0, settings.RUN_RANGE):
                        for restart in range(0, settings.SIM_LENGTH - ind_year):
                                file_name = get_file_names(file_dict=self.file_dict,
                                                           directory=self.temp_direc, final=False, year=year,
                                                           month=month,
                                                           run=run, restart=restart)
                                dataset_post = utils.load_obj(filename=file_name)
                                for beach_state in self.output_dict.keys():
                                    if beach_state != 'time':
                                        if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
                                            for size_class in range(settings.SIZE_CLASS_NUMBER):
                                                self.output_dict[beach_state][size_class] += dataset_post[beach_state][size_class]
                                        else:
                                            if month == 1:
                                                self.output_dict[beach_state] += dataset_post[beach_state]
                                utils.remove_file(file_name)
            # Saving the output
            file_name = get_file_names(file_dict=self.file_dict, directory=self.output_direc, final=True)
            utils.save_obj(filename=file_name, item=self.output_dict)
            print_statement = 'The timeseries have been saved'
            utils.print_statement(print_statement, to_print=True)

        else:
            ValueError('settings.PARALLEL_STEP can not have a value of {}'.format(self.parallel_step))


########################################################################################################################
"""
These following functions are used across all scenarios
"""


def get_directories(scenario_name):
    temp_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/{}/temporary/'.format(scenario_name)
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/{}/'.format(scenario_name)
    utils.check_direc_exist(temp_direc)
    utils.check_direc_exist(output_direc)
    return temp_direc, output_direc


def set_beach_label_dict():
    return {'beach': 1, 'adrift': 0, 'seabed': 3, 'removed': 2}


def create_output_dict_time_list(beach_label_dict):
    reference_time = datetime(2010, 1, 1, 12, 0)
    current_time, end_time = datetime(2010, 1, 1, 0), datetime(2010 + settings.SIM_LENGTH, 1, 1, 0)
    time_step, time_list = timedelta(hours=12), []
    while current_time < end_time:
        time_list.append((current_time - reference_time).total_seconds())
        current_time += time_step

    output_dict = {'time': time_list}
    base_array = np.zeros(time_list.__len__(), dtype=float)
    for beach_state in beach_label_dict.keys():
        if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
            output_dict['total'] = {}
            output_dict[beach_state] = {}
            for size_class in range(settings.SIZE_CLASS_NUMBER):
                output_dict[beach_state][size_class] = deepcopy(base_array)
                output_dict['total'][size_class] = deepcopy(base_array)

        else:
            output_dict[beach_state] = deepcopy(base_array)
            output_dict['total'] = deepcopy(base_array)
    return output_dict, time_list


def load_parcels_post_output(scenario_name, file_dict, year=settings.STARTYEAR, month=settings.STARTMONTH,
                             run=settings.RUN, restart=settings.RESTART):
    if scenario_name in ['FragmentationKaandorpPartial']:
        parcels_dataset = Dataset(file_dict['parcels'][year][month][run][restart])
        post_dataset = utils.load_obj(file_dict['postprocess'][year][month][run][restart])
    else:
        parcels_dataset = Dataset(file_dict[run][restart])
        post_dataset = None
    return parcels_dataset, post_dataset


def set_full_data_dict(parcels_dataset, post_dataset):
    full_data_dict = {}
    for variable in ['lon', 'lat', 'beach', 'time', 'size_class']:
        if variable in parcels_dataset.variables.keys():
            full_data_dict[variable] = parcels_dataset.variables[variable][:, :-1].flatten()
    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
        full_data_dict['weights'] = post_dataset['particle_number'][:, :-1].flatten()
    else:
        full_data_dict['weights'] = np.ones(full_data_dict['lon'].shape, dtype=float)
    # Only return the non-nan values
    is_not_nan = ~np.isnan(full_data_dict['lon'])
    for variable in full_data_dict.keys():
        full_data_dict[variable] = full_data_dict[variable][is_not_nan]
    return full_data_dict


def get_file_names(file_dict, directory, final, year=settings.STARTYEAR, month=settings.STARTMONTH,
                   run=settings.RUN, restart=settings.RESTART):
    split = {True: None, False: '.nc'}[final]
    prefix = 'timeseries'
    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict['postprocess'][year][month][run][restart],
        prefix=prefix, split=split)
    else:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix, split=split)
    return output_name
