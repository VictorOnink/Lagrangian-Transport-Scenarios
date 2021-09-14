from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
from progressbar import ProgressBar
import settings
import utils


class parcels_to_timeslicing:
    def __init__(self, file_dict: dict):
        self.parallel_step = settings.PARALLEL_STEP
        self.file_dict = file_dict
        self.temp_direc, self.output_direc = get_directories(scenario_name=settings.SCENARIO_NAME)
        self.reference_time = datetime(2010, 1, 1, 12, 0)
        self.time_difference_list = set_time_difference_list(ref_time=self.reference_time)
        self.variable_list = ['lon', 'lat', 'z', 'beach', 'size_class', 'time']

    def run(self):
        if self.parallel_step == 1:
            year, month, run, restart = settings.STARTYEAR, settings.STARTMONTH, settings.RUN, settings.RESTART
            utils.print_statement('Splitting {}-{}, run {} restart {}'.format(year, month, run, restart), to_print=True)
            # Loading the data
            parcels_dataset, _ = load_parcels_post_output(file_dict=self.file_dict)
            full_dict = {}
            for key in self.variable_list:
                if key in parcels_dataset.variables.keys():
                    full_dict[key] = parcels_dataset.variables[key][:, :-1]
            # Looping through the dates
            for time_value in self.time_difference_list:
                time_selection = full_dict['time'] == time_value
                if np.nansum(time_selection) > 0:
                    time_dict = {}
                    for variable in self.variable_list[:-1]:
                        if variable in full_dict.keys():
                            time_dict[variable] = full_dict[variable][time_selection]
                    # Setting the output name
                    date = (self.reference_time + timedelta(seconds=time_value)).strftime("%Y-%m-%d-%H-%M-%S")
                    prefix = 'timeslices_{}'.format(date)
                    output_name = get_file_names(file_dict=self.file_dict, prefix=prefix, directory=self.temp_direc,
                                                 final=False)
                    utils.save_obj(filename=output_name, item=time_dict)
        elif self.parallel_step == 2:
            pbar = ProgressBar()
            for ind_year, year in pbar(enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH))):
                for month in range(1, 13):
                    for run in range(0, settings.RUN_RANGE):
                        for restart in range(0, settings.SIM_LENGTH - ind_year):
                            for time_value in self.time_difference_list:
                                date = (self.reference_time + timedelta(seconds=time_value)).strftime(
                                    "%Y-%m-%d-%H-%M-%S")
                                prefix = 'timeslices_{}'.format(date)
                                file_name = get_file_names(file_dict=self.file_dict, directory=self.temp_direc,
                                                           final=False, year=year, month=month, run=run,
                                                           restart=restart, prefix=prefix)
                                date_dict = utils.load_obj(filename=file_name)
                                if utils.check_file_exist(file_name):
                                    output_name = get_file_names(file_dict=self.file_dict, directory=self.output_direc,
                                                                 final=True, prefix=prefix)
                                    if utils.check_file_exist(output_name):
                                        previous_dict = utils.load_obj(filename=output_name)
                                        for key in previous_dict.keys():
                                            previous_dict[key] = np.append(previous_dict[key], date_dict[key])
                                        utils.save_obj(filename=output_name, item=previous_dict)
                                    else:
                                        utils.save_obj(filename=output_name, item=date_dict)
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


def load_parcels_post_output(file_dict, year=settings.STARTYEAR, month=settings.STARTMONTH,
                             run=settings.RUN, restart=settings.RESTART):
    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
        parcels_dataset = Dataset(file_dict['parcels'][year][month][run][restart])
        post_dataset = utils.load_obj(file_dict['postprocess'][year][month][run][restart])
    else:
        parcels_dataset = Dataset(file_dict[run][restart])
        post_dataset = None
    return parcels_dataset, post_dataset


def get_file_names(file_dict, directory, final, prefix, year=settings.STARTYEAR, month=settings.STARTMONTH,
                   run=settings.RUN, restart=settings.RESTART):
    split = {True: None, False: '.nc'}[final]
    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict['postprocess'][year][month][run][restart],
        prefix=prefix, split=split)
    else:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix, split=split)
    return output_name


def set_time_difference_list(ref_time):
    current_time, end_time = datetime(2010, 1, 1, 0), datetime(2010 + settings.SIM_LENGTH, 1, 1, 0)
    time_step, time_list = timedelta(hours=12), []
    while current_time < end_time:
        time_list.append((current_time - ref_time).total_seconds())
        current_time += time_step
    return time_list
