from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
from progressbar import ProgressBar
import settings
import utils
import glob


class parcels_to_timeslicing:
    def __init__(self, file_dict: dict):
        self.parallel_step = settings.PARALLEL_STEP
        self.file_dict = file_dict
        self.temp_direc, self.output_direc = get_directories(scenario_name=settings.SCENARIO_NAME)
        self.reference_time = datetime(2010, 1, 1, 12, 0)
        self.time_difference_list = set_time_difference_list(ref_time=self.reference_time, parallel_step=self.parallel_step)
        self.variable_list = ['lon', 'lat', 'z', 'beach', 'size_class', 'time']

    def run(self):
        if self.parallel_step == 1:
            year, month, run, restart = settings.STARTYEAR, settings.STARTMONTH, settings.RUN, settings.RESTART
            utils.print_statement('Splitting {}-{}, run {} restart {}'.format(year, month, run, restart), to_print=True)
            # Loading the data
            parcels_dataset = load_parcels_post_output(file_dict=self.file_dict)
            full_dict = {}
            for key in self.variable_list:
                if key in parcels_dataset.variables.keys():
                    full_dict[key] = parcels_dataset.variables[key][:, :-1]
            # Looping through the dates
            for time_value in self.time_difference_list:
                date = (self.reference_time + timedelta(seconds=time_value)).strftime("%Y-%m-%d-%H-%M-%S")
                prefix = 'timeslices_{}'.format(date)
                output_name = get_file_names(file_dict=self.file_dict, prefix=prefix, directory=self.temp_direc,
                                             final=False)
                if not utils.check_file_exist(output_name, without_pkl=True):
                    time_selection = full_dict['time'] == time_value
                    # Setting the output name
                    if np.nansum(time_selection) > 0:
                        time_dict = {}
                        for variable in self.variable_list[:-1]:
                            if variable in full_dict.keys():
                                time_dict[variable] = full_dict[variable][time_selection]
                        utils.save_obj(filename=output_name, item=time_dict)
        elif self.parallel_step == 2:
            pbar = ProgressBar()
            for time_value in pbar(self.time_difference_list):
                date = (self.reference_time + timedelta(seconds=time_value)).strftime("%Y-%m-%d-%H-%M-%S")
                prefix = 'timeslices_{}'.format(date)
                # Setting the output name
                output_name = get_file_names(file_dict=self.file_dict, directory=self.output_direc, final=True,
                                             prefix=prefix)
                # Creating an output dict
                output_dict = {}
                for key in self.variable_list[:-1]:
                    output_dict[key] = np.array([], dtype=float)
                # Getting all the timeslice files for this date, and looping through them
                file_list = glob.glob(self.temp_direc + prefix + '*')
                for file_name in file_list:
                    time_file = utils.load_obj(file_name)
                    for key in output_dict.keys():
                        output_dict[key] = np.append(output_dict[key], time_file[key])
                    utils.remove_file(file_name)
                utils.save_obj(output_name, output_dict)
        else:
            ValueError('settings.PARALLEL_STEP can not have a value of {}'.format(self.parallel_step))


########################################################################################################################
"""
These following functions are used across all scenarios
"""


def get_directories(scenario_name):
    # temp_direc = utils.get_output_directory(server=settings.SERVER) + 'timeslices/{}/temporary/'.format(scenario_name)
    temp_direc = settings.SCRATCH_DIR
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeslices/{}/'.format(scenario_name)
    utils.check_direc_exist(temp_direc)
    utils.check_direc_exist(output_direc)
    return temp_direc, output_direc


def load_parcels_post_output(file_dict, year=settings.STARTYEAR, month=settings.STARTMONTH,
                             run=settings.RUN, restart=settings.RESTART, with_post=True):
    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
        parcels_dataset = Dataset(file_dict['parcels'][year][month][run][restart])
    else:
        parcels_dataset = Dataset(file_dict[run][restart])
    return parcels_dataset


def get_file_names(file_dict, directory, final, prefix, year=settings.STARTYEAR, month=settings.STARTMONTH,
                   run=settings.RUN, restart=settings.RESTART):
    split = {True: None, False: '.nc'}[final]
    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict['postprocess'][year][month][run][restart],
        prefix=prefix, split=split)
    else:
        output_name = directory + utils.analysis_save_file_name(input_file=file_dict[run][restart], prefix=prefix, split=split)
    return output_name


def set_time_difference_list(ref_time, parallel_step):
    if parallel_step == 1:
        current_time, end_time = datetime(2010 + settings.RESTART, 1, 1, 0), datetime(2010 + settings.RESTART + 1, 1, 1, 0)
    elif parallel_step == 2:
        current_time, end_time = datetime(2010, 1, 1, 0), datetime(2010 + settings.SIM_LENGTH, 1, 1, 0)
    time_step, time_list = timedelta(hours=24), []
    while current_time < end_time:
        time_list.append((current_time - ref_time).total_seconds())
        current_time += time_step
    return time_list
