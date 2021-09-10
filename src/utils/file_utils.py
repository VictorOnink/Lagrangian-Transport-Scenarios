import pickle
import settings as settings
from datetime import datetime
from netCDF4 import Dataset
import numpy as np
import os
import shutil

import utils


def get_data_directory(server: int) -> str:
    """

    :param server:
    :return:
    """
    return settings.DATA_DIR_SERVERS[server]


def get_input_directory(server: int) -> str:
    """

    :param server:
    :return:
    """
    return settings.DATA_INPUT_DIR_SERVERS[server]


def get_output_directory(server: int) -> str:
    """

    :param server:
    :return:
    """
    return settings.DATA_OUTPUT_DIR_SERVERS[server]


def get_start_end_time(time: str):
    if settings.RESTART == 0:
        start_time = datetime(settings.STARTYEAR, settings.STARTMONTH, settings.STARTDAY, 0, 0)
    else:
        start_time = datetime(settings.STARTYEAR + settings.RESTART, 1, 1, 0, 0)
    end_time = datetime(settings.STARTYEAR + settings.RESTART + 1, 1, 1, 0, 0)
    simulation_length = (end_time - start_time).days
    if time == 'start':
        return start_time
    elif time == 'end':
        return end_time
    elif time == 'length':
        return simulation_length


def restart_nan_removal(dataset: Dataset, variable: str, last_selec: np.array,
                        final_time: datetime, last_time_selec: datetime):
    """
    This function inputs a dataset object for the rfile. We then take the last
    non-masked value for each row, which we then use to initialise the new ofile
    run. However, in some cases a particle has been deleted during the rfile run,
    and while that particle does stay deleted, we do need to incorporate it in
    the ofile run so that we don't misalign the rows. Therefore, in cases where
    a particle is deleted, we return varSelec with 2 for all those particles,
    since particles where particle.beach==2 will not be advected or be resuspended.

    Parameters
    ----------
    dataset : netcdf4 dataset object
        the dataset object we get the variable field from.
    variable : string
        name of the variable we wish to examine.
    lastSelec : int
        For each row, the index of the last unmasked cell. If there are no
        masked cells, it indicates the last cell of the .
    finalTime : datetime object
        The last timestep of the previous restart file.
    lastTimeSelec : datetime object
        the time of the last unmasked call, as indicated by the index from
        lastSelec.

    Returns
    -------
    varSelec : array Nx1
        The restart array for the given variable to start up the ofile run.

    """
    var = np.array(dataset.variables[variable][:])
    var_selec = var[last_selec[0], last_selec[1]]
    var_selec[last_time_selec != final_time] = 2
    return var_selec


def check_direc_exist(direc: str):
    if not os.path.isdir(direc):
        os.makedirs(direc)


def remove_directory(direc: str):
    if check_direc_exist(direc):
        shutil.rmtree(direc)


def remove_file(File: str, conduct: bool=True):
    if conduct:
        if check_file_exist(File):
            os.remove(File)
        else:
            utils.print_statement('The file {} does not exist.'.format(File), to_print=True)


def check_file_exist(File: str):
    return os.path.isfile(File)


def save_obj(filename, item):
    if filename[-4:] != '.pkl':
        filename += '.pkl'
    with open(filename, 'wb') as f:
        pickle.dump(item, f, pickle.HIGHEST_PROTOCOL)


def load_obj(filename):
    if filename[-4:] != '.pkl':
        filename += '.pkl'
    with open(filename, 'rb') as f:
        return pickle.load(f)


def create_list(var, length):
    return [var] * length


def print_statement(statement, to_print=False):
    if settings.SUBMISSION in ['simulation'] or to_print is True:
        os.system('echo {}'.format(statement))
