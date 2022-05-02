import pickle
import settings as settings
from datetime import datetime
from netCDF4 import Dataset
import numpy as np
import os
import shutil
import utils


def get_start_end_time(time: str) -> datetime:
    """
    Depending the input variable time, returns either the datetime object for the beginning or end of the simulation, or
    the length of the simulation in days.
    :param time: specifies what datetime object to return
    :return:
    """
    if settings.RESTART == 0:
        start_time = datetime(settings.STARTYEAR, settings.STARTMONTH, settings.STARTDAY, 0, 0)
    else:
        start_time = datetime(settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT, 1, 1, 0, 0)
    end_time = datetime(settings.STARTYEAR + (settings.RESTART + 1) * settings.BACKWARD_MULT, 1, 1, 0, 0)
    simulation_length = abs((end_time - start_time).days)
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
    """
    Check if the directory specified by the path direc exists, and if it doesn't, create the directory.
    :param direc: path of the directory
    :return:
    """
    if not os.path.isdir(direc):
        os.makedirs(direc)


def remove_directory(direc: str):
    """
    Check if a directory exists and if it does, delete it.
    :param direc: path of the directory
    :return:
    """
    if check_direc_exist(direc):
        shutil.rmtree(direc)


def remove_file(File: str, conduct: bool = True):
    """
    Check if the file with path File exists, and if it does then delete it
    :param File: path of the file
    :param conduct: boolean that determines where the remove_file procedure is conducted
    :return:
    """
    if conduct:
        if check_file_exist(File):
            os.remove(File)
        else:
            utils.print_statement('The file {} does not exist.'.format(File), to_print=True)


def check_file_exist(File: str, without_pkl=False):
    """
    Check if a file exists
    :param File: path of the file
    :param without_pkl: if false, append '.pkl' to the end of the file path
    :return:
    """
    if without_pkl:
        File += '.pkl'
    if os.path.isfile(File):
        return os.path.isfile(File)
    else:
        utils.print_statement("{} doesn't exist".format(File))
        return os.path.isfile(File)


def save_obj(filename, item):
    """
    Pickles an object at the location specified by the filename path
    :param filename: path of the location where pickle object is to be saved
    :param item: object being pickled
    :return:
    """
    if filename[-4:] != '.pkl':
        filename += '.pkl'
    with open(filename, 'wb') as f:
        pickle.dump(item, f, pickle.HIGHEST_PROTOCOL)


def load_obj(filename):
    """
    Load the pickled object within the file specified by the filename path
    :param filename: path of the pickle file
    :return:
    """
    if filename[-4:] != '.pkl':
        filename += '.pkl'
    with open(filename, 'rb') as f:
        return pickle.load(f)


def print_statement(statement, to_print=False):
    """

    :param statement: string of the
    :param to_print:
    :return:
    """
    if settings.SUBMISSION in ['simulation'] or to_print is True:
        os.system('echo {}'.format(statement))


def flatten_list_of_lists(t):
    return [item for sublist in t for item in sublist]
