from operator import attrgetter
import pickle
import settings as settings
from datetime import datetime
from netCDF4 import Dataset
import numpy as np
from parcels import JITParticle, Variable
import os


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
    start_time = datetime(settings.START_YEAR + settings.RESTART, 1, 1, 0, 0)
    end_time = datetime(settings.START_YEAR + settings.RESTART + 1, 1, 1, 0, 0)
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


def get_repeat_dt():
    if settings.RESTART == 0:
        repeat_dt = settings.REPEAT_DT_R0
    else:
        repeat_dt = settings.REPEAT_DT_ELSE
    return repeat_dt


def add_particle_variable(particleType: JITParticle, name: str, other_name = None, other_value: str = None,
                          dtype=np.int32, set_initial: bool = True, to_write: bool = True):
    if set_initial:
        if other_name is None and other_value is None:
            init = attrgetter(name)
        elif other_name is None and other_value is not None:
            init = other_value
        else:
            init = attrgetter(other_name)
    else:
        init = 0
    var = Variable(name, dtype=dtype, initial=init, to_write=to_write)
    setattr(particleType, name, var)


def check_direc_exist(direc: str):
    if not os.path.isdir(direc):
        os.mkdir(direc)


def check_file_exist(File: str):
    return os.path.isfile(File)


def save_obj(filename, item):
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(item, f, pickle.HIGHEST_PROTOCOL)


def load_obj(filename):
    with open(filename + '.pkl', 'rb') as f:
        return pickle.load(f)
