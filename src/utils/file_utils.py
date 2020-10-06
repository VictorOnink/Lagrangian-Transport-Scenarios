from operator import attrgetter

import src.settings as settings
from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
from parcels import JITParticle, Variable

def _get_data_directory(server: int) -> str:
    """

    :param server:
    :return:
    """
    return settings.DATA_DIR_SERVERS[server]


def _get_input_directory(server: int) -> str:
    """

    :param server:
    :return:
    """
    return settings.DATA_INPUT_DIR_SERVERS[server]

def _get_start_end_time(time:str):
    start_time = datetime(settings.START_YEAR + settings.RESTART, 1, 1, 0, 0)
    end_time = datetime(settings.START_YEAR + settings.RESTART + 1, 1, 1, 0, 0)
    simulation_length = (end_time - start_time).days
    if time=='start':
        return start_time
    elif time=='end':
        return end_time
    elif time=='length':
        return simulation_length

def _nan_removal(dataset: Dataset, variable: str, last_selec: np.array,
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

def _get_repeat_dt():
    if settings.RESTART == 0:
        repeat_dt = timedelta(days=31)
    else:
        repeat_dt = None
    return repeat_dt

def _add_var_particle(particleType: JITParticle, name: str, dtype=np.int32,
                      set_initial: bool = True):
    if set_initial == True:
        init = attrgetter(name)
    else:
        init = 0
    var = Variable(name, dtype=dtype, initial=init)
    setattr(particleType, name, var)
