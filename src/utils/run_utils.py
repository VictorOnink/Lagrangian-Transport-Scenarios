from datetime import time
from parcels import ParcelsRandom, JITParticle, Variable
import numpy as np
import os
from operator import attrgetter
import settings


def set_random_seed(seed: str):
    """
    Setting the random seed of the run
    :param seed: Fixed = Then use the default random seed that I have defined
                 TimeSeed = Set the random seed based on the time when simulation starts
                 else: if not either of the previous two, then the input is a string of an
                       integer which will then be used as the random seed
    :return:
    """
    if seed == 'Fixed':
        seed_value = int(11235811)
    elif seed == 'TimeSeed':
        seed_value = int(time())*1000000
    else:
        seed_value = int(seed)
    os.system('echo "The random seed is "'+str(seed_value))
    ParcelsRandom.seed(seed_value)


def add_particle_variable(particleType: JITParticle, name: str, other_name=None, other_value: str = None,
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


def get_repeat_dt():
    if settings.RESTART == 0:
        repeat_dt = settings.REPEAT_DT_R0
    else:
        repeat_dt = settings.REPEAT_DT_ELSE
    return repeat_dt
