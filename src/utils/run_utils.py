from datetime import time
from parcels import ParcelsRandom as random
import os

def _set_random_seed(seed: str):
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
    random.seed(seed_value)
