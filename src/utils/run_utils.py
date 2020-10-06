from datetime import time
import parcels.rng as rng

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
        rng.seed(11235811)
    elif seed == 'TimeSeed':
        rng.seed(int(time())*1000000)
    else:
        rng.seed(int(seed))