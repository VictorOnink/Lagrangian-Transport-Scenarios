import src.settings as settings
from parcels import ErrorCode
from src.factories.fieldset_factory import _get_input_directory
from src.factories.particle_factory import _get_start_end_time
from datetime import timedelta
import parcels.rng as rng
from time import time
import src.scenarios.coastal_proximity as proximity
from src.scenarios.base_scenario import DeleteParticle
import os

class ExecuteFactory:
    """A factory class for the execution of the scenario run"""
    def execute_scenario(cls,server: int, pset: ParticleSet, seed: str='Fixed'):
        input_dir = _get_input_directory(server=server)

        os.system('echo "Setting the random seed"')
        _set_random_seed(seed=seed)
        os.system('echo "Defining the particle behavior"')
        behavior_kernel=_set_particle_behavior(pset=pset)
        os.system('echo "Setting the output file"')
        pfile = _get_pfile(pset=pset,input_dir=input_dir)
        os.system('echo "Determine the simulation length"')
        _, _, simulation_length = _get_start_end_time()
        os.system('echo "The actual execution of the run"')
        pset.execute(behavior_kernel,
                     runtime=timedelta(days=simulation_length),
                     dt=settings.TIME_STEP,
                     recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
                     output_file=pfile
                     )
        os.system('echo "Exporting the pfile at the end of the simulation"')
        pfile.export()


def _set_random_seed(seed: str):
    """
    Setting the random seed of the run
    :param seed: Fixed = Then use the default random seed that I have defined
                 TimeSeed = Set the random seed based on the time when simulation starts
                 else: if not either of the previous two, then the input is a string of an
                       integer which will then be used as the random seed
    :return:
    """
    if seed=='Fixed'
        rng.seed(11235811)
    elif seed=='TimeSeed':
        rng.seed(int(time())*1000000)
    else:
        rng.seed(int(seed))

def _set_particle_behavior(pset: ParticleSet):
    if settings.SCENARIO_NAME=='CoastalProximity':
        particle_behavior = proximity._particle_behavior_proximity(pset)
    return particle_behavior

def _get_pfile(pset: ParticleSet,input_dir: str):
    if settings.SCENARIO_NAME=='CoastalProximity':
        ofile,_=proximity._file_names_proximity(input_dir)

    pfile = pset.ParticleFile(name=ofile,
                              outputdt=timedelta(hours=24))
    return pfile