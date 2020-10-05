import os
import src.settings as settings
import numpy as np
from parcels import JITParticle, Variable
from operator import attrgetter
from src.factories.fieldset_factory import _get_input_directory
from datetime import timedelta, datetime
import src.scenarios.coastal_proximity as proximity

class ParticleFactory:
    """A factory class for the particles"""
    @classmethod
    def create_particle_set(cls,server: int, fieldset: FieldSet, distance=True):
        """

        :rtype: object
        """
        input_dir = _get_input_directory(server=server)

        #Create the particle class
        os.system('echo "Create the particle class"')
        particle_type = BaseParticle
        if settings.SCENARIO_NAME=='CoastalProximity':
            _add_var_particle(particle_type,'prox')
        if distance:
            _add_var_particle(particle_type,'distance',dtype=np.float32,set_initial=False)
        #Get the particle start positions and other relevant variables
        os.system('echo "Get the (re)start particle positions and other relevant variables"')
        rfile=_get_rfile(input_dir)
        if settings.RESTART==0:
            var_dict = _particle_input_new_simulation(input_dir)
        else:
            var_dict = _get_restart_variables(rfile)
        #Get timestep of repeated particle releases and the starting date of the simulation
        os.system('echo "Set the timestep of particle release and simulation start time"')
        repeat_dt=_get_repeat_dt()
        start_time, _, _ =_get_start_end_time()
        #Create the pset
        pset = _create_pset(fieldset, particle_type, var_dict, start_time, repeat_dt)
        return pset

class BaseParticle(JITParticle):
    #0=open ocean, 1=beached
    beach = Variable('beach', dtype=np.int32, initial=attrgetter('beach'))
    #Finally, I want to keep track of the age of the particle
    age = Variable('age', dtype=np.int32, initial=attrgetter('age'))
    #Weight of the particle in tons
    weights = Variable('weights', dtype=np.float32, initial=attrgetter('weights'))

 def _add_var_particle(particleType: JITParticle, name: str, dtype = np.int32,
                       set_initial: bool =True):
     if set_initial==True:
         init=attrgetter(name)
     else:
         init=0
     var = Variable(name, dtype = dtype, initial = init)
     setattr(particleType, name, var)

def _nan_removal(dataset: netcdf_dataset, variable: str, last_selec: array,
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

def _get_rfile(input_dir: str):
    if settings.SCENARIO_NAME=='CoastalProximity':
        _,rfile=proximity._file_names_proximity(input_dir)
    return rfile

def _get_restart_variables(rfile):
    """
    Parameters
    ----------
    rfile : string
        restart file name.
    restart : int
        Which restart this is
    scenario : int
        which beaching scenario. 0 = Advection & Diffusion Only
                                 1 = Coastal vicinity
                                 2 = Stochastic
                                 3 = Shore dependent resuspension

    Returns
    -------
    TYPE
        lon_reset,lat_reset,beach_reset,age_reset,prox_reset,weight_reset
        restart arrays for longitude, latitude, beached status, particle weight,
        particle age, particle proximity to shore (for scenario 1 only)

    """
    if settings.SCENARIO_NAME=='CoastalProximity':
        var_dict=proximity._get_restart_variables_proximity(rfile)
    #Returning the dictionary
    return var_dict

def _particle_input_new_simulation(input_dir: str):
    if settings.SCENARIO_NAME=='CoastalProximity':
        var_dict=proximity._particle_input_new_simulation_proximity(input_dir)
    return var_dict

def _get_repeat_dt():
    if settings.RESTART==0:
        repeat_dt=timedelta(days=31)
    else:
        repeat_dt=None
    return repeat_dt

def _get_start_end_time():
    start_time = datetime(settings.START_YEAR + settings.RESTART, 1, 1, 0, 0)
    end_time = datetime(settings.START_YEAR + settings.RESTART + 1, 1, 1, 0, 0)
    simulation_length = (end_time - start_time).days
    return start_time, end_time, simulation_length

def _create_pset(fieldset: FieldSet, particle_type: pclass, var_dict: dict,
                           start_time: datetime, repeat_dt: timedelta):
    if settings.SCENARIO_NAME=='CoastalProximity':
        pset=proximity._create_pset_proximity(fieldset, particle_type, var_dict, start_time, repeat_dt)
    return pset
