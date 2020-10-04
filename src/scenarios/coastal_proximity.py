from parcels import FieldSet, JITParticle

import src.settings as settings
import src.scenarios.base_scenario as base
from src.factories.fieldset_factory import FieldSetFactory
from src.factories.particle_factory import ParticleFactory


class CoastalProximity(base.BaseScenario):
    """Coastal proximity scenario"""

    def __init__(self, server, stokes):
        """Constructor for coastal_proximity"""
        super().__init__(server, stokes)

    def create_fieldset(self) -> FieldSet:
        fieldset = FieldSetFactory().create_fieldset(server=self.server, stokes=self.stokes,
                                                     border_current=True, diffusion=True,
                                                     landID=True,distance=True,vicinity=True)
        return fieldset


    def create_particle(self) -> JITParticle:
        particle = ParticleFactory().create_particle()
        return particle

    def run(self) -> None:
        pass

def _file_names_proximity(input_dir: str):
    prefix=base.prefix_determinant(settings.SCENARIO_NAME,settings.STOKES)
    odirec = input_dir + "coastal_v_" + str(settings.VICINITY) + "_e_" + str(settings.ENSEMBLE) + "/"
    ofile = odirec + prefix + "_v=" + str(settings.VICINITY) + "_y=" + str(settings.START_YEAR) + "_I=" + \
            str(settings.INPUT) + "_r=" + str(settings.RESTART) + "_run=" + str(settings.RESTART) + ".nc"
    rfile = odirec + prefix + "_v=" + str(settings.VICINITY) + "_y=" + str(settings.START_YEAR) + "_I=" + \
            str(settings.INPUT) + "_r=" + str(settings.RESTART-1) + "_run=" + str(settings.RESTART) + ".nc"
    return ofile, rfile


def _create_pset_proximity(fieldset: FieldSet, particle_type: pclass, var_dict: dict,
                           start_time: datetime, repeat_dt: timedelta):
    """

    :return:
    """
    pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                       lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                       age=var_dict['age'], prox=var_dict['prox'], weights=var_dict['weights'],
                       time=start_time, repeatdt=repeat_dt)
    return pset

def _get_restart_variables_proximity(rfile: str):
    dataset = Dataset(rfile)
    time = dataset.variables['time'][:]
    final_time = time[0, -1]
    last_selec = np.ma.notmasked_edges(time, axis=1)[1]
    last_time_selec = time[last_selec[0], last_selec[1]]
    # Creating the list of variables and the dictionary containing the arrays for the restart
    var_list = ['lon', 'lat', 'beach', 'age', 'weights','prox']
    var_dict = {}
    for var in var_list:
        var_dict[var] = _nan_removal(dataset, var, last_selec, final_time, last_time_selec)
    return var_dict

def _particle_input_new_simulation_proximity(input_dir: str):
    var_dict = {}
    if settings.INPUT == 'Jambeck':
        var_dict['lon'] = np.load(input_dir + 'Jambeck2010/Jam' + str(2010) + 'Lons' + str(settings.RUN) + '.npy')
        var_dict['lat'] = np.load(input_dir + 'Jambeck2010/Jam' + str(2010) + 'Lats' + str(settings.RUN) + '.npy')
        var_dict['weights'] = np.load(input_dir + 'Jambeck2010/Jam' + str(2010) + 'Weight' + str(settings.RUN) + '.npy')
    elif settings.INPUT == 'Lebreton':
        var_dict['lon'] = np.load(input_dir + 'Lebreton2010/Leb' + str(2010) + 'Lons' + str(settings.RUN) + '.npy')
        var_dict['lat'] = np.load(input_dir + 'Lebreton2010/Leb' + str(2010) + 'Lats' + str(settings.RUN) + '.npy')
        var_dict['weights'] = np.load(input_dir + 'Lebreton2010/Leb' + str(2010) + 'Weight' + str(settings.RUN) + '.npy')
    var_dict['beach'] = np.zeros(len(lons), dtype=np.int32)
    var_dict['age'] = np.zeros(len(lons), dtype=np.int32)
    var_dict['prox'] = np.zeros(len(lons), dtype=np.int32)
    return var_dict

def _beaching_kernel_proximity(particle,fieldset,time):
    if particle.beach==0:
        dist=fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
        #If a particle is within 10 km of the shore
        if dist<10:
            particle.prox+=particle.dt
        else:
            particle.prox=0.
        if particle.prox>86400*fieldset.vic:
            particle.beach=1
    #Update the age of the particle
    particle.age+=particle.dt

def _particle_behavior_proximity(pset: ParticleSet):
    base_behavior = pset.Kernel(base._initial_input) + pset.Kernel(base._floating_advection_rk4) + \
                 pset.Kernel(base._floating_2d_brownian_motion)
    total_behavior = base_behavior + pset.Kernel(base._anti_beach_nudging) + pset.Kernel(_beaching_kernel_proximity)
    return total_behavior