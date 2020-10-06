from operator import attrgetter

from parcels import FieldSet, ParticleSet, Variable, JITParticle
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
import utils as utils
from datetime import datetime, timedelta
import os

class CoastalProximity(base_scenario.BaseScenario):
    """Coastal proximity scenario"""
    def __init__(self, server, stokes):
        """Constructor for coastal_proximity"""
        super().__init__(server, stokes)
        self.prefix = "Prox"
        self.input_dir = utils._get_input_directory(server=self.server)

    var_list = ['lon', 'lat', 'beach', 'age', 'weights', 'prox']

    def create_fieldset(self) -> FieldSet:
        os.system('echo "Creating the fieldset"')
        fieldset = fieldset_factory.FieldSetFactory().create_fieldset(server=self.server, stokes=self.stokes,
                                                                      border_current=True, diffusion=True,
                                                                      landID=True, distance=True, vicinity=True)
        return fieldset

    def _get_pset(self, fieldset: FieldSet, particle_type: utils.BaseParticle, var_dict: dict,
                   start_time: datetime, repeat_dt: timedelta):
        """

        :return:
        """
        os.system('echo "Creating the particle set"')
        pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                           lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                           age=var_dict['age'], prox=var_dict['prox'], weights=var_dict['weights'],
                           time=start_time, repeatdt=repeat_dt)
        return pset

    def _get_pclass(self):
        os.system('echo "Creating the particle class"')
        # particle_type = utils.BaseParticle
        # utils._add_var_particle(particle_type, 'prox')
        # utils._add_var_particle(particle_type, 'distance', dtype=np.float32, set_initial=False)

    class particle_type(JITParticle):
        # First we keep track of how long a particle has been close to the shore
        prox = Variable('prox', dtype=np.int32, initial=attrgetter('prox'))
        # Now the beaching variables
        # 0=open ocean, 1=beached
        beach = Variable('beach', dtype=np.int32,
                         initial=attrgetter('beach'))
        # Finally, I want to keep track of the age of the particle
        age = Variable('age', dtype=np.int32, initial=attrgetter('age'))
        # Weight of the particle in tons
        weights = Variable('weights', dtype=np.float32, initial=attrgetter('weights'))
        # Distance of the particle to the coast
        distance = Variable('distance', dtype=np.float32, initial=0)

    def _file_names(self, new: bool = False):
        odirec = self.input_dir + "coastal_v_" + str(settings.VICINITY) + "_e_" + str(settings.ENSEMBLE) + "/"
        if new==True:
            os.system('echo "Set the output file name"')
            return odirec + self.prefix + "_v=" + str(settings.VICINITY) + "_y=" + str(settings.START_YEAR) + "_I=" + \
                    str(settings.INPUT) + "_r=" + str(settings.RESTART) + "_run=" + str(settings.RESTART) + ".nc"
        else:
            os.system('echo "Set the restart file name"')
            return odirec + self.prefix + "_v=" + str(settings.VICINITY) + "_y=" + str(settings.START_YEAR) + "_I=" + \
                    str(settings.INPUT) + "_r=" + str(settings.RESTART - 1) + "_run=" + str(settings.RESTART) + ".nc"

    # def _get_var_dict(self):
    #     if settings.RESTART==0:
    #         return self._get_var_dict()
    #     else:
    #         return self._get_restart_variables(rfile=self.self._file_names(new=False),var_list=self.var_list)

    def _beaching_kernel(particle, fieldset, time):
        if particle.beach == 0:
            dist = fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
            # If a particle is within 10 km of the shore
            if dist < 10:
                particle.prox += particle.dt
            else:
                particle.prox = 0.
            if particle.prox > 86400 * fieldset.vic:
                particle.beach = 1
        # Update the age of the particle
        particle.age += particle.dt


    def _get_particle_behavior(self, pset: ParticleSet):
        os.system('echo "Setting the particle behavior"')
        base_behavior = pset.Kernel(utils._initial_input) + pset.Kernel(utils._floating_advection_rk4) + \
                        pset.Kernel(utils._floating_2d_brownian_motion)
        total_behavior = base_behavior + pset.Kernel(utils._anti_beach_nudging) + pset.Kernel(self._beaching_kernel)
        return total_behavior
