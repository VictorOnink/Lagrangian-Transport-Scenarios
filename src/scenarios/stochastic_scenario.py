from parcels import FieldSet, ParticleSet
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
import utils as utils
from datetime import datetime, timedelta
import os
import parcels.ParcelsRandom as random
import math


class Stochastic(base_scenario.BaseScenario):
    """Stochastic beaching and resuspension scenario"""

    def __init__(self, server, stokes):
        """Constructor for coastal_proximity"""
        super().__init__(server, stokes)
        self.prefix = "Stochastic"
        self.input_dir = utils._get_input_directory(server=self.server)
        self.output_dir = utils._get_output_directory(server=self.server)

    var_list = ['lon', 'lat', 'beach', 'age', 'weights']

    def create_fieldset(self) -> FieldSet:
        os.system('echo "Creating the fieldset"')
        fieldset = fieldset_factory.FieldSetFactory().create_fieldset(server=self.server, stokes=self.stokes,
                                                                      border_current=True, diffusion=True, landID=True,
                                                                      distance=True,beach_timescale=True,resus_timescale=True
                                                                      )
        return fieldset

    def _get_pset(self, fieldset: FieldSet, particle_type: utils.BaseParticle, var_dict: dict,
                  start_time: datetime, repeat_dt: timedelta):
        """
        :return:
        """
        os.system('echo "Creating the particle set"')
        pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                           lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                           age=var_dict['age'], weights=var_dict['weights'],
                           time=start_time, repeatdt=repeat_dt)
        return pset

    def _get_pclass(self):
        os.system('echo "Creating the particle class"')
        particle_type = utils.BaseParticle
        utils._add_var_particle(particle_type, 'distance', dtype=np.float32, set_initial=False)
        return particle_type

    def _file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART):
        odirec = self.output_dir + "st_" + str(settings.SHORE_TIME) + "_rs_" + str(settings.RESUS_TIME) + \
                 "_e_" + str(settings.ENSEMBLE) + "/"
        if new == True:
            os.system('echo "Set the output file name"')
            return odirec + self.prefix + "_st=" + str(settings.SHORE_TIME) + "_rt=" + str(settings.RESUS_TIME) + \
                   "_y=" + str(settings.START_YEAR) + "_I=" + str(settings.INPUT) + "_r=" + str(restart) + \
                   "_run=" + str(run) + ".nc"
        else:
            os.system('echo "Set the restart file name"')
            return odirec + self.prefix + "_st=" + str(settings.SHORE_TIME) + "_rt=" + str(settings.RESUS_TIME) + \
                   "_y=" + str(settings.START_YEAR) + "_I=" + str(settings.INPUT) + "_r=" + str(restart - 1) + \
                   "_run=" + str(run) + ".nc"

    def _beaching_kernel(particle, fieldset, time):
        if particle.beach == 0:
            dist = fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
            if dist < fieldset.Coastal_Boundary:
                if random.random() > fieldset.p_beach:
                    particle.beach = 1
        # Now the part where we build in the resuspension
        elif particle.beach == 1:
            if random.random() > fieldset.p_resus:
                particle.beach = 0
        # Update the age of the particle
        particle.age += particle.dt

    def _get_particle_behavior(self, pset: ParticleSet):
        os.system('echo "Setting the particle behavior"')
        base_behavior = pset.Kernel(utils._initial_input) + pset.Kernel(utils._floating_advection_rk4) + \
                        pset.Kernel(utils._floating_2d_brownian_motion)
        total_behavior = base_behavior + pset.Kernel(utils._anti_beach_nudging) + pset.Kernel(self._beaching_kernel)
        return total_behavior
