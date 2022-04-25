from parcels import FieldSet, ParticleSet
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
from advection_scenarios import advection_files
import utils as utils
from datetime import datetime, timedelta
import os
import math


class AdvectionDiffusionOnly(base_scenario.BaseScenario):
    """Advection and Diffusion only scenario"""

    def __init__(self, server, stokes):
        """Constructor for AdvectionDiffusionOnly"""
        super().__init__(server, stokes)

    def set_prefix(self) -> str:
        """
        Set the scenario prefix
        :return:
        """
        return "AdvDifOnly"

    def set_var_list(self) -> list:
        """
        Set the var_list, which contains all the variables that need to be loaded during the restarts
        :return:
        """
        return ['lon', 'lat', 'weights', 'beach', 'age']

    def set_time_steps(self) -> tuple:
        """
        Set the integration, output and repeat timesteps
        :return: self.dt, self.output_time_step, self.repeat_dt
        """
        dt = timedelta(minutes=10 * settings.BACKWARD_MULT)
        output_time_step = timedelta(hours=12)
        if settings.RESTART == 0:
            repeat_dt = timedelta(days=31)
        else:
            repeat_dt = None
        return dt, output_time_step, repeat_dt

    def create_fieldset(self) -> FieldSet:
        utils.print_statement("Creating the fieldset")
        fieldset = fieldset_factory.FieldSetFactory().create_fieldset(file_dict=self.file_dict, stokes=self.stokes,
                                                                      diffusion=True, landID=True, distance=True,
                                                                      coastal_zone=False, time_step=self.dt)
        return fieldset

    def get_pset(self, fieldset: FieldSet, particle_type: utils.BaseParticle, var_dict: dict,
                 start_time: datetime, repeat_dt: timedelta):
        """
        :return:
        """
        utils.print_statement("Creating the particle set")
        pset = ParticleSet(fieldset=fieldset, pclass=particle_type, lon=var_dict['lon'], lat=var_dict['lat'],
                           beach=var_dict['beach'], age=var_dict['age'], weights=var_dict['weight'], time=start_time,
                           repeatdt=repeat_dt)
        return pset

    def get_pclass(self):
        utils.print_statement("Creating the particle class")
        particle_type = utils.BaseParticle
        utils.add_particle_variable(particle_type, 'distance', dtype=np.float32, set_initial=False)
        utils.add_particle_variable(particle_type, 'weights', dtype=np.float32, set_initial=True)
        return particle_type

    def file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART):
        odirec = self.output_dir + "AdvDifOnly_e_{}/".format(settings.ENSEMBLE)
        if new:
            os.system('echo "Set the output file name"')
            str_format = (settings.ADVECTION_DATA, settings.STARTYEAR, settings.INPUT, restart, run)
        else:
            os.system('echo "Set the restart file name"')
            str_format = (settings.ADVECTION_DATA, settings.STARTYEAR, settings.INPUT, restart - 1, run)
        return odirec + self.prefix + '_{}_y={}_I={}_r={}_run={}.nc'.format(*str_format)

    def beaching_kernel(particle, fieldset, time):
        # A particle is considered beached if it is within a land cell
        if math.floor(fieldset.landID[time, particle.depth, particle.lat, particle.lon]) == 1:
            particle.beach = 1
        # Update the age of the particle
        particle.age += particle.dt

    def get_particle_behavior(self, pset: ParticleSet):
        utils.print_statement("Setting the particle behavior")
        total_behavior = pset.Kernel(utils.initial_input) + \
                         pset.Kernel(utils.floating_advection_rk4) + \
                         pset.Kernel(utils.floating_2d_brownian_motion) + \
                         pset.Kernel(self.beaching_kernel)
        return total_behavior
