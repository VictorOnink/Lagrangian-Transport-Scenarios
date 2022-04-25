from parcels import FieldSet, ParticleSet
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
import utils as utils
from datetime import datetime, timedelta
from parcels import ParcelsRandom


class Stochastic(base_scenario.BaseScenario):
    """Stochastic beaching and resuspension scenario"""

    def __init__(self, server, stokes):
        """Constructor for coastal_proximity"""
        super().__init__(server, stokes)

    def set_prefix(self) -> str:
        """
        Set the scenario prefix
        :return:
        """
        return "Stochastic"

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
                                                                      border_current=True, diffusion=True, landID=True,
                                                                      distance=True, beach_timescale=True,
                                                                      resus_timescale=True, fixed_resus=True,
                                                                      time_step=self.dt
                                                                      )
        return fieldset

    def get_pset(self, fieldset: FieldSet, particle_type: utils.BaseParticle, var_dict: dict,
                 start_time: datetime, repeat_dt: timedelta):
        """
        :return:
        """
        utils.print_statement("Creating the particle set")
        pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                           lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                           age=var_dict['age'], weights=var_dict['weight'],
                           time=start_time, repeatdt=repeat_dt)
        return pset

    def get_pclass(self):
        utils.print_statement("Creating the particle class")
        particle_type = utils.BaseParticle
        utils.add_particle_variable(particle_type, 'distance', dtype=np.float32, set_initial=False)
        utils.add_particle_variable(particle_type, 'weights', dtype=np.float32, set_initial=True)
        return particle_type

    def file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART):
        odirec = self.output_dir + 'st_{}_rs_{}_e_{}/'.format(settings.SHORE_TIME, settings.RESUS_TIME,
                                                              settings.ENSEMBLE)
        if new:
            str_format = (settings.ADVECTION_DATA, settings.SHORE_TIME, settings.RESUS_TIME, settings.STARTYEAR,
                          settings.INPUT, restart, run)
        else:
            str_format = (settings.ADVECTION_DATA, settings.SHORE_TIME, settings.RESUS_TIME, settings.STARTYEAR,
                          settings.INPUT, restart - 1, run)
        return odirec + self.prefix + '_{}_st={}_rt={}_y={}_I={}_r={}_run={}.nc'.format(*str_format)

    def beaching_kernel(particle, fieldset, time):
        """
        The beaching and resuspension kernels for beaching on the coastline follows the procedure outlined in Onink et
        al. (2021) https://doi.org/10.1088/1748-9326/abecbd
        """
        if particle.beach == 0:
            dist = fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
            if dist < fieldset.Coastal_Boundary:
                if ParcelsRandom.uniform(0, 1) > fieldset.p_beach:
                    particle.beach = 1
        # Now the part where we build in the resuspension
        elif particle.beach == 1:
            if ParcelsRandom.uniform(0, 1) > fieldset.p_resus:
                particle.beach = 0
        # Update the age of the particle
        particle.age += particle.dt

    def get_particle_behavior(self, pset: ParticleSet):
        utils.print_statement("Setting the particle behavior")
        total_behavior = pset.Kernel(utils.initial_input) + \
                         pset.Kernel(utils.floating_advection_rk4) + \
                         pset.Kernel(utils.floating_2d_brownian_motion) + \
                         pset.Kernel(utils.anti_beach_nudging) + \
                         pset.Kernel(self.beaching_kernel)
        return total_behavior
