from parcels import FieldSet, ParticleSet
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
import utils
from datetime import datetime, timedelta


class CoastalProximity(base_scenario.BaseScenario):
    """Coastal proximity scenario"""

    def __init__(self):
        """Constructor for coastal_proximity"""
        super().__init__()

    def set_prefix(self) -> str:
        """
        Set the scenario advection_prefix
        :return:
        """
        return "Prox"

    def set_var_list(self) -> list:
        """
        Set the var_list, which contains all the variables that need to be loaded during the restarts
        :return:
        """
        return ['lon', 'lat', 'weights', 'beach', 'age', 'prox']

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
        fieldset = fieldset_factory.FieldSetFactory().create_fieldset(file_dict=self.file_dict,
                                                                      border_current=True, diffusion=True,
                                                                      landID=True, distance=True, vicinity=True,
                                                                      time_step=self.dt)
        return fieldset

    def get_pset(self, fieldset: FieldSet, particle_type: utils.BaseParticle, var_dict: dict,
                 start_time: datetime, repeat_dt: timedelta):
        """
        :return:
        """
        utils.print_statement("Creating the particle set")
        pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                           lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                           age=var_dict['age'], prox=var_dict['prox'], weights=var_dict['weight'],
                           time=start_time, repeatdt=repeat_dt)
        return pset

    def get_pclass(self):
        utils.print_statement("Creating the particle class")
        particle_type = utils.BaseParticle
        utils.add_particle_variable(particle_type, 'prox')
        utils.add_particle_variable(particle_type, 'distance', dtype=np.float32, set_initial=False)
        utils.add_particle_variable(particle_type, 'weights', dtype=np.float32, set_initial=True)
        return particle_type

    def file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART):
        odirec = self.output_dir + "coastal_v_" + str(settings.VICINITY) + "_e_" + str(settings.ENSEMBLE) + "/"
        if new:
            str_format = (settings.ADVECTION_DATA, settings.VICINITY, settings.STARTYEAR, settings.INPUT, restart, run)
        else:
            str_format = (settings.ADVECTION_DATA, settings.VICINITY, settings.STARTYEAR, settings.INPUT, restart - 1,
                          run)
        return odirec + self.prefix + '_{}_v={}_y={}_I={}_r={}_run={}.nc'.format(*str_format)

    def beaching_kernel(particle, fieldset, time):
        if particle.beach == 0:
            dist = fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
            # If a particle is within 10 km of the shore
            if dist < fieldset.Coastal_Boundary:
                particle.prox += particle.dt
            else:
                particle.prox = 0.
            if particle.prox > 86400 * fieldset.vic:
                particle.beach = 1
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
