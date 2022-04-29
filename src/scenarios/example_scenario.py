from parcels import FieldSet, ParticleSet
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
import utils as utils
from datetime import datetime, timedelta
from parcels import ParcelsRandom


class ExampleScenario(base_scenario.BaseScenario):
    """Example scenario"""

    def __init__(self):
        """Constructor for ExampleScenario"""
        super().__init__()

    def set_prefix(self) -> str:
        """
        Set the scenario advection_prefix
        :return:
        """
        return "ExampleScenario"

    def set_var_list(self) -> list:
        """
        Set the var_list, which contains all the variables that need to be loaded during the restarts
        :return:
        """
        return ['lon', 'lat', 'beach']

    def set_time_steps(self) -> tuple:
        """
        Set the integration, output and repeat timesteps
        :return: self.dt, self.output_time_step, self.repeat_dt
        """
        dt = timedelta(minutes=5 * settings.BACKWARD_MULT)
        output_time_step = timedelta(hours=12)
        repeat_dt = timedelta(days=1)
        return dt, output_time_step, repeat_dt

    def create_fieldset(self) -> FieldSet:
        utils.print_statement("Creating the fieldset")
        fieldset = fieldset_factory.FieldSetFactory().create_fieldset(file_dict=self.file_dict)
        return fieldset

    def get_pset(self, fieldset: FieldSet, particle_type: utils.BaseParticle, var_dict: dict,
                 start_time: datetime, repeat_dt: timedelta):
        """
        :return:
        """
        utils.print_statement("Creating the particle set")
        pset = ParticleSet(fieldset=fieldset, pclass=particle_type, lon=var_dict['lon'], lat=var_dict['lat'],
                           beach=var_dict['beach'], time=start_time, repeatdt=repeat_dt)
        return pset

    def get_pclass(self):
        utils.print_statement("Creating the particle class")
        particle_type = utils.BaseParticle
        utils.add_particle_variable(particle_type, 'distance', dtype=np.float32, set_initial=False, to_write=False)
        return particle_type

    def file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART,
                   backward: bool = settings.BACKWARD):
        odirec = self.output_dir + 'ExampleScenario/'
        if new:
            str_format = (settings.ADVECTION_DATA, backward, settings.SHORE_TIME, settings.RESUS_TIME,
                          settings.STARTYEAR, restart, run)
        else:
            str_format = (settings.ADVECTION_DATA, backward, settings.SHORE_TIME, settings.RESUS_TIME,
                          settings.STARTYEAR, restart - 1, run)
        return odirec + self.prefix + '_{}_{}_st={}_rt={}_y={}_r={}_run={}.nc'.format(*str_format)

    def beaching_kernel(particle, fieldset, time):
        pass

    def get_particle_behavior(self, pset: ParticleSet):
        utils.print_statement("Setting the particle behavior")
        total_behavior = pset.Kernel(utils.floating_advection_rk4)
        return total_behavior
