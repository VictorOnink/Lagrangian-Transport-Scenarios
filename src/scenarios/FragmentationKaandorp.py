from parcels import FieldSet, ParticleSet
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
from advection_scenarios import advection_files
import utils as utils
from datetime import datetime, timedelta
import os
from parcels import ParcelsRandom
import math


class FragmentationKaandorp(base_scenario.BaseScenario):
    """Fragmentation scenario based on the Kaandorp fragmentation model. Beaching based on the stochastic scenario"""

    def __init__(self, server, stokes):
        """Constructor for FragmentationKaandorp"""
        super().__init__(server, stokes)
        self.prefix = "Frag_Kaandorp"
        self.input_dir = utils._get_input_directory(server=self.server)
        self.output_dir = utils._get_output_directory(server=self.server)
        self.repeat_dt = None
        if settings.SUBMISSION == 'simulation':
            advection_scenario = advection_files.AdvectionFiles(server=self.server, stokes=self.stokes,
                                                                advection_scenario=settings.ADVECTION_DATA,
                                                                repeat_dt=self.repeat_dt)
            self.file_dict = advection_scenario.file_names
            self.field_set = self.create_fieldset()

    var_list = ['lon', 'lat', 'beach', 'age', 'size', 'rho_plastic']

    def create_fieldset(self) -> FieldSet:
        os.system('echo "Creating the fieldset"')
        fieldset = fieldset_factory.FieldSetFactory().create_fieldset(file_dict=self.file_dict, stokes=self.stokes,
                                                                      stokes_depth=True,
                                                                      border_current=True, diffusion=True, landID=True,
                                                                      distance=True, salinity=True, temperature=True,
                                                                      bathymetry=True, beach_timescale=True,
                                                                      resus_timescale=True, MLD=True, KPP_mixing=True,
                                                                      wind=True
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
                           age=var_dict['age'], weights=var_dict['weight'], size=var_dict['size'],
                           rho_plastic=var_dict['rho_plastic'], time=start_time, repeatdt=repeat_dt)
        return pset

    def _get_pclass(self):
        os.system('echo "Creating the particle class"')
        particle_type = utils.BaseParticle
        utils._add_var_particle(particle_type, 'distance', dtype=np.float32, set_initial=False)
        utils._add_var_particle(particle_type, 'density', dtype=np.float32, set_initial=False, to_write=False)
        utils._add_var_particle(particle_type, 'surface_density', dtype=np.float32, set_initial=False, to_write=False)
        utils._add_var_particle(particle_type, 'kinematic_viscosity', dtype=np.float32, set_initial=False,
                                to_write=False)
        utils._add_var_particle(particle_type, 'rise_velocity', dtype=np.float32, set_initial=False)
        utils._add_var_particle(particle_type, 'reynolds', dtype=np.float32, set_initial=False)
        utils._add_var_particle(particle_type, 'rho_plastic', dtype=np.float32, set_initial=True, to_write=False)
        utils._add_var_particle(particle_type, 'size', dtype=np.float32)
        return particle_type

    def _file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART):
        odirec = self.output_dir + "Kaandorp_Fragmentation/st_" + str(settings.SHORE_TIME) + "_rt_" + \
                 str(settings.RESUS_TIME) + "_e_" + str(settings.ENSEMBLE) + "/"
        if new == True:
            os.system('echo "Set the output file name"')
            return odirec + self.prefix + '_{}'.format(settings.ADVECTION_DATA) + "_st=" + str(settings.SHORE_TIME) + \
                   "_rt=" + str(settings.RESUS_TIME) + "_y=" + str(settings.START_YEAR) + "_I=" + str(settings.INPUT) \
                   + "_r=" + str(restart) + "_run=" + str(run) + ".nc"
        else:
            os.system('echo "Set the restart file name"')
            return odirec + self.prefix + '_{}'.format(settings.ADVECTION_DATA) + "_st=" + str(settings.SHORE_TIME) + \
                   "_rt=" + str(settings.RESUS_TIME) + "_y=" + str(settings.START_YEAR) + "_I=" + str(settings.INPUT) \
                   + "_r=" + str(restart - 1) + "_run=" + str(run) + ".nc"

    def _beaching_kernel(particle, fieldset, time):
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

    def _get_rising_velocity(particle, fieldset, time):
        Re = particle.reynolds  # Reynolds number
        rho_sw = particle.density  # sea water density (kg m^-3)
        rho_p = particle.rho_plastic  # plastic particle density (kg m^-3)
        L = particle.size  # particle size (m)
        g = 9.81  # gravitational acceleration (m s^-2)
        # Getting the equation according to Poulain et al (2019), equation 5 in the supplementary materials
        left = 240 / (math.pi * Re) * (1 + 0.138 * Re ** 0.792)
        right = 2. / 15. * L * (1. - rho_p/rho_sw) * g
        # Calculate the rise velocity
        particle.rise_velocity = - 1 * math.sqrt(right / left)

    def _get_reynolds_number(particle, fieldset, time):
        kin_visc = particle.kinematic_viscosity
        L = particle.size
        if particle.age == 0:
            # An initial starting value for the Reynolds number, used to calculate the first rising velocity
            particle.reynolds = 2.0
        else:
            w_b = math.fabs(particle.rise_velocity)
            particle.reynolds = L * w_b / kin_visc

    def _get_particle_behavior(self, pset: ParticleSet):
        os.system('echo "Setting the particle behavior"')
        base_behavior = pset.Kernel(utils._initial_input) + pset.Kernel(utils.PolyTEOS10_bsq) + \
                        pset.Kernel(utils._get_kinematic_viscosity) + \
                        pset.Kernel(self._get_reynolds_number) + \
                        pset.Kernel(self._get_rising_velocity) + \
                        pset.Kernel(utils._floating_AdvectionRK4DiffusionEM_stokes_depth) + \
                        pset.Kernel(utils.KPP_wind_mixing)
        total_behavior = base_behavior + pset.Kernel(utils._anti_beach_nudging) + pset.Kernel(self._beaching_kernel)
        return total_behavior
