from parcels import FieldSet, ParticleSet
from parcels.kernels.TEOSseawaterdensity import PolyTEOS10_bsq
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
import utils as utils
from datetime import datetime, timedelta
import os
from parcels import rng as random
import math


class FragmentationCozar(base_scenario.BaseScenario):
    """Fragmentation scenario based on the Cozar fragmentation model. Beaching based on the stochastic scenario"""

    def __init__(self, server, stokes):
        """Constructor for FragmentationCozar"""
        super().__init__(server, stokes)
        self.prefix = "Frag_cozar"
        self.input_dir = utils._get_input_directory(server=self.server)
        self.output_dir = utils._get_output_directory(server=self.server)
        self.repeat_dt = None

    var_list = ['lon', 'lat', 'beach', 'age', 'size', 'density']

    def create_fieldset(self) -> FieldSet:
        os.system('echo "Creating the fieldset"')
        fieldset = fieldset_factory.FieldSetFactory().create_fieldset(file_dict=self.file_dict, stokes=self.stokes,
                                                                      border_current=True, diffusion=True, landID=True,
                                                                      distance=True, salinity=True, temperature=True,
                                                                      bathymetry=True, beach_timescale=True,
                                                                      resus_timescale=True
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
                           time=start_time, repeatdt=repeat_dt)
        return pset

    def _get_pclass(self):
        os.system('echo "Creating the particle class"')
        particle_type = utils.BaseParticle
        utils._add_var_particle(particle_type, 'distance', dtype=np.float32, set_initial=False)
        utils._add_var_particle(particle_type, 'density', dtype=np.float32, set_initial=False, to_write=False)
        utils._add_var_particle(particle_type, 'kinematic_viscosity', dtype=np.float32, set_initial=False)
        utils._add_var_particle(particle_type, 'size', dtype=np.float32)
        return particle_type

    def _file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART):
        odirec = self.output_dir + "Cozar_Fragmentation/st_" + str(settings.SHORE_TIME) + "_rs_" + \
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
                if random.random() > fieldset.p_beach:
                    particle.beach = 1
        # Now the part where we build in the resuspension
        elif particle.beach == 1:
            if random.random() > fieldset.p_resus:
                particle.beach = 0
        # Update the age of the particle
        particle.age += particle.dt


    def _get_kinematic_viscosity(particle, fieldset, time):
        # Using equations 25 - 29 from Kooi et al. 2017
        # Salinity and Temperature at the particle position, where salinity is converted from g/kg -> kg/kg
        Sz = fieldset.abs_salinity[time, particle.depth, particle.lat, particle.lon]/1000
        Tz = fieldset.cons_temperature[time, particle.depth, particle.lat, particle.lon]
        # The constants A and B
        A = 1.541 + 1.998 * 10**-2 * Tz - 9.52 * 10**-5 * math.pow(Tz, 2)
        B = 7.974 - 7.561 * 10**-2 * Tz + 4.724 * 10**-4 * math.pow(Tz, 2)
        # Calculating the water dynamic viscosity
        mu_wz = 4.2844 * 10**-5 + math.pow(0.156 * math.pow(Tz + 64.993, 2) - 91.296, -1)
        # Calculating the sea water kinematic viscosity
        particle.kinematic_viscosity = mu_wz*(1 + A * Sz + B * math.pow(Sz, 2)) / particle.density




    def _get_particle_behavior(self, pset: ParticleSet):
        os.system('echo "Setting the particle behavior"')
        base_behavior = pset.Kernel(utils._initial_input) + pset.Kernel(PolyTEOS10_bsq) + \
                        pset.Kernel(self._get_kinematic_viscosity) + \
                        pset.Kernel(utils._floating_advection_rk4) + \
                        pset.Kernel(utils._floating_2d_brownian_motion)
        total_behavior = base_behavior + pset.Kernel(utils._anti_beach_nudging) + pset.Kernel(self._beaching_kernel)
        return total_behavior
