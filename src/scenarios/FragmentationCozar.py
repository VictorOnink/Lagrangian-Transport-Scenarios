from parcels import FieldSet, ParticleSet
from parcels.kernels.TEOSseawaterdensity import PolyTEOS10_bsq
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
from advection_scenarios import advection_files
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
                           rho_plastic=var_dict['rho_plastic'], time=start_time, repeatdt=repeat_dt)
        return pset

    def _get_pclass(self):
        os.system('echo "Creating the particle class"')
        particle_type = utils.BaseParticle
        utils._add_var_particle(particle_type, 'distance', dtype=np.float32, set_initial=False)
        utils._add_var_particle(particle_type, 'density', dtype=np.float32, set_initial=False, to_write=False)
        utils._add_var_particle(particle_type, 'kinematic_viscosity', dtype=np.float32, set_initial=False,
                                to_write=False)
        utils._add_var_particle(particle_type, 'rise_velocity', dtype=np.float32, set_initial=False)
        utils._add_var_particle(particle_type, 'rho_plastic', dtype=np.float32, set_initial=True, to_write=False)
        utils._add_var_particle(particle_type, 'size', dtype=np.float32)
        return particle_type

    def _file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART):
        odirec = self.output_dir + "Cozar_Fragmentation/st_" + str(settings.SHORE_TIME) + "_rt_" + \
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
        # We are assuming that
        # Salinity and Temperature at the particle position, where salinity is converted from g/kg -> kg/kg
        Sz = fieldset.abs_salinity[time, particle.depth, particle.lat, particle.lon] / 1000
        Tz = fieldset.cons_temperature[time, particle.depth, particle.lat, particle.lon]
        # The constants A and B
        A = 1.541 + 1.998 * 10 ** -2 * Tz - 9.52 * 10 ** -5 * math.pow(Tz, 2)
        B = 7.974 - 7.561 * 10 ** -2 * Tz + 4.724 * 10 ** -4 * math.pow(Tz, 2)
        # Calculating the water dynamic viscosity
        mu_wz = 4.2844 * 10 ** -5 + math.pow(0.156 * math.pow(Tz + 64.993, 2) - 91.296, -1)
        # Calculating the sea water kinematic viscosity
        particle.kinematic_viscosity = mu_wz * (1 + A * Sz + B * math.pow(Sz, 2)) / particle.density

    def _get_rising_velocity(particle, fieldset, time):
        """
        Kernel to compute the vertical velocity (Vs) of particles due to their different sizes and densities
        This is very heavily based on the Kernel "Kooi_no_biofouling" written by Delphine Lobelle
        https://github.com/dlobelle/TOPIOS/blob/master/scripts/Kooi%2BNEMO_3D_nobiofoul.py
        """

        # ------ Profiles from MEDUSA or Kooi theoretical profiles -----
        z = particle.depth  # [m]
        # if particle.age == 0:
        #     particle.depth = fieldset.bathymetry[time, particle.depth, particle.lat, particle.lon] - 1
        kin_visc = particle.kinematic_viscosity  # kinematic viscosity[m2 s-1]
        rho_sw = particle.density  # seawater density[kg m-3]
        rise = particle.rise_velocity  # vertical velocity[m s-1]

        # ------ Constants -----
        g = 7.32e10 / (86400. ** 2.)  # gravitational acceleration (m d-2), now [s-2]

        # ------ Volumes -----
        v_pl = (4. / 3.) * math.pi * particle.size ** 3.  # volume of plastic [m3]
        theta_pl = 4. * math.pi * particle.size ** 2.  # surface area of plastic particle [m2]

        # ------ Diffusivity -----
        r_tot = particle.size                             # total radius [m]
        rho_tot = (particle.size ** 3. * particle.rho_plastic) / (particle.size) ** 3.  # total density [kg m-3]

        dn = 2. * (r_tot)  # equivalent spherical diameter [m]
        delta_rho = (rho_tot - rho_sw) / rho_sw  # normalised difference in density between total plastic+bf and seawater[-]
        dstar = ((rho_tot - rho_sw) * g * dn ** 3.) / (rho_sw * kin_visc ** 2.)  # dimensional diameter[-]

        # Getting the dimensionless settling velocity
        if dstar > 5e9:
            w = 1000.
        elif dstar < 0.05:
            w = (dstar ** 2.) * 1.71E-4
        else:
            w = 10. ** (-3.76715 + (1.92944 * math.log10(dstar)) - (0.09815 * math.log10(dstar) ** 2.) - (
                        0.00575 * math.log10(dstar) ** 3.) + (0.00056 * math.log10(dstar) ** 4.))
        # ------ Settling of particle -----

        if delta_rho > 0:  # sinks
            vs = (g * kin_visc * w * delta_rho) ** (1. / 3.)
        else:  # rises
            a_del_rho = delta_rho * -1.
            vs = -1. * (g * kin_visc * w * a_del_rho) ** (1. / 3.)  # m s-1

        z0 = z + vs * particle.dt
        # # 1.472102
        # if z0 <= 1.472102 or z0 >= fieldset.bathymetry[time, particle.depth, particle.lat, particle.lon]:
        #     # particle.depth = 1.472102
        #     vs = 0
            # particle.depth = z0

        particle.rise_velocity = math.floor(fieldset.landID[time, particle.depth, particle.lat, particle.lon])

    def _get_particle_behavior(self, pset: ParticleSet):
        os.system('echo "Setting the particle behavior"')
        base_behavior = pset.Kernel(utils._initial_input) + pset.Kernel(PolyTEOS10_bsq) + \
                        pset.Kernel(self._get_kinematic_viscosity) + \
                        pset.Kernel(utils._floating_advection_rk4) + \
                        pset.Kernel(utils._floating_2d_brownian_motion) + \
                        pset.Kernel(self._get_rising_velocity)
        total_behavior = base_behavior + pset.Kernel(utils._anti_beach_nudging) + pset.Kernel(self._beaching_kernel)
        return total_behavior
