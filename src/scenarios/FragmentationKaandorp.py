from parcels import FieldSet, ParticleSet
# from parcels.kernels.TEOSseawaterdensity import PolyTEOS10_bsq
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
        self.input_dir = utils.get_input_directory(server=self.server)
        self.output_dir = utils.get_output_directory(server=self.server)
        self.repeat_dt = None
        if settings.SUBMISSION in ['simulation', 'visualization']:
            advection_scenario = advection_files.AdvectionFiles(server=self.server, stokes=self.stokes,
                                                                advection_scenario=settings.ADVECTION_DATA,
                                                                repeat_dt=self.repeat_dt)
            self.file_dict = advection_scenario.file_names
            self.field_set = self.create_fieldset()

    var_list = ['lon', 'lat', 'weights', 'beach', 'age', 'size', 'rho_plastic']

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

    def get_pset(self, fieldset: FieldSet, particle_type: utils.BaseParticle, var_dict: dict,
                 start_time: datetime, repeat_dt: timedelta):
        """
        :return:
        """
        os.system('echo "Creating the particle set"')
        if settings.RESTART == 0:
            pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                               lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                               age=var_dict['age'], weights=var_dict['weight'], size=var_dict['size'],
                               rho_plastic=var_dict['rho_plastic'],
                               rise_velocity=utils.initial_estimate_particle_rise_velocity(L=var_dict['size']),
                               time=start_time, repeatdt=repeat_dt)
        else:
            pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                               lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                               age=var_dict['age'], weights=var_dict['weight'], size=var_dict['size'],
                               rho_plastic=var_dict['rho_plastic'], rise_velocity=var_dict['rise_velocity'],
                               time=start_time, repeatdt=repeat_dt)
        return pset

    def get_pclass(self):
        os.system('echo "Creating the particle class"')
        particle_type = utils.BaseParticle
        utils.add_particle_variable(particle_type, 'distance', dtype=np.float32, set_initial=False)
        utils.add_particle_variable(particle_type, 'density', dtype=np.float32, set_initial=False, to_write=True)
        utils.add_particle_variable(particle_type, 'surface_density', dtype=np.float32, set_initial=False,
                                    to_write=True)
        utils.add_particle_variable(particle_type, 'kinematic_viscosity', dtype=np.float32, set_initial=False,
                                    to_write=False)
        utils.add_particle_variable(particle_type, 'rise_velocity', dtype=np.float32, set_initial=True)
        utils.add_particle_variable(particle_type, 'reynolds', dtype=np.float32, set_initial=False)
        utils.add_particle_variable(particle_type, 'rho_plastic', dtype=np.float32, set_initial=True, to_write=False)
        utils.add_particle_variable(particle_type, 'size', dtype=np.float32)
        utils.add_particle_variable(particle_type, 'weights', dtype=np.float32, set_initial=True)
        return particle_type

    def file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART,
                   shore_time=settings.SHORE_TIME, resus_time=settings.RESUS_TIME, ensemble=settings.ENSEMBLE,
                   advection_data=settings.ADVECTION_DATA, start_year=settings.START_YEAR, input=settings.INPUT):
        odirec = self.output_dir + "Kaandorp_Fragmentation/st_{}_rt_{}_e_{}/".format(shore_time,
                                                                                     resus_time,
                                                                                     ensemble)
        if new:
            os.system('echo "Set the output file name"')
            str_format = (advection_data, shore_time, resus_time, start_year, input, restart, run)
        else:
            os.system('echo "Set the restart file name"')
            str_format = (advection_data, shore_time, resus_time, start_year, input, restart - 1, run)
        return odirec + self.prefix + '_{}_st={}_rt={}_y={}_I={}_r={}_run={}.nc'.format(*str_format)

    def beaching_kernel(particle, fieldset, time):
        """
        The basic beaching and resuspension procedure Onink et al. (2021). However, since Hinata et al. (2017) showed
        that the resuspension  varies with
        Onink et al. (2021) = https://doi.org/10.1088/1748-9326/abecbd
        Hinata et al. (2017) = https://doi.org/10.1016/j.marpolbul.2017.05.012

        For particles on the seabed, we follow the resuspension procedure outlined in Carvajalino-Fernandez et al.
        (2020), where a particle at the sea bed gets resuspended if the estimated sea floor sea stress is greater than a
        critical threshold. Particles can get stuck on the seabed if the potential depth due to KPP or internal tide
        mixing is below the bathymetry depth

        The bottom sea stress is calculated using a quadratic drag extrapolation according to Warner et al. (2008)
        Carvajalino-Fernandez et al. (2020) = https://doi.org/10.1016/j.marpolbul.2020.111685
        Warner et al. (2008) = https://doi.org/10.1016/j.cageo.2008.02.012
        """
        # First, the beaching of particles on the coastline
        if particle.beach == 0:
            dist = fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
            if dist < fieldset.Coastal_Boundary:
                if ParcelsRandom.uniform(0, 1) > fieldset.p_beach:
                    particle.beach = 1
        # Next the resuspension of particles on the coastline
        elif particle.beach == 1:
            lambda_resus = 2.6e2 * math.fabs(particle.rise_velocity) + 7.1
            prob_resus = math.exp(-particle.dt / (lambda_resus * 86400.))
            if ParcelsRandom.uniform(0, 1) > prob_resus:
                particle.beach = 0
        # Finally, the resuspension of particles on the seabed
        elif particle.beach == 3:
            dWx = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
            dWy = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)

            bx = math.sqrt(2 * fieldset.SEABED_KH)

            # Getting the current strength at the particle position at the sea bed, and converting it to m/s
            U_bed, V_bed = fieldset.U[time, particle.depth, particle.lat, particle.lon], fieldset.V[time, particle.depth, particle.lat, particle.lon]
            U_bed, V_bed = U_bed * 1852. * 60. * math.cos(40. * math.pi / 180.), V_bed * 1852. * 60.
            U_bed, V_bed = U_bed + bx * dWx, V_bed + bx * dWy
            # Getting the bottom shear stress
            tau_bss = 0.003 * (math.pow(U_bed, 2) + math.pow(V_bed, 2))
            # if tau_bss is greater than fieldset.SEABED_CRIT, then the particle gets resuspended
            if tau_bss > fieldset.SEABED_CRIT:
                particle.beach = 0
        # Update the age of the particle
        particle.age += particle.dt

    def get_particle_behavior(self, pset: ParticleSet):
        os.system('echo "Setting the particle behavior"')
        total_behavior = pset.Kernel(utils.PolyTEOS10_bsq) + \
                        pset.Kernel(utils.get_kinematic_viscosity) + \
                        pset.Kernel(utils.get_reynolds_number) + \
                        pset.Kernel(utils.floating_AdvectionRK4DiffusionEM_stokes_depth) + \
                        pset.Kernel(utils.anti_beach_nudging) + \
                        pset.Kernel(utils.get_rising_velocity) + \
                        pset.Kernel(utils.KPP_TIDAL_mixing) + \
                        pset.Kernel(self.beaching_kernel)
        return total_behavior

    def run(self):
        os.system('echo "Creating the particle set"')
        pset = self.get_pset(fieldset=self.field_set, particle_type=self.particle,
                             var_dict=self.get_var_dict(), start_time=utils.get_start_end_time(time='start'),
                             repeat_dt=self.repeat_dt)
        pfile = pset.ParticleFile(name=self.file_names(new=True),
                                  outputdt=settings.OUTPUT_TIME_STEP)
        os.system('echo "Setting the random seed"')
        utils.set_random_seed(seed=settings.SEED)
        os.system('echo "Defining the particle behavior"')
        behavior_kernel = self.get_particle_behavior(pset=pset)
        os.system('echo "The actual execution of the run"')
        # pset.execute(behavior_kernel,
        #              runtime=timedelta(days=get_start_end_time(time='length')),
        #              dt=settings.TIME_STEP,
        #              # recovery={ErrorCode.ErrorOutOfBounds: delete_particle},
        #              output_file=pfile
        #              )
        # pfile.export()
        os.system('echo "Run completed"')
