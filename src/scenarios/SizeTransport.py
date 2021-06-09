from parcels import FieldSet, ParticleSet
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
from advection_scenarios import advection_files
import utils
from datetime import datetime, timedelta
import os
from parcels import ParcelsRandom
import math


class SizeTransport(base_scenario.BaseScenario):
    """Fragmentation scenario based on the Kaandorp fragmentation model. Beaching based on the stochastic scenario"""

    def __init__(self, server, stokes):
        """Constructor for FragmentationKaandorp"""
        super().__init__(server, stokes)
        self.prefix = "Size_Transport"
        self.input_dir = utils.get_input_directory(server=self.server)
        self.output_dir = utils.get_output_directory(server=self.server)
        self.repeat_dt = None
        if settings.SUBMISSION in ['simulation', 'visualization']:
            advection_scenario = advection_files.AdvectionFiles(server=self.server, stokes=self.stokes,
                                                                advection_scenario=settings.ADVECTION_DATA,
                                                                repeat_dt=self.repeat_dt)
            self.file_dict = advection_scenario.file_names
            if settings.SUBMISSION in ['simulation']:
                self.field_set = self.create_fieldset()

    var_list = ['lon', 'lat', 'beach', 'age', 'distance_horizontal', 'distance_vertical', 'z']

    def create_fieldset(self) -> FieldSet:
        os.system('echo "Creating the fieldset"')
        fieldset = fieldset_factory.FieldSetFactory().create_fieldset(file_dict=self.file_dict, stokes=self.stokes,
                                                                      stokes_depth=True,
                                                                      border_current=True, diffusion=True, landID=True,
                                                                      distance=True, salinity=True, temperature=True,
                                                                      bathymetry=True, beach_timescale=True,
                                                                      resus_timescale=True, MLD=True, KPP_mixing=True,
                                                                      wind=True, TIDAL_mixing=True,
                                                                      seabed_resuspension=True)
        return fieldset

    def _get_pset(self, fieldset: FieldSet, particle_type: utils.BaseParticle, var_dict: dict,
                  start_time: datetime, repeat_dt: timedelta):
        """
        :return:
        """
        os.system('echo "Creating the particle set"')
        if settings.RESTART == 0:
            # pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
            #                    lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
            #                    age=var_dict['age'], time=start_time, repeatdt=repeat_dt)
            pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                               lon=var_dict['lon'][:10], lat=var_dict['lat'][:10], beach=3 * np.ones(10),
                               age=var_dict['age'][:10], time=start_time, repeatdt=repeat_dt)

        else:
            pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                               lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                               age=var_dict['age'], time=start_time,
                               distance_horizontal=var_dict['distance_horizontal'],
                               distance_vertical=var_dict['distance_vertical'], depth=var_dict['z'],
                               repeatdt=repeat_dt)
        return pset

    def _get_pclass(self):
        os.system('echo "Creating the particle class"')
        particle_type = utils.BaseParticle
        if settings.RESTART == 0:
            utils.add_particle_variable(particle_type, 'distance_horizontal', dtype=np.float32, set_initial=False,
                                        to_write=True)
            utils.add_particle_variable(particle_type, 'distance_vertical', dtype=np.float32, set_initial=False,
                                        to_write=True)
        else:
            utils.add_particle_variable(particle_type, 'distance_horizontal', dtype=np.float32, set_initial=True)
            utils.add_particle_variable(particle_type, 'distance_vertical', dtype=np.float32, set_initial=True)
        utils.add_particle_variable(particle_type, 'prev_lon', dtype=np.float32, set_initial=True, to_write=False,
                                    other_name='lon')
        utils.add_particle_variable(particle_type, 'prev_lat', dtype=np.float32, set_initial=True, to_write=False,
                                    other_name='lat')
        utils.add_particle_variable(particle_type, 'prev_depth', dtype=np.float32, set_initial=True, to_write=False,
                                    other_name='depth')
        utils.add_particle_variable(particle_type, 'density', dtype=np.float32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'surface_density', dtype=np.float32, set_initial=False,
                                    to_write=False)
        utils.add_particle_variable(particle_type, 'kinematic_viscosity', dtype=np.float32, set_initial=False,
                                    to_write=False)
        utils.add_particle_variable(particle_type, 'rise_velocity', dtype=np.float32, set_initial=True,
                                    other_value=utils.initial_estimate_particle_rise_velocity())
        utils.add_particle_variable(particle_type, 'reynolds', dtype=np.float32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'rho_plastic', dtype=np.float32, set_initial=True, to_write=False,
                                    other_value=settings.INIT_DENSITY)
        utils.add_particle_variable(particle_type, 'size', dtype=np.float32, to_write=False,
                                    other_value=settings.INIT_SIZE)
        return particle_type

    def _file_names(self, new: bool = False, advection_data: str = settings.ADVECTION_DATA,
                    shore_time: int = settings.SHORE_TIME, init_size: float = settings.INIT_SIZE,
                    init_density: int = settings.INIT_DENSITY, start_year: int = settings.START_YEAR,
                    input: str = settings.INPUT, run: int = settings.RUN, restart: int = settings.RESTART,
                    seabed_crit: float = settings.SEABED_CRIT):
        odirec = self.output_dir + "SizeTransport/size_{:.1E}/".format(settings.INIT_SIZE)
        if new:
            os.system('echo "Set the output file name"')
            str_format = (
                advection_data, shore_time, utils.get_resuspension_timescale(L=init_size), init_size, init_density,
                seabed_crit, start_year, input, restart, run)
        else:
            os.system('echo "Set the restart file name"')
            str_format = (
                advection_data, shore_time, utils.get_resuspension_timescale(L=init_size), init_size, init_density,
                seabed_crit, start_year, input, restart - 1, run)
        return odirec + self.prefix + '_{}_st={}_rt={:.6f}_size={:.1E}_rho={}_taubss={:.2E}_y={}_I={}_r={}_run={}.nc'.format(*str_format)

    def _beaching_kernel(particle, fieldset, time):
        """
        The beaching and resuspension kernels for beaching on the coastline follows the procedure outlined in Onink et
        al. (2021).
        Onink et al. (2021) = https://doi.org/10.1088/1748-9326/abecbd

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
            if ParcelsRandom.uniform(0, 1) > fieldset.p_resus:
                particle.beach = 0
        # Finally, the resuspension of particles on the seabed
        elif particle.beach == 3:
            # Getting the current strength at the particle position at the sea bed, and converting it to m/s
            bath = fieldset.bathymetry[time, particle.depth, particle.lat, particle.lon]
            # U_bed, V_bed = fieldset.U[time, particle.depth, particle.lat, particle.lon], fieldset.V[time, particle.depth, particle.lat, particle.lon]
            U_bed, V_bed = fieldset.U[time, bath, particle.lat, particle.lon], fieldset.V[time, bath, particle.lat, particle.lon]
            id = particle.id
            U_bed, V_bed = U_bed * 1852. * 60. * math.cos(40. * math.pi / 180.), V_bed * 1852. * 60.
            # Getting the bottom shear stress
            tau_bss = 0.003 * (math.pow(U_bed, 2) + math.pow(V_bed, 2))
            print(id)
            print(U_bed)
            print(V_bed)
            print(tau_bss)
            # if tau_bss is greater than fieldset.SEABED_CRIT, then the particle gets resuspended
            if tau_bss > fieldset.SEABED_CRIT:
                particle.beach = 0
        # Update the age of the particle
        particle.age += particle.dt

    def _get_rising_velocity(particle, fieldset, time):
        rho_sw = particle.density  # sea water density (kg m^-3)
        rho_p = particle.rho_plastic  # plastic particle density (kg m^-3)
        left = (1. - rho_p / rho_sw) * 8. / 3. * particle.size * fieldset.G
        right = 24. / particle.reynolds + 5. / math.sqrt(particle.reynolds) + 2. / 5.
        particle.rise_velocity = - 1 * math.sqrt(left / right)

    def _get_reynolds_number(particle, fieldset, time):
        w_b = math.fabs(particle.rise_velocity)
        particle.reynolds = particle.size * w_b / particle.kinematic_viscosity

    def _TotalDistance(particle, fieldset, time):
        """
        Calculating the cumulative distance travelled by the particle in vertical and horizontal directions
        """
        # Calculate the distance in latitudinal direction (using 1.11e2 kilometer per degree latitude)
        lat_dist = (particle.lat - particle.prev_lat) * 1.11e2
        # Calculate the distance in longitudinal direction, using cosine(latitude) - spherical earth
        lon_dist = (particle.lon - particle.prev_lon) * 1.11e2 * math.cos(particle.lat * math.pi / 180)
        # Calculate the total Euclidean distance travelled by the particle
        particle.distance_horizontal += math.sqrt(math.pow(lon_dist, 2) + math.pow(lat_dist, 2))
        particle.distance_vertical += math.fabs(particle.depth - particle.prev_depth)

        particle.prev_lon = particle.lon  # Set the stored values for next iteration.
        particle.prev_lat = particle.lat
        particle.prev_depth = particle.depth

    def _get_particle_behavior(self, pset: ParticleSet):
        os.system('echo "Setting the particle behavior"')
        base_behavior = pset.Kernel(utils.PolyTEOS10_bsq) + \
                        pset.Kernel(utils._get_kinematic_viscosity) + \
                        pset.Kernel(self._get_reynolds_number) + \
                        pset.Kernel(utils._floating_AdvectionRK4DiffusionEM_stokes_depth) + \
                        pset.Kernel(utils._anti_beach_nudging) + \
                        pset.Kernel(self._get_rising_velocity) + \
                        pset.Kernel(utils.KPP_wind_mixing) + \
                        pset.Kernel(self._TotalDistance) + \
                        pset.Kernel(utils.internal_tide_mixing) + \
                        pset.Kernel(self._beaching_kernel)
        return base_behavior
