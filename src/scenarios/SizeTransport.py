from parcels import FieldSet, ParticleSet
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
import utils
from datetime import datetime, timedelta
from parcels import ParcelsRandom


class SizeTransport(base_scenario.BaseScenario):
    """Fragmentation scenario based on the Kaandorp fragmentation model. Beaching based on the stochastic scenario"""

    def __init__(self):
        """Constructor for SizeTransport"""
        super().__init__()

    def set_prefix(self) -> str:
        """
        Set the scenario advection_prefix
        :return:
        """
        return "Size_Transport"

    def set_var_list(self) -> list:
        """
        Set the var_list, which contains all the variables that need to be loaded during the restarts
        :return:
        """
        return ['lon', 'lat', 'beach', 'age', 'distance2coast', 'z']

    def set_time_steps(self) -> tuple:
        """
        Set the integration, output and repeat timesteps
        :return: self.dt, self.output_time_step, self.repeat_dt
        """
        dt = timedelta(minutes=0.5 * settings.BACKWARD_MULT)
        output_time_step = timedelta(hours=12)
        repeat_dt = None
        return dt, output_time_step, repeat_dt

    def create_fieldset(self) -> FieldSet:
        utils.print_statement("Creating the fieldset")
        fieldset = fieldset_factory.FieldSetFactory().UV_interpolation(file_dict=self.file_dict,
                                                                      stokes_depth=True, border_current=True,
                                                                      diffusion=True,
                                                                      distance=True, salinity=True, temperature=True,
                                                                      bathymetry=True, beach_timescale=True,
                                                                      fixed_resus=settings.FIXED_RESUS,
                                                                      resus_timescale=True, MLD=True,
                                                                      physics_constants=True,
                                                                      wind=True, TIDAL_mixing=True,
                                                                      time_step=self.dt
                                                                      )
        return fieldset

    def get_pset(self, fieldset: FieldSet, particle_type: utils.BaseParticle, var_dict: dict,
                 start_time: datetime, repeat_dt: timedelta):
        """
        :return:
        """
        utils.print_statement("Creating the particle set")
        if settings.RESTART == 0:
            pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                               lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                               age=var_dict['age'], time=start_time, repeatdt=repeat_dt)
        else:
            pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                               lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                               age=var_dict['age'], time=start_time,
                               distance2coast=var_dict['distance2coast'], depth=var_dict['z'],
                               repeatdt=repeat_dt)
        return pset

    def get_pclass(self):
        utils.print_statement("Creating the particle class")
        particle_type = utils.BaseParticle
        if settings.RESTART == 0:
            utils.add_particle_variable(particle_type, 'distance2coast', dtype=np.float32, set_initial=False,
                                        to_write=True)
        else:
            utils.add_particle_variable(particle_type, 'distance2coast', dtype=np.float32, set_initial=True)
        utils.add_particle_variable(particle_type, 'prev_lon', dtype=np.float32, set_initial=True, to_write=False,
                                    other_name='lon')
        utils.add_particle_variable(particle_type, 'prev_lat', dtype=np.float32, set_initial=True, to_write=False,
                                    other_name='lat')
        utils.add_particle_variable(particle_type, 'prev_depth', dtype=np.float32, set_initial=True, to_write=False,
                                    other_name='depth')
        utils.add_particle_variable(particle_type, 'potential', dtype=np.float32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'density', dtype=np.float32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'surface_density', dtype=np.float32, set_initial=False,
                                    to_write=False)
        utils.add_particle_variable(particle_type, 'kinematic_viscosity', dtype=np.float32, set_initial=False,
                                    to_write=False)
        utils.add_particle_variable(particle_type, 'rise_velocity', dtype=np.float32, set_initial=True,
                                    other_value=utils.initial_estimate_particle_rise_velocity(), to_write=False)
        utils.add_particle_variable(particle_type, 'reynolds', dtype=np.float32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'rho_plastic', dtype=np.float32, set_initial=True, to_write=False,
                                    other_value=settings.INIT_DENSITY)
        utils.add_particle_variable(particle_type, 'size', dtype=np.float32, to_write=False,
                                    other_value=settings.INIT_SIZE)
        return particle_type

    def file_names(self, new: bool = False, advection_data: str = settings.ADVECTION_DATA,
                   shore_time: int = settings.SHORE_TIME, init_size: float = settings.INIT_SIZE,
                   init_density: int = settings.INIT_DENSITY, start_year: int = settings.STARTYEAR,
                   input: str = settings.INPUT, run: int = settings.RUN, restart: int = settings.RESTART,
                   seabed_crit: float = settings.SEABED_CRIT, fixed_resus: bool = settings.FIXED_RESUS,
                   resus_time: str = settings.RESUS_TIME):
        odirec = self.output_dir + "SizeTransport/size_{:.1E}/".format(init_size)
        if not fixed_resus:
            resus_time = utils.get_resuspension_timescale(L=init_size, rho_p=init_density)
        if new:
            str_format = (advection_data, shore_time, resus_time, init_size, init_density, seabed_crit, start_year,
                          input, restart, run)
        else:
            str_format = (advection_data, shore_time, resus_time, init_size, init_density, seabed_crit, start_year,
                          input, restart - 1, run)
        return odirec + self.prefix + '_{}_st={}_rt={:.6f}_size={:.1E}_rho={}_taubss={:.2E}_y={}_I={}_r={}_run={}.nc'.format(*str_format)

    def beaching_kernel(particle, fieldset, time):
        """
        The beaching and resuspension kernels for beaching on the coastline follows the procedure outlined in Onink et
        al. (2021).
        Onink et al. (2021) = https://doi.org/10.1088/1748-9326/abecbd
        """
        # First, the beaching of particles on the coastline
        if particle.beach == 0:
            particle.distance2coast = fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
            if particle.distance2coast < fieldset.Coastal_Boundary:
                if ParcelsRandom.uniform(0, 1) > fieldset.p_beach:
                    particle.beach = 1
        # Next the resuspension of particles on the coastline
        elif particle.beach == 1:
            if ParcelsRandom.uniform(0, 1) > fieldset.p_resus:
                particle.beach = 0
        # Update the age of the particle
        particle.age += particle.dt

    def get_particle_behavior(self, pset: ParticleSet):
        utils.print_statement("Setting the particle behavior")
        base_behavior = pset.Kernel(utils.PolyTEOS10_bsq) + \
                        pset.Kernel(utils.get_kinematic_viscosity) + \
                        pset.Kernel(utils.get_reynolds_number) + \
                        pset.Kernel(utils.floating_AdvectionRK4DiffusionEM_stokes_depth) + \
                        pset.Kernel(utils.anti_beach_nudging) + \
                        pset.Kernel(utils.get_rising_velocity) + \
                        pset.Kernel(utils.KPP_TIDAL_mixing) + \
                        pset.Kernel(utils.vertical_reflecting_boundary) + \
                        pset.Kernel(self.beaching_kernel)
        return base_behavior
