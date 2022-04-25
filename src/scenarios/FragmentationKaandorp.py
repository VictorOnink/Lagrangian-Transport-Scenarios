from parcels import FieldSet, ParticleSet, ErrorCode
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
import utils
from datetime import datetime, timedelta
from parcels import ParcelsRandom
import math


class FragmentationKaandorp(base_scenario.BaseScenario):
    """Fragmentation scenario based on the Kaandorp fragmentation model. Beaching based on the stochastic scenario"""

    def __init__(self, server, stokes):
        """Constructor for FragmentationKaandorp"""
        super().__init__(server, stokes)

    def set_prefix(self) -> str:
        """
        Set the scenario prefix
        :return:
        """
        return "Frag_Kaandorp"

    def set_var_list(self) -> list:
        """
        Set the var_list, which contains all the variables that need to be loaded during the restarts
        :return:
        """
        return ['lon', 'lat', 'beach', 'age', 'size', 'rho_plastic', 'parent', 'rise_velocity', 'size_class']

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
        fieldset = fieldset_factory.FieldSetFactory().create_fieldset(file_dict=self.file_dict, stokes=self.stokes,
                                                                      stokes_depth=True,
                                                                      border_current=True, diffusion=True, landID=True,
                                                                      distance=True, salinity=True, temperature=True,
                                                                      bathymetry=True, beach_timescale=True,
                                                                      resus_timescale=True, MLD=True,
                                                                      physics_constants=True, wind=True,
                                                                      TIDAL_mixing=True, fragmentation_timescale=True,
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
            step = 100
            rise_velocity = utils.initial_estimate_particle_rise_velocity(L=var_dict['size'][::step])
            pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                               lon=var_dict['lon'][::step], lat=var_dict['lat'][::step], beach=var_dict['beach'][::step],
                               age=var_dict['age'][::step], size=var_dict['size'][::step],
                               rho_plastic=var_dict['rho_plastic'][::step], parent=range(len(var_dict['lon'][::step])),
                               rise_velocity=rise_velocity, size_class=np.zeros(rise_velocity.shape, dtype=np.float32),
                               prob_resus=utils.resuspension_probability(w_rise=rise_velocity, time_step=self.dt),
                               time=start_time, repeatdt=repeat_dt)
        else:
            pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                               lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'], parent=var_dict['parent'],
                               age=var_dict['age'], size=var_dict['size'], size_class=var_dict['size_class'],
                               rho_plastic=var_dict['rho_plastic'], rise_velocity=var_dict['rise_velocity'],
                               prob_resus=utils.resuspension_probability(w_rise=var_dict['rise_velocity'], time_step=self.dt),
                               time=start_time, repeatdt=repeat_dt)
        return pset

    def get_pclass(self):
        utils.print_statement("Creating the particle class")
        particle_type = utils.BaseParticle
        utils.add_particle_variable(particle_type, 'distance', dtype=np.float32, set_initial=False)
        utils.add_particle_variable(particle_type, 'density', dtype=np.float32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'surface_density', dtype=np.float32, set_initial=False,
                                    to_write=False)
        utils.add_particle_variable(particle_type, 'kinematic_viscosity', dtype=np.float32, set_initial=False,
                                    to_write=False)
        utils.add_particle_variable(particle_type, 'rise_velocity', dtype=np.float32, set_initial=True)
        utils.add_particle_variable(particle_type, 'reynolds', dtype=np.float32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'rho_plastic', dtype=np.float32, set_initial=True, to_write=False)
        utils.add_particle_variable(particle_type, 'size', dtype=np.float32)
        utils.add_particle_variable(particle_type, 'to_split', dtype=np.int32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'to_delete', dtype=np.int32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'parent', dtype=np.int32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'prob_resus', dtype=np.int32, set_initial=True, to_write=False)
        utils.add_particle_variable(particle_type, 'size_class', dtype=np.int32, set_initial=True, to_write=True)
        return particle_type

    def file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART,
                   shore_time=settings.SHORE_TIME, ensemble=settings.ENSEMBLE, advection_data=settings.ADVECTION_DATA,
                   start_year=settings.STARTYEAR, input=settings.INPUT,
                   p_frag=settings.P_FRAG, dn=settings.DN, size_class_number=settings.SIZE_CLASS_NUMBER,
                   lambda_frag=settings.LAMBDA_FRAG, density=settings.INIT_DENSITY):
        odirec = self.output_dir + "Kaandorp_Fragmentation/st_{}_e_{}/".format(shore_time, ensemble)
        if new:
            str_format = (advection_data, shore_time, p_frag, lambda_frag, dn, size_class_number, density,
                          start_year, input, restart, run)
        else:
            str_format = (advection_data, shore_time, p_frag, lambda_frag, dn, size_class_number, density,
                          start_year, input, restart - 1, run)
        return odirec + self.prefix + '_{}_st={}_pfrag={}_lambdafrag={}_dn={}_sizeclasses={}_rho={}_y={}_I={}_r={}_run={}.nc'.format(*str_format)

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
            if ParcelsRandom.uniform(0, 1) > particle.prob_resus:
                particle.beach = 0
        # Finally, the resuspension of particles on the seabed
        elif particle.beach == 3:
            dWx = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
            dWy = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)

            bx = math.sqrt(2 * fieldset.SEABED_KH)

            # Getting the current strength at the particle position at the sea bed, and converting it to m/s
            U_bed = fieldset.U[time, particle.depth, particle.lat, particle.lon]
            V_bed = fieldset.V[time, particle.depth, particle.lat, particle.lon]
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
        utils.print_statement("Setting the particle behavior")
        total_behavior = pset.Kernel(utils.PolyTEOS10_bsq) + \
                         pset.Kernel(utils.get_kinematic_viscosity) + \
                         pset.Kernel(utils.get_reynolds_number) + \
                         pset.Kernel(utils.floating_AdvectionRK4DiffusionEM_stokes_depth) + \
                         pset.Kernel(utils.anti_beach_nudging) + \
                         pset.Kernel(utils.get_rising_velocity) + \
                         pset.Kernel(utils.KPP_TIDAL_mixing) + \
                         pset.Kernel(self.beaching_kernel) + \
                         pset.Kernel(self.fragmentation_kernel)
        return total_behavior

    def fragmentation_kernel(particle, fieldset, time):
        if particle.to_delete == 1:
            particle.delete()
        else:
            if particle.beach == 1:
                if ParcelsRandom.uniform(0, 1) > fieldset.p_frag:
                    particle.to_split = 1

    def particle_splitter(self, fieldset, pset, size_limit):
        for particle in pset:
            if particle.to_split == 1:
                # First, we set the split condition statement back to 0
                particle.to_split = 0
                # Then, we calculate the new size of the particle
                parent_size = particle.size
                parent_size_class = particle.size_class
                particle.size = parent_size * (1 - settings.P_FRAG)
                # If the particle is in the smallest size class, then mark the particle for deletion since it is then
                # too small for us to follow further fragmentation
                if parent_size_class == (settings.SIZE_CLASS_NUMBER - 1):
                    particle.to_delete = 1
                # Otherwise, if the particle is not in the smallest size class then we can figure out how many fragments
                # are created in smaller size classes
                else:
                    remaining_classes = settings.SIZE_CLASS_NUMBER - parent_size_class - 1
                    for k in range(0, remaining_classes):
                        new_particle_size = parent_size * settings.P_FRAG ** (k + 1)
                        particle_number = int(np.round(self.particle_number_per_size_class(k + 1)))
                        pset_new = ParticleSet(fieldset=fieldset, pclass=self.particle,
                                               lon=utils.create_list(particle.lon, particle_number),
                                               lat=utils.create_list(particle.lat, particle_number),
                                               depth=utils.create_list(particle.depth, particle_number),
                                               size=utils.create_list(new_particle_size, particle_number),
                                               parent=utils.create_list(particle.id, particle_number),
                                               age=utils.create_list(0, particle_number),
                                               beach=utils.create_list(particle.beach, particle_number),
                                               time=utils.create_list(particle.time, particle_number),
                                               rho_plastic=utils.create_list(particle.rho_plastic, particle_number),
                                               rise_velocity=utils.create_list(utils.initial_estimate_particle_rise_velocity(L=new_particle_size), particle_number),
                                               reynolds=utils.create_list(0, particle_number),
                                               size_class=utils.create_list(parent_size_class + k + 1, particle_number),
                                               prob_resus=utils.create_list(utils.resuspension_probability(w_rise=utils.initial_estimate_particle_rise_velocity(L=new_particle_size), time_step=self.dt), particle_number),
                                               repeatdt=None)
                        pset.add(pset_new)
        return pset

    def mass_per_size_class(self, k, f, p=settings.P_FRAG):
        gamma_ratio = math.gamma(k + f) / (math.gamma(k + 1) * math.gamma(f))
        return gamma_ratio * p ** k * (1 - p) ** f

    def particle_number_per_size_class(self, k, f=1, p=settings.P_FRAG, Dn=settings.DN):
        return self.mass_per_size_class(k, f, p) * 2 ** (Dn * k)

    def size_class_limits(self, k_range=settings.SIZE_CLASS_NUMBER, init_size=settings.INIT_SIZE,
                          p_frag=settings.P_FRAG):
        return np.array([init_size * p_frag ** k for k in range(k_range)])

    def run(self):
        # Creating the particle set and output file
        pset = self.get_pset(fieldset=self.field_set, particle_type=self.particle,
                             var_dict=self.get_var_dict(), start_time=utils.get_start_end_time(time='start'),
                             repeat_dt=self.repeat_dt)
        pfile = pset.ParticleFile(name=self.file_names(new=True),
                                  outputdt=self.output_time_step)
        # Setting the random seed and defining the particle behavior
        utils.print_statement("Setting the random seed")
        utils.set_random_seed(seed=settings.SEED)
        utils.print_statement("Defining the particle behavior")
        behavior_kernel = self.get_particle_behavior(pset=pset)
        # Carrying out the execution of the simulation
        utils.print_statement("The actual execution of the run")
        time = utils.get_start_end_time(time='start')
        while time <= utils.get_start_end_time(time='end'):
            pset.execute(behavior_kernel, runtime=self.output_time_step, dt=self.dt,
                         recovery={ErrorCode.ErrorOutOfBounds: utils.delete_particle},
                         output_file=pfile)
            time += self.output_time_step
            pset = self.particle_splitter(self.field_set, pset)
            utils.print_statement('time = {}'.format(time))
        pfile.export()
        utils.print_statement("Run completed")