from parcels import FieldSet, ParticleSet, ErrorCode
import numpy as np
import settings as settings
import scenarios.base_scenario as base_scenario
import factories.fieldset_factory as fieldset_factory
from advection_scenarios import advection_files
import utils
from datetime import datetime, timedelta
import Analysis
from parcels import ParcelsRandom
import math

from utils.physics_utils import particle_number_per_size_class


class FragmentationKaandorpPartial(base_scenario.BaseScenario):
    """Fragmentation scenario based on the Kaandorp fragmentation model. Beaching based on the stochastic scenario"""

    def __init__(self, server, stokes):
        """Constructor for FragmentationKaandorp"""
        super().__init__(server, stokes)
        self.prefix = "Frag_Kaandorp"
        self.input_dir = utils.get_input_directory(server=self.server)
        self.output_dir = utils.get_output_directory(server=self.server)
        self.repeat_dt = None
        self.OCEAN_FRAG = settings.OCEAN_FRAG
        if settings.SUBMISSION in ['simulation', 'visualization']:
            advection_scenario = advection_files.AdvectionFiles(server=self.server, stokes=self.stokes,
                                                                advection_scenario=settings.ADVECTION_DATA,
                                                                repeat_dt=self.repeat_dt)
            self.file_dict = advection_scenario.file_names
            if settings.SUBMISSION in ['simulation'] and not settings.POST_PROCESS:
                self.field_set = self.create_fieldset()

    var_list = ['lon', 'lat', 'beach', 'age', 'size', 'parent', 'beach_time', 'size_class', 'ocean_time']

    def create_fieldset(self) -> FieldSet:
        utils.print_statement("Creating the fieldset")
        fieldset = fieldset_factory.FieldSetFactory().create_fieldset(file_dict=self.file_dict, stokes=self.stokes,
                                                                      stokes_depth=True,
                                                                      border_current=True, diffusion=True, landID=True,
                                                                      distance=True, salinity=True, temperature=True,
                                                                      bathymetry=True, beach_timescale=True,
                                                                      resus_timescale=True, MLD=True,
                                                                      physics_constants=True, wind=True,
                                                                      TIDAL_mixing=True, fragmentation_timescale=True
                                                                      )
        return fieldset

    def get_pset(self, fieldset: FieldSet, particle_type: utils.BaseParticle, var_dict: dict,
                 start_time: datetime, repeat_dt: timedelta):
        """
        :return:
        """
        utils.print_statement("Creating the particle set")
        if settings.RESTART == 0:
            rise_velocity = utils.initial_estimate_particle_rise_velocity(L=var_dict['size'])
            rho_plastic = np.ones(rise_velocity.shape, dtype=np.float32) * settings.INIT_DENSITY
            pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                               lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                               age=var_dict['age'], size=var_dict['size'],
                               rho_plastic=rho_plastic, parent=range(len(rise_velocity)),
                               rise_velocity=rise_velocity, beach_time=np.zeros(rise_velocity.shape, dtype=np.int32),
                               prob_resus=utils.resuspension_probability(w_rise=rise_velocity),
                               size_class=np.zeros(rise_velocity.shape, dtype=np.float32),
                               ocean_time=np.zeros(rise_velocity.shape, dtype=np.int32),
                               time=start_time, repeatdt=repeat_dt)
        else:
            rise_velocity = utils.initial_estimate_particle_rise_velocity(L=var_dict['size'])
            rho_plastic = np.ones(rise_velocity.shape, dtype=np.float32) * settings.INIT_DENSITY
            pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                               lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'], parent=var_dict['parent'],
                               age=var_dict['age'], size=var_dict['size'], beach_time=var_dict['beach_time'],
                               rho_plastic=rho_plastic, rise_velocity=rise_velocity,
                               prob_resus=utils.resuspension_probability(w_rise=rise_velocity),
                               size_class=var_dict['size_class'], ocean_time=var_dict['ocean_time'],
                               time=start_time, repeatdt=repeat_dt)
        return pset

    def get_pclass(self):
        utils.print_statement("Creating the particle class")
        particle_type = utils.BaseParticle
        utils.add_particle_variable(particle_type, 'density', dtype=np.float32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'surface_density', dtype=np.float32, set_initial=False,
                                    to_write=False)
        utils.add_particle_variable(particle_type, 'kinematic_viscosity', dtype=np.float32, set_initial=False,
                                    to_write=False)
        utils.add_particle_variable(particle_type, 'rise_velocity', dtype=np.float32, set_initial=True, to_write=False)
        utils.add_particle_variable(particle_type, 'reynolds', dtype=np.float32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'rho_plastic', dtype=np.float32, set_initial=True, to_write=False)
        utils.add_particle_variable(particle_type, 'size', dtype=np.float32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'to_split', dtype=np.int32, set_initial=False, to_write=True)
        utils.add_particle_variable(particle_type, 'parent', dtype=np.int32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'prob_resus', dtype=np.float32, set_initial=True, to_write=False)
        utils.add_particle_variable(particle_type, 'beach_time', dtype=np.int32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'ocean_time', dtype=np.int32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'adrift_time', dtype=np.float32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'size_class', dtype=np.int32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'potential', dtype=np.float32, set_initial=False, to_write=False)
        return particle_type

    def file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART,
                   shore_time=settings.SHORE_TIME, ensemble=settings.ENSEMBLE,
                   advection_data=settings.ADVECTION_DATA, year=settings.START_YEAR + settings.RESTART,
                   month=settings.START_MONTH, input=settings.INPUT,
                   p_frag=settings.P_FRAG, dn=settings.DN, size_class_number=settings.SIZE_CLASS_NUMBER,
                   lambda_frag=settings.LAMBDA_FRAG, density=settings.INIT_DENSITY, postprocess=settings.POST_PROCESS):
        odirec = self.output_dir + "Kaandorp_Fragmentation_Partial/st_{}_e_{}/".format(shore_time, ensemble)
        if new:
            str_format = (advection_data, shore_time, p_frag, lambda_frag, dn, size_class_number, density, input, year,
                          month, restart, run)
        else:
            str_format = (advection_data, shore_time, p_frag, lambda_frag, dn, size_class_number, density, input,
                          year - 1,month, restart - 1, run)
        if postprocess:
            prefix = self.prefix + '_PP'
        else:
            prefix = self.prefix
        return odirec + prefix + '_{}_st={}_pfrag={}_lambdafrag={}_dn={}_sizeclasses={}_rho={}_I={}_y={}-{}_r={}_run={}.nc'.format(*str_format)

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
            particle.ocean_time += particle.dt
            dist = fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
            if dist < fieldset.Coastal_Boundary:
                if ParcelsRandom.uniform(0, 1) > fieldset.p_beach:
                    particle.beach = 1
        # Next the resuspension of particles on the coastline
        elif particle.beach == 1:
            particle.beach_time += particle.dt
            if ParcelsRandom.uniform(0, 1) > particle.prob_resus:
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
                         pset.Kernel(utils.vertical_reflecting_boundary) + \
                         pset.Kernel(self.beaching_kernel) + \
                         pset.Kernel(self.fragmentation_kernel)
        return total_behavior

    if settings.OCEAN_FRAG:
        def fragmentation_kernel(particle, fieldset, time):
            if particle.beach == 0:
                if particle.ocean_time >= (60 * 86400):
                    particle.to_split = 2
                    particle.ocean_time = 0
            elif particle.beach == 1:
                if particle.beach_time >= (60 * 86400):
                    particle.to_split = 1
                    particle.beach_time = 0
    else:
        def fragmentation_kernel(particle, fieldset, time):
            if particle.beach == 1:
                if particle.beach_time >= (60 * 86400):
                    particle.to_split = 1
                    particle.beach_time = 0

    def particle_splitter(self, fieldset, pset):
        for particle in pset:
            if particle.to_split > 0:
                # First, we set the split condition statement back to 0
                particle.to_split = 0
                # Getting the properties of the parent particle
                parent_size = particle.size
                parent_size_class = particle.size_class
                # Looping through the new particles being created, where new particles are only being created if the
                # parent particle
                remaining_classes = settings.SIZE_CLASS_NUMBER - parent_size_class - 1
                if remaining_classes > 0:
                    for k in range(0, remaining_classes):
                        new_particle_size = parent_size * 0.5 ** (k + 1)
                        particle_w_rise = utils.initial_estimate_particle_rise_velocity(L=new_particle_size)
                        pset_new = ParticleSet(fieldset=fieldset, pclass=self.particle,
                                               lon=particle.lon,
                                               lat=particle.lat,
                                               depth=particle.depth,
                                               size=new_particle_size,
                                               parent=particle.id,
                                               age=0,
                                               beach=particle.beach,
                                               time=particle.time,
                                               rho_plastic=particle.rho_plastic,
                                               rise_velocity=particle_w_rise,
                                               reynolds=0,
                                               prob_resus=utils.resuspension_probability(w_rise=particle_w_rise),
                                               size_class=parent_size_class + (k + 1),
                                               beach_time=0,
                                               ocean_time=0,
                                               repeatdt=None)
                        pset.add(pset_new)
        return pset

    def run(self):
        if settings.POST_PROCESS:
            utils.print_statement("Running postprocessing for LAMBDA_FRAG = {}".format(settings.LAMBDA_FRAG))
            base_file = self.file_names(new=True, lambda_frag=388, postprocess=False)
            output_file = self.file_names(new=True)
            Analysis.parcels_to_particle_number(base_file=base_file, output_file=output_file)
            utils.print_statement("Postprocessing is complete")
        else:
            # Creating the particle set and output file
            pset = self.get_pset(fieldset=self.field_set, particle_type=self.particle,
                                 var_dict=self.get_var_dict(), start_time=utils.get_start_end_time(time='start'),
                                 repeat_dt=self.repeat_dt)
            pfile = pset.ParticleFile(name=self.file_names(new=True),
                                      outputdt=settings.OUTPUT_TIME_STEP)
            # Setting the random seed and defining the particle behavior
            utils.print_statement("Setting the random seed")
            utils.set_random_seed(seed=settings.SEED)
            utils.print_statement("Defining the particle behavior")
            behavior_kernel = self.get_particle_behavior(pset=pset)
            # Carrying out the execution of the simulation
            utils.print_statement("The actual execution of the run")
            time = utils.get_start_end_time(time='start')
            while time <= utils.get_start_end_time(time='end'):
                pset.execute(behavior_kernel, runtime=settings.OUTPUT_TIME_STEP, dt=settings.TIME_STEP,
                             recovery={ErrorCode.ErrorOutOfBounds: utils.delete_particle},
                             output_file=pfile)
                time += settings.OUTPUT_TIME_STEP
                pset = self.particle_splitter(self.field_set, pset)
                utils.print_statement('time = {}'.format(time))
            pfile.export()
            utils.print_statement("Run completed")