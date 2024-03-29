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


def var_dict_expansion(var_dict: dict):
    if settings.INPUT in ['LebretonKaandorpInit']:
        array_size = var_dict['lon'].size
        for variable in var_dict.keys():
            if variable == 'size_class':
                for size_class in range(1, settings.SIZE_CLASS_NUMBER):
                    var_dict[variable] = np.append(var_dict[variable], np.ones(array_size, dtype=np.float32) * size_class)
            else:
                var_dict[variable] = np.tile(var_dict[variable], settings.SIZE_CLASS_NUMBER)
    return var_dict


class FragmentationKaandorpPartial(base_scenario.BaseScenario):
    """Fragmentation scenario based on the Kaandorp fragmentation model. Beaching based on the stochastic scenario"""

    def __init__(self):
        """Constructor for FragmentationKaandorp"""
        super().__init__()
        # Whether to include ocean fragmentation
        self.OCEAN_FRAG = settings.OCEAN_FRAG
        # Beached time cutoff for fragmentation to occur (days)
        self.T_frag = 90
        self.advection_scenario = advection_files.AdvectionFiles(repeat_dt=self.repeat_dt)
        self.file_dict = self.advection_scenario.file_names
        if settings.SUBMISSION in ['simulation'] and not settings.POST_PROCESS:
            self.field_set = self.create_fieldset()

    def set_prefix(self) -> str:
        """
        Set the scenario advection_prefix
        :return:
        """
        return "Frag_Kaandorp"

    def set_var_list(self) -> list:
        """
        Set the var_list, which contains all the variables that need to be loaded during the restarts
        :return:
        """
        return ['lon', 'lat', 'beach', 'age', 'parent', 'beach_time', 'size_class', 'ocean_time', 'at_seafloor',
                'distance2coast']

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
        fieldset = fieldset_factory.FieldSetFactory().create_fieldset(file_dict=self.file_dict,
                                                                      stokes_depth=True,
                                                                      border_current=True, diffusion=True, landID=True,
                                                                      distance=True, salinity=True, temperature=True,
                                                                      bathymetry=True, beach_timescale=True,
                                                                      resus_timescale=True,
                                                                      fixed_resus=settings.FIXED_RESUS,
                                                                      MLD=True, physics_constants=True, wind=True,
                                                                      TIDAL_mixing=True,
                                                                      fragmentation_period=self.T_frag,
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
            var_dict = var_dict_expansion(var_dict=var_dict)
            particle_size = settings.INIT_SIZE * np.power(2, -1 * var_dict['size_class'])
            rise_velocity = utils.initial_estimate_particle_rise_velocity(L=particle_size)
            rho_plastic = np.ones(rise_velocity.shape, dtype=np.float32) * settings.INIT_DENSITY
            pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                               lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                               age=var_dict['age'], size=particle_size,
                               rho_plastic=rho_plastic, parent=range(len(rise_velocity)),
                               rise_velocity=rise_velocity, beach_time=np.zeros(rise_velocity.shape, dtype=np.int32),
                               prob_resus=utils.resuspension_probability(w_rise=rise_velocity, time_step=self.dt),
                               size_class=var_dict['size_class'],
                               ocean_time=np.zeros(rise_velocity.shape, dtype=np.int32),
                               at_seafloor=np.zeros(rise_velocity.shape, dtype=np.int32),
                               distance2coast=np.zeros(rise_velocity.shape, dtype=np.float32),
                               time=start_time, repeatdt=repeat_dt)
        else:
            particle_size = settings.INIT_SIZE * np.power(2, -1 * var_dict['size_class'])
            rise_velocity = utils.initial_estimate_particle_rise_velocity(L=particle_size)
            rho_plastic = np.ones(rise_velocity.shape, dtype=np.float32) * settings.INIT_DENSITY
            pset = ParticleSet(fieldset=fieldset, pclass=particle_type,
                               lon=var_dict['lon'], lat=var_dict['lat'], beach=var_dict['beach'],
                               parent=var_dict['parent'],
                               age=var_dict['age'], size=particle_size, beach_time=var_dict['beach_time'],
                               rho_plastic=rho_plastic, rise_velocity=rise_velocity,
                               prob_resus=utils.resuspension_probability(w_rise=rise_velocity, time_step=self.dt),
                               size_class=var_dict['size_class'], ocean_time=var_dict['ocean_time'],
                               at_seafloor=var_dict['at_seafloor'], distance2coast=var_dict['distance2coast'],
                               time=start_time, repeatdt=repeat_dt)
        return pset

    def get_pclass(self):
        utils.print_statement("Creating the particle class")
        particle_type = utils.BaseParticle
        utils.add_particle_variable(particle_type, 'density', dtype=np.float32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'distance2coast', dtype=np.float32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'surface_density', dtype=np.float32, set_initial=False,
                                    to_write=False)
        utils.add_particle_variable(particle_type, 'kinematic_viscosity', dtype=np.float32, set_initial=False,
                                    to_write=False)
        utils.add_particle_variable(particle_type, 'rise_velocity', dtype=np.float32, set_initial=True, to_write=False)
        utils.add_particle_variable(particle_type, 'reynolds', dtype=np.float32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'rho_plastic', dtype=np.float32, set_initial=True, to_write=False)
        utils.add_particle_variable(particle_type, 'size', dtype=np.float32, set_initial=True, to_write=False)
        utils.add_particle_variable(particle_type, 'to_split', dtype=np.int32, set_initial=False, to_write=True)
        utils.add_particle_variable(particle_type, 'parent', dtype=np.int32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'prob_resus', dtype=np.float32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'beach_time', dtype=np.int32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'ocean_time', dtype=np.int32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'size_class', dtype=np.int32, set_initial=True, to_write=True)
        utils.add_particle_variable(particle_type, 'potential', dtype=np.float32, set_initial=False, to_write=False)
        utils.add_particle_variable(particle_type, 'at_seafloor', dtype=np.int32, set_initial=True, to_write=True)
        return particle_type

    def file_names(self, new: bool = False, run: int = settings.RUN, restart: int = settings.RESTART,
                   shore_time=settings.SHORE_TIME, ensemble=settings.ENSEMBLE,
                   advection_data=settings.ADVECTION_DATA, year=settings.STARTYEAR,
                   month=settings.STARTMONTH, input=settings.INPUT, ocean_frag=settings.OCEAN_FRAG, ocean_lambda=settings.LAMBDA_OCEAN_FRAG,
                   p_frag=settings.P_FRAG, dn=settings.DN, size_class_number=settings.SIZE_CLASS_NUMBER,
                   lambda_frag=settings.LAMBDA_FRAG, density=settings.INIT_DENSITY, postprocess=settings.POST_PROCESS):
        odirec = self.output_dir + "Kaandorp_Fragmentation_Partial/st_{}_e_{}/".format(shore_time, ensemble)
        if postprocess:
            odirec += 'lambda_f={}/'.format(lambda_frag)
        else:
            odirec += 'parcels_output/'
        if new:
            str_format = (advection_data, shore_time, p_frag, lambda_frag, dn, size_class_number, density, input, year,
                          month, restart, run)
        else:
            str_format = (advection_data, shore_time, p_frag, lambda_frag, dn, size_class_number, density, input,
                          year, month, restart - 1, run)
        prefix = {True: self.prefix + '_PP', False: self.prefix}[postprocess]
        if ocean_frag:
            prefix += '_OFRAG_{}'.format(ocean_lambda)
        file_type = {True: '.pkl', False: '.nc'}[postprocess]
        return odirec + prefix + '_{}_st={}_pfrag={}_lambdafrag={}_dn={}_sizeclasses={}_rho={}_I={}_y={}-{}_r={}_run={}'.format(*str_format) + file_type

    def beaching_kernel(particle, fieldset, time):
        """
        The basic beaching and resuspension procedure Onink et al. (2021). However, since Hinata et al. (2017) showed
        that the resuspension varies with size, we have expande
        Onink et al. (2021) = https://doi.org/10.1088/1748-9326/abecbd
        Hinata et al. (2017) = https://doi.org/10.1016/j.marpolbul.2017.05.012
        """
        # First, the beaching of particles on the coastline
        if particle.beach == 0:
            particle.ocean_time += particle.dt
            particle.distance2coast = fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
            if particle.distance2coast < fieldset.Coastal_Boundary:
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
                         pset.Kernel(self.fragmentation_kernel) + \
                         pset.Kernel(self.at_seafloor)
        return total_behavior

    if settings.OCEAN_FRAG:
        def fragmentation_kernel(particle, fieldset, time):
            if particle.beach == 0:
                if particle.ocean_time >= (fieldset.T_frag * 86400):
                    particle.to_split += 2
                    particle.ocean_time = 0
            elif particle.beach == 1:
                if particle.beach_time >= (fieldset.T_frag * 86400):
                    particle.to_split += 1
                    particle.beach_time = 0
    else:
        def fragmentation_kernel(particle, fieldset, time):
            if particle.beach == 1:
                if particle.beach_time >= (fieldset.T_frag * 86400):
                    particle.to_split += 1
                    particle.beach_time = 0

    def at_seafloor(particle, fieldset, time):
        if particle.beach == 0:
            local_bathymetry = fieldset.bathymetry[time, fieldset.SURF_Z, particle.lat, particle.lon]
            if math.fabs(local_bathymetry - particle.depth) < 1:
                particle.at_seafloor += particle.dt

    def particle_splitter(self, fieldset, pset):
        # Initialize the dictionary containing the new particles
        new_dict = {'count': 0, 'lon': np.array([]), 'lat': np.array([]), 'depth': np.array([]), 'size': np.array([]),
                    'parent': np.array([]), 'age': np.array([]), 'beach': np.array([]), 'time': np.array([]),
                    'rho_plastic': np.array([]), 'rise_velocity': np.array([]), 'reynolds': np.array([]),
                    'prob_resus': np.array([]), 'size_class': np.array([]), 'beach_time': np.array([]),
                    'ocean_time': np.array([]), 'at_seafloor': np.array([])}
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
                        new_dict['count'] += 1
                        new_dict['lon'] = np.append(new_dict['lon'], particle.lon)
                        new_dict['lat'] = np.append(new_dict['lat'], particle.lat)
                        new_dict['depth'] = np.append(new_dict['depth'], particle.depth)
                        new_dict['size'] = np.append(new_dict['size'], parent_size * 0.5 ** (k + 1))
                        new_dict['parent'] = np.append(new_dict['parent'], particle.id)
                        new_dict['age'] = np.append(new_dict['age'], 0)
                        new_dict['beach'] = np.append(new_dict['beach'], particle.beach)
                        new_dict['time'] = np.append(new_dict['time'], particle.time)
                        new_dict['rho_plastic'] = np.append(new_dict['rho_plastic'], particle.id)
                        new_dict['rise_velocity'] = np.append(new_dict['rise_velocity'], utils.initial_estimate_particle_rise_velocity(L=new_dict['size'][-1]))
                        new_dict['reynolds'] = np.append(new_dict['reynolds'], 0)
                        new_dict['prob_resus'] = np.append(new_dict['prob_resus'], utils.resuspension_probability(w_rise=new_dict['rise_velocity'][-1], time_step=self.dt))
                        new_dict['size_class'] = np.append(new_dict['size_class'], parent_size_class + (k + 1))
                        new_dict['beach_time'] = np.append(new_dict['beach_time'], 0)
                        new_dict['ocean_time'] = np.append(new_dict['ocean_time'], 0)
                        new_dict['at_seafloor'] = np.append(new_dict['at_seafloor'], 0)

        # if there are new particles, create a new particle set and add this to the pset object
        if new_dict['count'] > 0:
            pset_new = ParticleSet(fieldset=fieldset, pclass=self.particle, lon=new_dict['lon'], lat=new_dict['lat'],
                                   depth=new_dict['depth'], size=new_dict['size'], parent=new_dict['parent'],
                                   age=new_dict['age'], beach=new_dict['beach'], time=new_dict['time'],
                                   rho_plastic=new_dict['rho_plastic'], rise_velocity=new_dict['rise_velocity'],
                                   reynolds=new_dict['reynolds'], prob_resus=new_dict['prob_resus'],
                                   size_class=new_dict['size_class'], beach_time=new_dict['beach_time'],
                                   ocean_time=new_dict['ocean_time'], at_seafloor=new_dict['at_seafloor'],
                                   repeatdt=None)
            pset.add(pset_new)
        return pset

    def run(self):
        if settings.POST_PROCESS:
            utils.print_statement("Running postprocessing for LAMBDA_FRAG = {}".format(settings.LAMBDA_FRAG))
            base_file = self.file_names(new=True, lambda_frag=388, postprocess=False)
            output_file, restart_file = self.file_names(new=True), self.file_names(new=False)
            Analysis.parcels_to_particle_number(base_file=base_file, output_file=output_file, restart_file=restart_file,
                                                time_step=self.dt).run()
            utils.print_statement("Postprocessing is complete")
        else:
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

    def return_full_run_directory(self) -> dict:
        """
        Return a directory with all file names depending on the restart and run variables, and also for the different
        starting months and years
        We also have a separation to account that we can have both parcels files and post-processing output files
        :return:
        """
        file_types = ['parcels', 'postprocess']
        file_dict = {'parcels': {}, 'postprocess': {}}

        for ind_type, file_key in enumerate(file_types):
            for ind_year, year in enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH)):
                file_dict[file_key][year] = {}
                for month in range(1, 13):
                    file_dict[file_key][year][month] = {}
                    for run in range(settings.RUN_RANGE):
                        file_dict[file_key][year][month][run] = {}
                        for restart in range(settings.SIM_LENGTH - ind_year):
                            if file_key == 'parcels':
                                lambda_frag = 388
                            else:
                                lambda_frag = settings.LAMBDA_FRAG
                            file_dict[file_key][year][month][run][restart] = self.file_names(new=True, run=run, restart=restart,
                                                                                             year=year, month=month,
                                                                                             postprocess=settings.BOOLEAN_DICT[ind_type],
                                                                                             lambda_frag=lambda_frag)
        return file_dict
