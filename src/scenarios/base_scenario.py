from abc import ABC, abstractmethod
from datetime import timedelta
from parcels import FieldSet, JITParticle, ParticleSet, ErrorCode
import numpy as np
from netCDF4 import Dataset
from utils import set_random_seed, delete_particle, restart_nan_removal, get_start_end_time
import utils
import settings as settings
from factories.pset_variable_factory import PsetVariableFactory as pvf
from advection_scenarios import advection_files


class BaseScenario(ABC):
    """A base class for the different scenarios"""

    def __init__(self, server, stokes):
        self.server: int = server
        self.stokes: int = stokes
        self.input_dir: str = settings.DATA_INPUT_DIREC
        self.output_dir: str = settings.DATA_OUTPUT_DIREC
        self.particle: ParticleSet = self.get_pclass()
        self.prefix: str = self.set_prefix()
        self.dt, self.output_time_step, self.repeat_dt = self.set_time_steps()
        self.var_list: list = self.set_var_list()
        if settings.SUBMISSION in ['simulation', 'visualization']:
            advection_scenario = advection_files.AdvectionFiles(server=self.server, stokes=self.stokes,
                                                                advection_scenario=settings.ADVECTION_DATA,
                                                                repeat_dt=self.repeat_dt)
            self.file_dict = advection_scenario.file_names
            if settings.SUBMISSION in ['simulation']:
                self.field_set = self.create_fieldset()

    @abstractmethod
    def set_prefix(self) -> str:
        pass

    @abstractmethod
    def set_time_steps(self) -> tuple:
        pass

    @abstractmethod
    def set_var_list(self) -> list:
        pass

    @abstractmethod
    def create_fieldset(self) -> FieldSet:
        pass

    @abstractmethod
    def get_pset(self) -> ParticleSet:
        pass

    @abstractmethod
    def get_pclass(self) -> ParticleSet:
        pass

    @abstractmethod
    def file_names(self, input_dir: str, new: bool) -> str:
        pass

    @abstractmethod
    def beaching_kernel(self) -> ParticleSet.Kernel:
        pass

    @abstractmethod
    def get_particle_behavior(self) -> ParticleSet.Kernel:
        pass

    def get_restart_variables(self) -> dict:
        dataset = Dataset(self.file_names(new=False))
        time = dataset.variables['time'][:]
        final_time = time[0, -1]
        last_selec = np.ma.notmasked_edges(time, axis=1)[1]
        last_time_selec = time[last_selec[0], last_selec[1]]
        var_dict = {}
        for var in self.var_list:
            var_dict[var] = restart_nan_removal(dataset, var, last_selec, final_time, last_time_selec)
        return var_dict

    def get_var_dict(self) -> dict:
        if settings.RESTART == 0:
            return pvf.initialize_variable_dict_from_varlist(var_list=self.var_list,
                                                             start_files=self.file_dict['STARTFILES_filename'])
        else:
            return self.get_restart_variables()

    def run(self) -> object:
        utils.print_statement("Creating the particle set")
        pset = self.get_pset(fieldset=self.field_set, particle_type=self.particle,
                             var_dict=self.get_var_dict(), start_time=get_start_end_time(time='start'),
                             repeat_dt=self.repeat_dt)
        pfile = pset.ParticleFile(name=self.file_names(new=True),
                                  outputdt=self.output_time_step)
        utils.print_statement("Setting the random seed")
        set_random_seed(seed=settings.SEED)
        utils.print_statement("Defining the particle behavior")
        behavior_kernel = self.get_particle_behavior(pset=pset)
        utils.print_statement("The actual execution of the run")
        pset.execute(behavior_kernel,
                     runtime=timedelta(days=get_start_end_time(time='length')),
                     dt=self.dt,
                     recovery={ErrorCode.ErrorOutOfBounds: delete_particle},
                     output_file=pfile
                     )
        pfile.export()
        utils.print_statement("Run completed")

    def return_full_run_directory(self) -> dict:
        """
        Return a directory with all file names depending on the restart and run variables
        :return:
        """
        file_dict = {}
        for run in range(settings.RUN_RANGE):
            restart_direc = {}
            for restart in range(settings.SIM_LENGTH):
                restart_direc[restart] = self.file_names(new=True, run=run, restart=restart)
            file_dict[run] = restart_direc
        return file_dict
