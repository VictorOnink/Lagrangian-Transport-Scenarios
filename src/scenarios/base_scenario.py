from abc import ABC, abstractmethod
from datetime import timedelta
from parcels import FieldSet, JITParticle, ParticleSet, ErrorCode
import numpy as np
from netCDF4 import Dataset
from utils import set_random_seed, delete_particle, restart_nan_removal, get_start_end_time
import settings as settings
import os
from factories.pset_variable_factory import PsetVariableFactory as pvf


class BaseScenario(ABC):
    server: int
    stokes: int
    field_set: FieldSet
    particle: JITParticle
    prefix: str

    """A base class for the different scenarios"""

    def __init__(self, server, stokes):
        self.server = server
        self.stokes = stokes
        self.particle = self.get_pclass()

    @property
    def var_list(self):
        raise NotImplementedError

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
    def beaching_kernel(self) -> ParticleSet:
        pass

    @abstractmethod
    def get_particle_behavior(self):
        pass

    def get_restart_variables(self):
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
        os.system('echo "Creating the particle set"')
        pset = self.get_pset(fieldset=self.field_set, particle_type=self.particle,
                             var_dict=self.get_var_dict(), start_time=get_start_end_time(time='start'),
                             repeat_dt=self.repeat_dt)
        pfile = pset.ParticleFile(name=self.file_names(new=True),
                                  outputdt=settings.OUTPUT_TIME_STEP)
        os.system('echo "Setting the random seed"')
        set_random_seed(seed=settings.SEED)
        os.system('echo "Defining the particle behavior"')
        behavior_kernel = self.get_particle_behavior(pset=pset)
        os.system('echo "The actual execution of the run"')
        pset.execute(behavior_kernel,
                     runtime=timedelta(days=get_start_end_time(time='length')),
                     dt=settings.TIME_STEP,
                     recovery={ErrorCode.ErrorOutOfBounds: delete_particle},
                     output_file=pfile
                     )
        pfile.export()
        os.system('echo "Run completed"')

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
