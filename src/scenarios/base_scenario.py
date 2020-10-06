from abc import ABC, abstractmethod
from datetime import timedelta

from parcels import FieldSet, JITParticle, ParticleSet, ErrorCode
import numpy as np
from netCDF4 import Dataset
from utils import _set_random_seed, _delete_particle,_nan_removal,_get_start_end_time,_get_repeat_dt
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
        self.field_set = self.create_fieldset()
        self.particle = self._get_pclass()

    @property
    def var_list(self):
        raise NotImplementedError

    @abstractmethod
    def create_fieldset(self) -> FieldSet:
        pass

    # @abstractmethod
    # def create_particle_set(self) -> JITParticle:
    #     pass

    @abstractmethod
    def _get_pset(self) -> ParticleSet:
        pass

    @abstractmethod
    def _get_pclass(self)-> ParticleSet:
        pass

    @abstractmethod
    def _file_names(self, input_dir: str, new: bool) -> str:
        pass

    @abstractmethod
    def _beaching_kernel(self) -> ParticleSet:
        pass

    @abstractmethod
    def _get_particle_behavior(self):
        pass

    def _get_var_dict(self) -> dict:
        if settings.RESTART == 0:
            return pvf.initialize_variable_dict_from_varlist(self.var_list)
        else:
            return self._get_restart_variables(rfile=self.self._file_names(new=False), var_list=self.var_list)

    def run(self) -> object:
        pset = self._get_pset(fieldset=self.field_set, particle_type=self.particle,
                              var_dict=self._get_var_dict(),start_time=_get_start_end_time(time='start'),
                              repeat_dt=_get_repeat_dt())
        pfile = pset.ParticleFile(name=self._file_names(new=True),
                                  outputdt=settings.OUTPUT_TIME_STEP)
        os.system('echo "Setting the random seed"')
        _set_random_seed(seed=settings.SEED)
        os.system('echo "Defining the particle behavior"')
        behavior_kernel = self._get_particle_behavior(pset)
        os.system('echo "Setting the output file"')
        os.system('echo "Determine the simulation length"')
        _, _, simulation_length = _get_start_end_time()
        os.system('echo "The actual execution of the run"')
        pset.execute(behavior_kernel,
                     runtime=timedelta(days=_get_start_end_time(time='length')),
                     dt=settings.TIME_STEP,
                     recovery={ErrorCode.ErrorOutOfBounds: _delete_particle},
                     output_file=pfile
                     )
        pfile.export()

    def _get_restart_variables(self):
        dataset = Dataset(self._file_names(new=self._file_names(new=False)))
        time = dataset.variables['time'][:]
        final_time = time[0, -1]
        last_selec = np.ma.notmasked_edges(time, axis=1)[1]
        last_time_selec = time[last_selec[0], last_selec[1]]
        var_dict = {}
        for var in self.var_list:
            var_dict[var] = _nan_removal(dataset, var, last_selec, final_time, last_time_selec)
        return var_dict
