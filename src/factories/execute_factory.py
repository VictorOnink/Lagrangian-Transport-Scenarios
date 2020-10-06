# import src.settings as settings
# from parcels import ErrorCode, ParticleSet
# import src.utils as utils
# from datetime import timedelta
# import parcels.rng as rng
# from time import time
# # from src.scenarios.coastal_proximity import _particle_behavior_proximity
# import src.scenarios.coastal_proximity as proximity
# from src.scenarios.base_scenario import DeleteParticle
# import os
#
# class ExecuteFactory:
#     """A factory class for the execution of the scenario run"""
#     def execute_scenario(cls,server: int, pset: ParticleSet, seed: str='Fixed'):
#         input_dir = utils._get_input_directory(server=server)
#
#         os.system('echo "Setting the random seed"')
#         _set_random_seed(seed=seed)
#         os.system('echo "Defining the particle behavior"')
#         behavior_kernel=_set_particle_behavior(pset=pset)
#         os.system('echo "Setting the output file"')
#         pfile = _get_pfile(pset=pset,input_dir=input_dir)
#         os.system('echo "Determine the simulation length"')
#         _, _, simulation_length = utils._get_start_end_time()
#         os.system('echo "The actual execution of the run"')
#         pset.execute(behavior_kernel,
#                      runtime=timedelta(days=simulation_length),
#                      dt=settings.TIME_STEP,
#                      recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
#                      output_file=pfile
#                      )
#         os.system('echo "Exporting the pfile at the end of the simulation"')
#         pfile.export()
#
#
#
#
# def _set_particle_behavior(pset: ParticleSet):
#     if settings.SCENARIO_NAME=='CoastalProximity':
#         particle_behavior = proximity._particle_behavior_proximity(pset)
#     return particle_behavior
#
# def _get_pfile(pset: ParticleSet,input_dir: str):
#     if settings.SCENARIO_NAME=='CoastalProximity':
#         ofile,_=proximity._file_names_proximity(input_dir)
#
#     pfile = pset.ParticleFile(name=ofile,
#                               outputdt=timedelta(hours=24))
#     return pfile