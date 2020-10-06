# import os
# import settings as settings
# import numpy as np
# from parcels import JITParticle, Variable, FieldSet
# from operator import attrgetter
# from datetime import timedelta, datetime
#
# from netCDF4 import Dataset
# from utils import _add_var_particle, _get_input_directory, _get_repeat_dt, _get_start_end_time
# from utils import BaseParticle
#
#
# class ParticleFactory:
#     """A factory class for the particles"""
#
#     @classmethod
#     def create_particle_set(cls, server: int, fieldset: FieldSet, rfile: str, distance=True,
#                             var_dict: dict):
#         """
#
#         :rtype: object
#         """
#         input_dir = _get_input_directory(server=server)
#
#         # Create the particle class
#         os.system('echo "Create the particle class"')
#         particle_type = BaseParticle
#         if settings.SCENARIO_NAME == 'CoastalProximity':
#             _add_var_particle(particle_type, 'prox')
#         if distance:
#             _add_var_particle(particle_type, 'distance', dtype=np.float32, set_initial=False)
#         # Get the particle start positions and other relevant variables
#         os.system('echo "Get the (re)start particle positions and other relevant variables"')
#         # Get timestep of repeated particle releases and the starting date of the simulation
#         os.system('echo "Set the timestep of particle release and simulation start time"')
#         repeat_dt = _get_repeat_dt()
#         start_time, _, _ = _get_start_end_time()
#         # Create the pset
#         pset = _create_pset(fieldset, particle_type, var_dict, start_time, repeat_dt)
#         return pset
#
#
