import settings as settings
import utils
from advection_scenarios import advection_files
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar
from copy import deepcopy


# class parcels_to_concentration:
#     def __init__(self, file_dict: dict):
#         self.parallel_step = settings.PARALLEL_STEP
#         self.file_dict = file_dict
#         self.LON, self.LAT, self.GRID, self.hexgrid = create_hex_grid()
#         self.beach_label_dict = set_beach_label_dict(scenario_name=settings.SCENARIO_NAME)
#         self.temp_direc, self.output_direc = get_directories(scenario_name=settings.SCENARIO_NAME)
#         self.weight_list = ['particle_mass_sink', 'particle_number_sink']
#         self.output_dict = create_output_file_dict(scenario_name=settings.SCENARIO_NAME, grid=self.GRID,
#                                                    beach_states=self.beach_label_dict.keys(), lon=self.LON,
#                                                    lat=self.LAT, weight_list=self.weight_list)
