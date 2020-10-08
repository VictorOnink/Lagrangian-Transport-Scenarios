"""
Converting the trajectory data from the parcels output to timeslices representing each date

"""
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import os
from progressbar import ProgressBar
import utils

class timeslicing():
    def __init__(self, server, stokes):
        self.server = server
        self.stokes = stokes
        self.input_dir = utils._get_input_directory(server=self.server)

