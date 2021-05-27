import settings as settings
import utils
from advection_scenarios import advection_files
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import progressbar
import os
from copy import deepcopy


def parcels_to_maxminmean(file_dict: dict):
    """
    This is a general function that we can use to calculate max, min and mean values of particle
    :param file_dict:
    :return:
    """
