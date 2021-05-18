from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
from progressbar import ProgressBar
import os
import settings
import utils


def parcels_to_timeslicing(file_dict: dict):
    """
    For easy animation making, we make slices of the particle lon, lat, depth and beach status so that we can make
    animations even when loading all the files in their entirety would overwhelm my computer RAM
    :param file_dict: a dictionary containing all file names for a specific simulation
    :return:
    """
    # Getting the directory saving the output files
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeslices/{}/'.format(
        settings.SCENARIO_NAME)
    utils.check_direc_exist(output_direc)


    # loop through the runs
    pbar = ProgressBar
    for run in pbar(range(settings.RUN_RANGE)):
        # Loop through the restart files
        for restart in range(settings.SIM_LENGTH):
            # Getting the load file
            parcels_file = file_dict[run][restart]
            dataset = Dataset(parcels_file)



