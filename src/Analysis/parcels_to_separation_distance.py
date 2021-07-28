from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
from progressbar import ProgressBar
import settings
import utils


def parcels_to_separation_distance(file_dict: dict, scenario):
    """
    This code is to calculate the average separation distance over time between two sets of particles
    This is currently being written to work just for the SizeTransport scenario, but future versions might be able to
    adapt that to be more flexible
    :param file_dict: dictionary containing all file runs
    :param scenario: this is the scenario object, from which we can get the file_names function
    :return:
    """
    # Getting the directory saving the output files
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'separation_distance/{}/'.format(
        settings.SCENARIO_NAME)
    utils.check_direc_exist(output_direc)
    # Creating a dictionary to save the average seperation distance
    time_steps = len(Dataset(scenario.file_names(new=True, run=0, restart=0)).variables['time'][0,:])
    output_dict = dict.fromkeys(range(time_steps))
    print(output_dict)

    # # loop through the runs
    # for run in range(settings.RUN_RANGE):
    #     # Loop through the restart files
    #     for restart in range(settings.SIM_LENGTH):
    #         # Getting the parcels output file
    #         parcels_file = file_dict[run][restart]
    #         dataset = Dataset(parcels_file)