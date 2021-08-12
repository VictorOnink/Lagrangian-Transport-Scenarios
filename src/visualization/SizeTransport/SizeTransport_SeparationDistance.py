import settings
import factories.scenario_factory as scenario_factory
import utils
from netCDF4 import Dataset
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import numpy as np
import string
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from advection_scenarios import advection_files
from copy import deepcopy
import os


def SizeTransport_SeparationDistance(scenario, figure_direc, size_selection, rho_selection, tau_selection, size_list,
                                     figsize=(18, 8), fontsize=14, legendsize=12):
    # Setting the folder within which we have the output, and where we have the saved data
    output_direc = figure_direc + 'separation_distance/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'separation_distance/{}/'.format('SizeTransport')
    utils.check_direc_exist(output_direc)

    # Loading the separation distances for the particle size we've selected
    prefix = 'separation_distance'
    data_dict = vUtils.SizeTransport_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                               size=size_selection, rho=rho_selection, tau=tau_selection)
    utils.print_statement(str(data_dict['MEAN'][0].keys()), to_print=True)
