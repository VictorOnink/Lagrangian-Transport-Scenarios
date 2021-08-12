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
                                     fig_size=(18, 8), fontsize=14, legendsize=12, y_label='Distance (km)',
                                     x_label='Time'):
    # Setting the folder within which we have the output, and where we have the saved data
    output_direc = figure_direc + 'separation_distance/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'separation_distance/{}/'.format('SizeTransport')
    utils.check_direc_exist(output_direc)

    # Loading the separation distances for the particle size we've selected
    prefix = 'separation_distance'
    data_dict = vUtils.SizeTransport_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                               size=size_selection, rho=rho_selection, tau=tau_selection)
    mean_dict = data_dict['MEAN']
    median_dict = data_dict['MEDIAN']
    time = data_dict['TIME']

    # Setting parameters for the time axis
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y')

    # Getting the datetime objects for all of the time arrays
    startdate = datetime(settings.START_YEAR, 1, 1, 12, 0)
    time_list = []
    for t in time:
        time_list.append(startdate + timedelta(seconds=t))

    # Creating the figure
    ax_range = datetime(settings.START_YEAR, 1, 1), datetime(settings.START_YEAR + 1, 1, 1), 4000, 0
    ax = vUtils.base_figure(fig_size=fig_size, ax_range=ax_range, x_label=x_label, y_label=y_label,
                            ax_label_size=fontsize, shape=(1, 2), plot_num=2, legend_axis=True, log_yscale=True,
                            width_ratios=[1, 1, 0.5])

    # Saving the figure
    file_name = output_direc + 'Separation_distance_size={}.png'.format(size_selection)
    plt.savefig(file_name, bbox_inches='tight')

