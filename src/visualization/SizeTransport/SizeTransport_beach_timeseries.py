import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import FixedLocator, FixedFormatter


def SizeTransport_beach_timeseries(scenario, figure_direc, size_list, rho_list, figsize=(10, 10), fontsize=14):
    # Setting the folder within which we have the output, and where we have the saved timeslices
    output_direc = figure_direc + 'concentrations/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/{}/'.format('SizeTransport')
    utils.check_direc_exist(output_direc)

    # Loading in the data
    prefix = 'timeseries'
    timeseries_dict = {}
    for index, size in enumerate(size_list):
        data_dict = vUtils.SizeTransport_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                                   size=size, rho=rho_list[index])
        timeseries_dict[size] = {}
        for beach_state in ['beach', 'afloat', 'seabed', 'removed', 'total']:
            timeseries_dict[size][beach_state] = data_dict[beach_state]
    time = data_dict['time']

    # Setting parameters for the time axis
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y')

    # Creating the figure
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=5, ncols=1)

    ax_list = []
    for row in range(gs.nrows):
        ax = fig.add_subplot([row, 0])
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_major_formatter(yearsFmt)
        ax.xaxis.set_minor_locator(months)
        ax.set_ylabel(r'Particle Count', fontsize=fontsize)
        ax.set_xlim(datetime(2010, 1, 1), datetime(2010 + settings.RUN_RANGE, 1, 1))
        ax.tick_params(which='major', length=7)
        ax.tick_params(which='minor', length=3)
        ax.grid(True)
        ax_list.append(ax)
    ax_list[0].set_xlabel('Time (yr)', fontsize=fontsize)

    file_name = output_direc + 'test.jpg'
    plt.savefig(file_name, bbox_inches='tight')
