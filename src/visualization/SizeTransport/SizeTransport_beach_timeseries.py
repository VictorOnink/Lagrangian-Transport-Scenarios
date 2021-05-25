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
    output_direc = figure_direc + 'timeseries/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/{}/'.format('SizeTransport')
    utils.check_direc_exist(output_direc)

    # Loading in the data
    prefix = 'timeseries'
    timeseries_dict = {}
    beach_state_list = ['beach', 'afloat', 'seabed', 'removed', 'total']
    for index, size in enumerate(size_list):
        data_dict = vUtils.SizeTransport_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                                   size=size, rho=rho_list[index])
        timeseries_dict[size] = {}
        for beach_state in beach_state_list:
            timeseries_dict[size][beach_state] = data_dict[beach_state]
    time = data_dict['time']

    # Setting parameters for the time axis
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y')

    # Creating the figure
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=5, ncols=1)
    gs.update(wspace=1, hspace=1)

    ax_list = []
    for row in range(gs.nrows):
        ax = fig.add_subplot(gs[row, 0])
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_major_formatter(yearsFmt)
        ax.xaxis.set_minor_locator(months)
        ax.set_ylabel(r'Particle Count', fontsize=fontsize)
        ax.set_xlim(datetime(2010, 1, 1), datetime(2013, 1, 1))
        ax.tick_params(which='major', length=7)
        ax.tick_params(which='minor', length=3)
        ax.grid(True)
        ax_list.append(ax)
    ax_list[0].set_xlabel('Time (yr)', fontsize=fontsize)
    for index, beach_state in enumerate(beach_state_list):
        ax_list[index].set_title(subfigure_title(index, beach_state), fontsize=fontsize)

    file_name = output_direc + 'test.jpg'
    plt.savefig(file_name, bbox_inches='tight')


def subfigure_title(index, beach_state):
    """
    setting the title of the subfigure
    :param index:
    :param size:
    :param rho:
    :return:
    """
    alphabet = string.ascii_lowercase
    return '({}) {}'.format(alphabet[index], beach_state)
