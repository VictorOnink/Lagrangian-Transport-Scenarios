import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import numpy as np
import string
from datetime import datetime, timedelta
import matplotlib.dates as mdates


def SizeTransport_beach_timeseries(scenario, figure_direc, size_list, rho_list, figsize=(10, 10), fontsize=12):
    # Setting the folder within which we have the output, and where we have the saved timeslices
    output_direc = figure_direc + 'timeseries/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/{}/'.format('SizeTransport')
    utils.check_direc_exist(output_direc)

    # Loading in the data
    prefix = 'timeseries'
    timeseries_dict = {}
    # beach_state_list = ['beach', 'afloat', 'seabed', 'removed', 'total']
    beach_state_list = ['beach', 'afloat', 'seabed']
    for index, size in enumerate(size_list):
        data_dict = vUtils.SizeTransport_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                                   size=size, rho=rho_list[index])
        timeseries_dict[size] = {}
        for beach_state in beach_state_list:
            timeseries_dict[size][beach_state] = data_dict[beach_state]
    time = data_dict['time']
    total = float(data_dict['total'][0])

    # Normalizing all the particle counts with the total number of particles, and then multiplying by 100 to get a
    # percentage
    for size in size_list:
        for beach_state in beach_state_list:
            timeseries_dict[size][beach_state] /= total
            timeseries_dict[size][beach_state] *= 100.

    # Setting parameters for the time axis
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y')

    # Getting the datetime objects for all of the time
    time_list = []
    startdate = datetime(settings.START_YEAR, 1, 1, 12, 0)
    for seconds in time:
        time_list.append(startdate + timedelta(seconds=seconds))

    # Creating the figure
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=3, ncols=1)
    gs.update(wspace=0.2, hspace=0.2)

    ax_list = []
    for row in range(gs.nrows):
        ax = fig.add_subplot(gs[row, 0])
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_minor_locator(months)
        ax.xaxis.set_major_formatter(yearsFmt)
        ax.set_ylabel(r'Fraction of Total (%)', fontsize=fontsize)
        ax.set_xlim(datetime(2010, 1, 1), datetime(2013, 1, 1))
        ax.set_ylim([0, 100])
        ax.tick_params(which='major', length=7)
        ax.tick_params(which='minor', length=3)
        if row != (gs.nrows - 1):
            ax.set_xticklabels([])
        ax_list.append(ax)
    ax_list[-1].set_xlabel('Time (yr)', fontsize=fontsize)

    for index, beach_state in enumerate(beach_state_list):
        ax_list[index].set_title(subfigure_title(index, beach_state), fontsize=fontsize)

    # Now, adding in the actual data
    for index_size, size in enumerate(size_list):
        for index_beach, beach_state in enumerate(beach_state_list):
            ax_list[index_beach].plot(time_list, timeseries_dict[size][beach_state], linestyle='-',
                                      color=vUtils.discrete_color_from_cmap(index_size, subdivisions=len(size_list)),
                                      label=size_label(size))
    # And adding in a legend
    ax_list[0].legend(fontsize=fontsize, loc='upper right')

    file_name = output_direc + 'SizeTransport_beach_state_timeseries.jpg'
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

def size_label(size):
    return r'r = {} mm'.format(size * 1e4)