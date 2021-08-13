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
                                     fig_size=(18, 8), ax_ticklabel_size=12, ax_label_size=14, legend_size=12,
                                     y_label='Distance (km)', x_label='Time'):
    # Setting the folder within which we have the output, and where we have the saved data
    output_direc = figure_direc + 'separation_distance/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'separation_distance/{}/'.format('SizeTransport')
    utils.check_direc_exist(output_direc)

    # Loading the separation distances for the particle size we've selected
    prefix = 'separation_distance'
    data_dict = vUtils.SizeTransport_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                               size=size_selection, rho=rho_selection, tau=tau_selection)

    # Getting the datetime objects for all of the time arrays
    startdate = datetime(settings.START_YEAR, 1, 1, 12, 0)
    time_list = []
    for t in data_dict['TIME']:
        time_list.append(startdate + timedelta(seconds=t))

    # Creating the figure
    ax_range = datetime(settings.START_YEAR + 1, 1, 1), datetime(settings.START_YEAR, 1, 1), 1000, 1
    ax = vUtils.base_figure(fig_size=fig_size, ax_range=ax_range, x_label=x_label, y_label=y_label,
                            ax_ticklabel_size=ax_ticklabel_size, ax_label_size=ax_label_size, shape=(1, 2),
                            plot_num=2, legend_axis=True, log_yscale=True, x_time_axis=True,
                            width_ratios=[1, 1, 0.3], all_x_labels=True)

    # Plotting the figure
    for index_size, size in enumerate(size_list):
        if size != size_selection:
            size_key = utils.init_size_key(size=size)
            color = vUtils.discrete_color_from_cmap(index=index_size, subdivisions=15)
            ax[0].plot(time_list, data_dict['MEAN'][size_key], color=color, label=size_label(size))
            ax[1].plot(time_list, data_dict['MEDIAN'][size_key], color=color)

    # adding a legend
    handles, labels = ax[0].get_legend_handles_labels()
    ax[-1].legend(handles=handles, labels=labels, fontsize=legend_size)

    # Adding subplot titles
    ax[0].set_title('(a) Mean separation to r = {:.3f} mm'.format(size * 1e3))
    ax[1].set_title('(b) Median separation to r = {:.3f} mm'.format(size * 1e3))

    # Saving the figure
    file_name = output_direc + 'Separation_distance_size={}.png'.format(size_selection)
    plt.savefig(file_name, bbox_inches='tight')


def size_label(size):
    return r'r = {:.3f} mm'.format(size * 1e3)