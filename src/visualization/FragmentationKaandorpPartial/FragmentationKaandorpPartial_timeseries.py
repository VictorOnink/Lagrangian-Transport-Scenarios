import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import numpy as np
import string
from datetime import datetime, timedelta
import matplotlib.dates as mdates


def FragmentationKaandorpPartial_timeseries(scenario, figure_direc, shore_time, lambda_frag, rho, simulation_length,
                                            fig_size=(12, 10), ax_label_size=14, ax_ticklabel_size=12,
                                            y_label='Particles', x_label='Time'):
    # Setting the folder within which we have the output, and where we have the saved timeslices
    output_direc = figure_direc + 'timeseries/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/FragmentationKaandorpPartial/'
    utils.check_direc_exist(output_direc)

    # Loading in the data
    prefix = 'timeseries'
    timeseries_dict = {}
    # beach_state_list = ['beach', 'afloat', 'seabed', 'removed', 'total']
    beach_state_list = ['beach', 'afloat', 'seabed', 'total']
    for size_class in range(settings.SIZE_CLASS_NUMBER):
        timeseries_dict[size_class] = {}
        data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                                                  shore_time=shore_time, lambda_frag=lambda_frag, rho=rho)
        for beach_state in beach_state_list:
            timeseries_dict[size_class][beach_state] = data_dict[beach_state][size_class]
        timeseries_dict[size_class]['time_raw'] = data_dict['time'][size_class]

    ax_range = datetime(settings.START_YEAR + simulation_length, 1, 1), datetime(settings.START_YEAR, 1, 1), 1e-2, 1e2
    figure_shape = (4, 1)
    ax = vUtils.base_figure(fig_size=fig_size, ax_range=ax_range, y_label=y_label, x_label=x_label,
                            ax_label_size=ax_label_size, ax_ticklabel_size=ax_ticklabel_size, shape=figure_shape,
                            plot_num=4, legend_axis=True, log_yscale=True, x_time_axis=True, width_ratios=[1, 0.2])

    # Creating a legend
    size_colors = [plt.plot([], [], c=vUtils.discrete_color_from_cmap(size_class, subdivisions=settings.SIZE_CLASS_NUMBER),
                            label=size_label(size_class), linestyle='-')[0] for size_class in range(settings.SIZE_CLASS_NUMBER)]
    ax[-1].legend(handles=size_colors, fontsize=ax_label_size, loc='upper right')

    # Saving the figure
    file_name = output_direc + 'FragmentationKaandorpPartial_beach_state_timeseries_ST={}_lamf={}.png'.format(shore_time,
                                                                                                              lambda_frag)
    plt.savefig(file_name, bbox_inches='tight')


def size_label(size_class):
    particle_size = settings.INIT_SIZE * settings.P_FRAG ** size_class
    return 'Size class {}, d = {:.3f} mm'.format(size_class, particle_size * 1e3)