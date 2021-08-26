import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import string
import numpy as np


def FragmentationKaandorpPartial_SizeSpectrumTime(figure_direc, scenario, shore_time, lambda_frag_list, density,
                                                  reservoir='total',
                                                  fig_size=(18, 12), x_label='Size (m)', y_label='Number of Particles',
                                                  ax_ticklabel_size=12, ax_label_size=14, legend_size=14):
    # Setting the folder within which we have the output, and where we have the saved data
    output_direc = figure_direc + 'size_distribution/'
    utils.check_direc_exist(output_direc)
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/FragmentationKaandorpPartial/'

    # Loading in the data
    prefix = 'size_distribution'
    size_dict = {}
    for lambda_frag in lambda_frag_list:
        data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=scenario, prefix=prefix,
                                                                  data_direc=data_direc, shore_time=shore_time,
                                                                  lambda_frag=lambda_frag, rho=density)
        size_dict[lambda_frag] = {}
        for month in data_dict[reservoir].keys():
            size_dict[lambda_frag][month] = data_dict[reservoir][month]

    size_bins = data_dict['size_bins'][:-1]

    # Creating the figure
    ax_range = 1e-2, 1e-5, 1e7, 1e-2
    plot_num = 6
    ax = vUtils.base_figure(fig_size=fig_size, ax_range=ax_range, x_label=x_label, y_label=y_label,
                            ax_ticklabel_size=ax_ticklabel_size, ax_label_size=ax_label_size, shape=(2, 3),
                            plot_num=plot_num, legend_axis=True, log_yscale=True, log_xscale=True, x_time_axis=False,
                            width_ratios=[1, 1, 1, 0.3], all_x_labels=True)

    # Labelling the subfigures
    for index_ax in range(plot_num):
        ax[index_ax].set_title(subfigure_title(index_ax, lambda_frag_list), fontsize=ax_label_size)

    # Plotting the size distributions one figure (so fragmentation timescale) at a time
    month_step = 12
    for index_ax in range(plot_num - 1, -1, -1):
        for index_month, month in enumerate(list(size_dict[lambda_frag_list[index_ax]].keys())[::month_step]):
            ax[index_ax].plot(size_bins, size_dict[lambda_frag_list[index_ax]][month], linestyle='-',
                              color=vUtils.discrete_color_from_cmap(index_month, 12 // month_step))

    # Adding in a legend
    size_colors = [plt.plot([], [], c=vUtils.discrete_color_from_cmap(month, subdivisions=12 // month_step),
                            label='t = {} months'.format(month * month_step), linestyle='-')[0] for month in range(12 // month_step)]

    ax[-1].legend(handles=size_colors, fontsize=legend_size)
    ax[-1].axis('off')

    file_name = output_direc + 'SizeSpectrumTimePartial_{}-ST={}-rho={}.png'.format(reservoir, shore_time, density)
    plt.savefig(file_name, bbox_inches='tight')


def subfigure_title(index, lambda_frag_list):
    """
    setting the title of the subfigure
    :param index:
    :param size:
    :param rho:
    :return:
    """
    alphabet = string.ascii_lowercase
    return '({}) '.format(alphabet[index]) + r'$\lambda_f$ = ' + '{} days'.format(lambda_frag_list[index])
