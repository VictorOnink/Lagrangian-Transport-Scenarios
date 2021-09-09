import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import string
import numpy as np


def FragmentationKaandorpPartial_SizeSpectrumBeach(figure_direc, scenario, shore_time, lambda_frag, rho,
                                                   fig_size=(10, 14), x_label='Size (m)', y_label='Number of Particles',
                                                   ax_ticklabel_size=12, ax_label_size=14, legend_size=14):
    # Setting the folder within which we have the output, and where we have the saved data
    output_direc = figure_direc + 'size_distribution/'
    utils.check_direc_exist(output_direc)
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/FragmentationKaandorpPartial/'

    # Loading in the data
    prefix = 'size_distribution'
    data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=scenario, prefix=prefix,
                                                              data_direc=data_direc, shore_time=shore_time,
                                                              lambda_frag=lambda_frag, rho=rho, postprocess=True)
    time_index = data_dict['final_index']
    size_bins = data_dict['size_bins'][:-1]
    beach_state_list = ['total', 'beach', 'seabed', 'adrift', 'adrift_2m']

    # Creating the figure
    ax_range = 1e-2, 1e-5, 1e5, 1e0
    plot_num = 5
    ax = vUtils.base_figure(fig_size=fig_size, ax_range=ax_range, x_label=x_label, y_label=y_label,
                            ax_ticklabel_size=ax_ticklabel_size, ax_label_size=ax_label_size, shape=(5, 1),
                            plot_num=plot_num, log_yscale=True, log_xscale=True, all_x_labels=True,
                            all_y_labels=False)

    # Labelling the subfigures
    for index_ax in range(plot_num):
        ax[index_ax].set_title(subfigure_title(index_ax, beach_state_list[index_ax]), fontsize=ax_label_size)

    # Plotting the size distributions
    for index_ax, beach_state in enumerate(beach_state_list):
        ax[index_ax].plot(size_bins, data_dict[beach_state][time_index], linestyle='-', color='r')

    file_name = output_direc + 'SizeSpectrumBeach-ST={}-rho={}-lamf={}.png'.format(shore_time, rho, lambda_frag)
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
