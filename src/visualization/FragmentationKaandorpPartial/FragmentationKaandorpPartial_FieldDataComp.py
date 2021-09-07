import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import string
import numpy as np


def FragmentationKaandorpPartial_FieldDataComp(figure_direc, scenario, shore_time, lambda_frag_list, density,
                                               fig_size=(10, 14), x_label='Size (mm)',
                                               y_label=r'Normalized Particle Number (n mm$^{-1}$)',
                                               ax_ticklabel_size=12, ax_label_size=14, legend_size=12):
    # Setting the folder within which we have the output, and where we have the saved data
    output_direc = figure_direc + 'size_distribution/'
    utils.check_direc_exist(output_direc)
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/FragmentationKaandorpPartial/'

    # Getting the sizes of the size classes, and we convert from meters to mm
    sizes_class = utils.size_range(size_class_number=settings.SIZE_CLASS_NUMBER, units='mm')

    # Loading the data
    prefix = 'size_distribution'
    data_dict = {}
    for lambda_frag in lambda_frag_list:
        data_dict[lambda_frag] = vUtils.FragmentationKaandorpPartial_load_data(scenario=scenario, prefix=prefix,
                                                                  data_direc=data_direc, shore_time=shore_time,
                                                                  lambda_frag=lambda_frag, rho=density, postprocess=True)
    time_index = data_dict[lambda_frag]['final_index']
    beach_state_list = ['adrift_open_surf', 'adrift_10km_surf', 'beach']
    field_dict = utils.load_obj(vUtils.FragmentationKaandorpPartial_fielddata_filename())

    # Creating the figure
    ax_range = 2e2, 1e-1, 2e5, 1e-4
    plot_num = 3
    ax = vUtils.base_figure(fig_size=fig_size, ax_range=ax_range, x_label=x_label, y_label=y_label,
                            ax_ticklabel_size=ax_ticklabel_size, ax_label_size=ax_label_size, shape=(3, 1),
                            plot_num=plot_num, log_yscale=True, log_xscale=True, all_x_labels=True,
                            all_y_labels=False)

    # Labelling the subfigures
    for index_ax in range(plot_num):
        ax[index_ax].set_title(subfigure_title(index_ax, beach_state_list[index_ax]), fontsize=ax_label_size)

    # Plotting the model distributions
    for ax_index, sub_ax in enumerate(ax):
        for lambda_index, lambda_frag in enumerate(lambda_frag_list):
            norm_factor = data_dict[lambda_frag][beach_state_list[ax_index]][time_index][0]
            sub_ax.plot(sizes_class, data_dict[lambda_frag][beach_state_list[ax_index]][time_index] / norm_factor,
                        linestyle='-',
                        color=vUtils.discrete_color_from_cmap(index=lambda_index, subdivisions=len(lambda_frag_list)),
                        label=r'$\lambda_f=$' + '{} days'.format(lambda_frag))

    # Adding the field data
    # First for the open ocean
    ax[0].plot(field_dict['Cozar']['bin_midpoint'], field_dict['Cozar']['pdf_counts'] / field_dict['Cozar']['pdf_counts'][14],
               marker='x', linestyle='-', color='tab:red', label='Cozar et al. (2015)')
    # Then for coastal waters
    ax[1].plot(field_dict['RuizOrejon']['bin_midpoint'],
               field_dict['RuizOrejon']['pdf_counts'] / field_dict['RuizOrejon']['pdf_counts'][6],
               marker='x', linestyle='-', color='tab:red', label=r'Ruiz-Orej\`on et al. (2018)')
    # and finally for on the beach
    ax[2].plot(field_dict['Fok']['bin_midpoint'], field_dict['Fok']['pdf_counts'] / field_dict['Fok']['pdf_counts'][5],
               marker='x', linestyle='-', color='tab:red', label='Fok et al. (2017)')
    ax[2].plot(field_dict['Constant1']['bin_midpoint'], field_dict['Constant1']['pdf_counts'] / field_dict['Constant1']['pdf_counts'][-2],
               marker='x', linestyle='-', color='tab:blue', label='Constant et al. (2019), site 1')
    ax[2].plot(field_dict['Constant2']['bin_midpoint'], field_dict['Constant2']['pdf_counts'] / field_dict['Constant2']['pdf_counts'][-2],
               marker='x', linestyle='-', color='tab:orange', label='Constant et al. (2019), site 2')

    # Adding legends
    for sub_ax in ax:
        sub_ax.legend(fontsize=legend_size, loc='upper right')

    # Saving the figure
    file_name = output_direc + 'SizeSpectrumFieldData-ST={}-rho={}-lamf={}.png'.format(shore_time, density, lambda_frag)
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
    title_dict = {'adrift_open_surf': 'Adrift - open ocean', 'adrift_10km_surf': 'Adrift - coastal', 'beach': 'Beach'}
    return '({}) {}'.format(alphabet[index], title_dict[beach_state])
