import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
from copy import deepcopy
import string
import numpy as np
import cmocean.cm as cmo


def FragmentationKaandorpPartial_vertical_profile(figure_direc, scenario, shore_time, lambda_frag, rho, simulation_year,
                                                  fig_size=(16, 10), x_label='Number of Particles', y_label='Depth (m)',
                                                  ax_ticklabel_size=12, ax_label_size=14, legend_size=11):
    # Setting the folder within which we have the output, and where we have the saved data
    output_direc = figure_direc + 'vertical_profile/'
    utils.check_direc_exist(output_direc)
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/FragmentationKaandorpPartial/'

    # Loading in the data
    prefix = 'vertical_concentration'
    data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=scenario, prefix=prefix,
                                                              data_direc=data_direc, shore_time=shore_time,
                                                              lambda_frag=lambda_frag, rho=rho, postprocess=True)
    fine_depth_bins = data_dict['depth']
    year_dict = deepcopy(data_dict[utils.analysis_simulation_year_key(simulation_year)])

    # Next, we are going to average out the bins, because they are currently too small
    data_dict = {}
    depth_bins = np.arange(start=np.min(fine_depth_bins), stop=np.max(fine_depth_bins), step=0.5)
    print(depth_bins)
    for month in year_dict.keys():
        data_dict[month] = {'concentration': None, 'counts': 0}
        for size_class in year_dict[month].keys():
            data_dict[month][size_class]['concentration'] = np.histogram(year_dict[month][size_class]['concentration'],
                                                                         bins=depth_bins)
            data_dict[month][size_class]['counts'] = year_dict[month][size_class]['counts']
    depth_bins = -0.5 * (depth_bins[1:] + depth_bins[:-1])

    # Creating the figure
    ax_range = 1e0, 1e-10, 0, -25
    plot_num = 4
    ax = vUtils.base_figure(fig_size=fig_size, ax_range=ax_range, x_label=x_label, y_label=y_label,
                            ax_ticklabel_size=ax_ticklabel_size, ax_label_size=ax_label_size, shape=(2, 2),
                            plot_num=plot_num, log_yscale=False, log_xscale=True, all_x_labels=True, all_y_labels=True,
                            legend_axis=True, width_ratios=[1, 1, 0.5])

    # Labelling the subfigures
    for index_ax in range(plot_num):
        ax[index_ax].set_title(subfigure_title(index_ax, simulation_year), fontsize=ax_label_size)

    # Adding in a legend
    cmap_list, label_list, line_list = ['k'], ['Total'], ['--']
    for size_class in range(settings.SIZE_CLASS_NUMBER):
        cmap_list.append(vUtils.discrete_color_from_cmap(size_class, subdivisions=settings.SIZE_CLASS_NUMBER, cmap='viridis'))
        label_list.append('Size class {}, d = {:.3f} mm'.format(size_class, utils.size_range(single_size_class=size_class,
                                                                                             units='mm')))
        line_list.append('-')

    size_colors = [plt.plot([], [], c=cmap_list[ind], label=label_list[ind], linestyle=line_list[ind])[0] for ind in range(cmap_list.__len__())]
    ax[-1].legend(handles=size_colors, fontsize=legend_size)
    ax[-1].axis('off')

    # And finally, the actual plotting:
    for ind_month, month in enumerate(np.arange(0, 12, 3)):
        total_count = np.zeros(data_dict[month][0]['concentration'].shape)
        for size_class in range(settings.SIZE_CLASS_NUMBER):
            total_count += data_dict[month][size_class]['concentration']
            norm_conc = data_dict[month][size_class]['concentration'] / data_dict[month][size_class]['counts']
            ax[ind_month].plot(norm_conc, depth_bins, linestyle='-',
                               c=vUtils.discrete_color_from_cmap(size_class, subdivisions=settings.SIZE_CLASS_NUMBER,
                                                                 cmap='viridis'))
        # ax[ind_month].plot(total_count, depth_bins, linestyle='-', c='k', zorder=10000)

    # Saving the figure
    str_format = lambda_frag, shore_time, rho, simulation_year
    file_name = output_direc + 'VerticalProfile-lamf={}-ST={}-rho={}_simyear={}.png'.format(*str_format)
    plt.savefig(file_name, bbox_inches='tight')


def subfigure_title(index, simulation_year):
    """
    setting the title of the subfigure
    :param index:
    :param size:
    :param rho:
    :return:
    """
    alphabet = string.ascii_lowercase
    return '({}) 1-{}-{}'.format(alphabet[index], index * 3 + 1, settings.STARTYEAR + simulation_year)

