import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
from copy import deepcopy
import string
import numpy as np
import cmocean.cm as cmo


class FragmentationKaandorpPartial_vertical_profile:
    def __init__(self, figure_direc, scenario, shore_time, lambda_frag, rho, simulation_year):
        # Figure Parameters
        self.fig_size = (16, 10)
        self.fig_shape = (2, 2)
        self.x_label = 'Number of Particles'
        self.y_label = 'Depth (m)'
        self.ax_ticklabel_size = 12
        self.ax_label_size = 14
        self.legend_size = 11
        self.xmin, self.xmax = 1e-10, 1e0
        self.ymin, self.ymax = -25, 0
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.number_of_plots = 4
        # Data parameters
        self.output_direc = figure_direc + 'vertical_profile/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/FragmentationKaandorpPartial/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'vertical_concentration'
        # Simulation parameters
        self.scenario = scenario
        self.shore_time = shore_time
        self.lambda_frag = lambda_frag
        self.rho = rho
        self.simulation_year = simulation_year
        self.size_classes = settings.SIZE_CLASS_NUMBER

    def plot(self):
        # Loading the data
        data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=self.scenario, prefix=self.prefix,
                                                                  data_direc=self.data_direc,
                                                                  shore_time=self.shore_time,
                                                                  lambda_frag=self.lambda_frag,
                                                                  rho=self.rho, postprocess=True)
        depth_bins = -0.5 * (data_dict['depth'][1:] + data_dict['depth'][:-1])

        # Creating the figure
        ax = vUtils.base_figure(fig_size=self.fig_size, ax_range=self.ax_range, x_label=self.x_label,
                                y_label=self.y_label, ax_ticklabel_size=self.ax_ticklabel_size,
                                ax_label_size=self.ax_label_size, shape=self.fig_shape, plot_num=self.number_of_plots,
                                log_yscale=False, log_xscale=True, all_x_labels=True, all_y_labels=True,
                                legend_axis=True, width_ratios=[1, 1, 0.5])

        # Labelling the subfigures
        for index_ax in range(self.number_of_plots):
            ax[index_ax].set_title(subfigure_title(index_ax, self.simulation_year), fontsize=self.ax_label_size)

        # Adding in a legend
        cmap_list, label_list, line_list = ['k'], ['Total'], ['--']
        for size_class in range(self.size_classes):
            cmap_list.append(vUtils.discrete_color_from_cmap(size_class, subdivisions=self.size_classes, cmap='viridis'))
            label_list.append(legend_label(size_class))
            line_list.append('-')
        size_colors = [plt.plot([], [], c=cmap_list[i], label=label_list[i], linestyle=line_list[i])[0] for i in
                       range(cmap_list.__len__())]
        ax[-1].legend(handles=size_colors, fontsize=self.legend_size)
        ax[-1].axis('off')

        # And finally, the actual plotting:
        for ind_month, month in enumerate(np.arange(0, 12, 3)):
            total_count = np.zeros(data_dict[month][0]['concentration'].shape)
            for size_class in range(settings.SIZE_CLASS_NUMBER):
                total_count += data_dict[month][size_class]['concentration']
                norm_conc = data_dict[month][size_class]['concentration'] / data_dict[month][size_class]['counts']
                c = vUtils.discrete_color_from_cmap(size_class, subdivisions=settings.SIZE_CLASS_NUMBER, cmap='viridis')
                ax[ind_month].plot(norm_conc, depth_bins, linestyle='-', c=c)
            ax[ind_month].plot(total_count, depth_bins, linestyle='-', c='k', zorder=10000)

        # Saving the figure
        str_format = self.lambda_frag, self.shore_time, self.rho, self.simulation_year
        file_name = self.output_direc + 'VerticalProfile-lamf={}-ST={}-rho={}_simyear={}.png'.format(*str_format)
        plt.savefig(file_name, bbox_inches='tight')


def legend_label(size_class):
    str_format = size_class, utils.size_range(single_size_class=size_class, units='mm')
    return 'Size class {}, d = {:.3f} mm'.format(*str_format)


def subfigure_title(index, simulation_year):
    alphabet = string.ascii_lowercase
    return '({}) 1-{}-{}'.format(alphabet[index], index * 3 + 1, settings.STARTYEAR + simulation_year)

