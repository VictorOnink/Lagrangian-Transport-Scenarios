import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string
import cmocean.cm as cmo


class SizeTransport_VerticalProfile:
    def __init__(self, scenario, figure_direc, size_list, time_selection):
        # Figure Parameters
        self.fig_size = (16, 10)
        self.fig_shape = (2, 2)
        self.x_label = 'Number of Particles'
        self.y_label = 'Depth (m)'
        self.ax_ticklabel_size = 12
        self.ax_label_size = 14
        self.legend_size = 11
        self.xmin, self.xmax = 1e-10, 1e0
        self.ymin, self.ymax = -200, 0
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.number_of_plots = 4

        # Data parameters
        self.output_direc = figure_direc + 'vertical_profile/'
        self.data_direc = utils.get_output_directory(
            server=settings.SERVER) + 'concentrations/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'vertical_concentration'
        # Simulation parameters
        self.scenario = scenario
        self.size_list = size_list
        self.time_selection = time_selection
        self.rho = 920
        self.tau = 0.0

    def plot(self):
        # Loading the data
        output_dict = {}
        for size in self.size_list:
            data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix, data_direc=self.data_direc,
                                                       size=size, rho=self.rho, tau=self.tau)
            output_dict[size] = data_dict[utils.analysis_simulation_year_key(self.time_selection)]
        depth_bins = -0.5 * (data_dict['depth'][1:] + data_dict['depth'][:-1])

        # Normalizing by the counts
        for size in self.size_list:
            for month in range(0, 12):
                output_dict[size][month]['concentration'] /= output_dict[size][month]['counts']

        # Creating the figure
        # Creating the figure
        ax = vUtils.base_figure(fig_size=self.fig_size, ax_range=self.ax_range, x_label=self.x_label,
                                y_label=self.y_label, ax_ticklabel_size=self.ax_ticklabel_size,
                                ax_label_size=self.ax_label_size, shape=self.fig_shape, plot_num=self.number_of_plots,
                                log_yscale=False, log_xscale=True, all_x_labels=True, all_y_labels=True,
                                legend_axis=True, width_ratios=[1, 1, 0.5])

        # Labelling the subfigures
        for index_ax in range(self.number_of_plots):
            ax[index_ax].set_title(subfigure_title(index_ax, self.time_selection), fontsize=self.ax_label_size)

        # Adding in a legend
        cmap_list, label_list = [], []
        for size_class in range(self.size_classes):
            cmap_list.append(vUtils.discrete_color_from_cmap(size_class, subdivisions=self.size_list.__len__()))
            label_list.append(legend_label(size_class))
        size_colors = [plt.plot([], [], c=cmap_list[i], label=label_list[i], linestyle='-')[0] for i in
                       range(cmap_list.__len__())]
        ax[-1].legend(handles=size_colors, fontsize=self.legend_size)
        ax[-1].axis('off')

        file_name = self.output_direc + 'SizeTransport_vertical_profile_{}.png'.format(self.time_selection)
        plt.savefig(file_name, bbox_inches='tight')


def legend_label(size):
    return r'r = {:.3f} mm'.format(size * 1e3)


def subfigure_title(index, simulation_year):
    alphabet = string.ascii_lowercase
    return '({}) 1-{}-{}'.format(alphabet[index], index * 3 + 1, settings.STARTYEAR + simulation_year)
