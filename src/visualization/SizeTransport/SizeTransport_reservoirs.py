import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import string
from datetime import datetime, timedelta


class SizeTransport_reservoirs:
    def __init__(self, scenario, figure_direc, size_list, rho_list=[920, 980]):
        utils.print_statement('Creating the SizeTransport reservoirs figure', to_print=True)
        # Simulation parameters
        self.scenario = scenario
        self.rho_list = rho_list
        self.size_list = size_list
        self.tau = 0.0
        # Data parameters
        self.output_direc = figure_direc + 'timeseries/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'timeseries'
        # Figure parameters
        self.figure_size = (12, 12)
        self.figure_shape = (1, 1)
        self.ax_label_size = 12
        self.ax_ticklabel_size = 12
        self.xmax, self.xmin = 1e-3, 1e1
        self.ymax, self.ymin = 100, 0
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.beach_state_list = ['beach', 'adrift']
        self.y_label = 'Fraction of Total (%)'
        self.x_label = 'Size (mm)'
        self.rho_marker_dict = {920: 'o', 980: 's'}
        self.state_color = {'beach': 'r', 'adrift': 'b'}
        self.number_of_plots = 1

    def plot(self):
        # Loading the data
        timeseries_dict = dict.fromkeys(self.rho_list)
        for rho in self.rho_list:
            timeseries_dict[rho] = {}
            for size in self.size_list:
                timeseries_dict[rho][size] = {}
                data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix, data_direc=self.data_direc,
                                                           size=size, rho=rho, tau=self.tau)
                for beach_state in self.beach_state_list:
                    timeseries_dict[rho][size][beach_state] = data_dict[beach_state]
                timeseries_dict[rho][size]['total_divide'] = data_dict['total'][0]


        # Normalizing all the particle counts with the total number of particles, and then multiplying by 100 to get a
        # percentage
        for rho in self.rho_list:
            for size in self.size_list:
                for beach_state in self.beach_state_list:
                    timeseries_dict[rho][size][beach_state] /= timeseries_dict[rho][size]['total_divide']
                    timeseries_dict[rho][size][beach_state] *= 100.
                    timeseries_dict[rho][size][beach_state] = timeseries_dict[rho][size][beach_state][-1]

        # Creating the figure
        ax = vUtils.base_figure(fig_size=self.figure_size, ax_range=self.ax_range, y_label=self.y_label,
                                x_label=self.x_label, ax_label_size=self.ax_label_size,
                                ax_ticklabel_size=self.ax_ticklabel_size, shape=self.figure_shape,
                                plot_num=self.number_of_plots, legend_axis=True, log_yscale=False, log_xscale=True,
                                all_x_labels=True)

        # Now, adding in the actual data
        for rho in self.rho_list:
            for index_size, size in enumerate(self.size_list):
                for index_beach, beach_state in enumerate(self.beach_state_list):
                    line_color = vUtils.discrete_color_from_cmap(index_size, subdivisions=len(self.size_list))
                    ax[index_beach].plot(size * 1e3, timeseries_dict[rho][size][beach_state],
                                         linestyle=None, marker=self.rho_marker_dict[rho], facecolor=None,
                                         color=line_color, label=size_label(size))
        # Creating a legend
        # rho_lines = [plt.plot([], [], c='k', label=r'$\rho=$' + str(rho) + r' kg m$^{-3}$', linestyle=self.rho_line_dict[rho])[0]
        #              for rho in self.rho_list]
        # size_number = self.size_list.__len__()
        # size_colors = [plt.plot([], [], c=vUtils.discrete_color_from_cmap(index_size, subdivisions=size_number),
        #                         label=size_label(size), linestyle='-')[0] for index_size, size in enumerate(self.size_list)]
        # ax[-1].legend(handles=rho_lines + size_colors, fontsize=self.ax_label_size, loc='upper right')
        # ax[-1].axis('off')

        file_name = self.output_direc + 'SizeTransport_reservoirs.jpg'
        plt.savefig(file_name, bbox_inches='tight')


def size_label(size):
    return r'r = {:.3f} mm'.format(size * 1e3)


