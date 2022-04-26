import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
from datetime import datetime, timedelta


class SizeTransport_SeparationDistance:
    def __init__(self, scenario, figure_direc, size_selection, size_list, rho=920, tau=0.0, ):
        # Simulation parameters
        self.scenario = scenario
        self.rho = rho
        self.size_list = size_list
        self.size_selection = size_selection
        self.tau = tau
        # Data parameters
        self.output_direc = figure_direc + 'separation_distance/'
        self.data_direc = settings.DATA_OUTPUT_DIREC + 'separation_distance/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'separation_distance'
        # Figure parameters
        self.figure_size = (18, 8)
        self.figure_shape = (1, 2)
        self.ax_label_size = 14
        self.ax_ticklabel_size = 12
        self.legend_size=12
        self.xmax, self.xmin = datetime(2010 + 1, 1, 1), datetime(2010, 1, 1)
        self.ymax, self.ymin = 1000, 1
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.y_label = 'Distance (km)'
        self.x_label = 'Time (yr)'
        self.number_of_plots = 2

    def plot(self):
        # Loading the separation distances for the particle size we've selected
        data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix, data_direc=self.data_direc,
                                                   size=self.size_selection, rho=self.rho, tau=self.tau)

        # Getting the datetime objects for all of the time arrays
        reference_data = datetime(settings.STARTYEAR, 1, 1, 12, 0)
        time_list = []
        for t in data_dict['TIME']:
            time_list.append(reference_data + timedelta(seconds=t))

        # Creating the figure
        ax = vUtils.base_figure(fig_size=self.figure_size, ax_range=self.ax_range, x_label=self.x_label, y_label=self.y_label,
                                ax_ticklabel_size=self.ax_ticklabel_size, ax_label_size=self.ax_label_size, shape=self.figure_shape,
                                plot_num=self.number_of_plots, legend_axis=True, log_yscale=True, x_time_axis=True,
                                width_ratios=[1, 1, 0.3], all_x_labels=True)

        # Plotting the figure
        for index_size, size in enumerate(self.size_list):
            if size != self.size_selection:
                size_key = utils.init_size_key(size=size)
                color = vUtils.discrete_color_from_cmap(index=index_size, subdivisions=15)
                ax[0].plot(time_list, data_dict['MEAN'][size_key], color=color, label=size_label(size))
                ax[1].plot(time_list, data_dict['MEDIAN'][size_key], color=color)

        # adding a legend
        handles, labels = ax[0].get_legend_handles_labels()
        ax[-1].legend(handles=handles, labels=labels, fontsize=self.legend_size)

        # Adding subplot titles
        ax[0].set_title('(a) Mean separation to r = {:.3f} mm'.format(self.size_selection * 1e3))
        ax[1].set_title('(b) Median separation to r = {:.3f} mm'.format(self.size_selection * 1e3))

        # Saving the figure
        file_name = self.output_direc + 'Separation_distance_size={:.3f}.png'.format(self.size_selection * 1e3)
        plt.savefig(file_name, bbox_inches='tight')


def size_label(size):
    return r'r = {:.3f} mm'.format(size * 1e3)