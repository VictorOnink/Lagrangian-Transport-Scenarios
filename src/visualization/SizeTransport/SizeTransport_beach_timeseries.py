import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import string
from datetime import datetime, timedelta


class SizeTransport_beach_timeseries:
    def __init__(self, scenario, figure_direc, size_list, simulation_years, rho=920, without_seabed=True):
        # Simulation parameters
        self.scenario = scenario
        self.rho = rho
        self.size_list = size_list
        self.tau = 0.0
        self.simulation_years = simulation_years
        self.without_seabed = without_seabed
        # Data parameters
        self.output_direc = figure_direc + 'timeseries/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'timeseries'
        # Figure parameters
        self.figure_size = (12, 10)
        self.figure_shape = (2, 1)
        self.ax_label_size = 12
        self.ax_ticklabel_size = 12
        self.xmax, self.xmin = datetime(2010 + self.simulation_years, 1, 1), datetime(2010, 1, 1)
        self.ymax, self.ymin = 100, 0
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.beach_state_list = ['beach', 'adrift', 'seabed']
        self.y_label = 'Fraction of Total (%)'
        self.x_label = 'Time (yr)'
        if self.without_seabed:
            _ = self.beach_state_list.pop(-1)
        self.number_of_plots = self.beach_state_list.__len__()

    def plot(self):
        # Loading the data
        timeseries_dict = {}
        for size in self.size_list:
            timeseries_dict[size] = {}
            data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix, data_direc=self.data_direc,
                                                       size=size, rho=self.rho, tau=self.tau)
            for beach_state in self.beach_state_list:
                timeseries_dict[size][beach_state] = data_dict[beach_state]
            timeseries_dict[size]['time_raw'] = data_dict['time']
            timeseries_dict[size]['total_divide'] = data_dict['total'][0]
        time_list = create_time_list(data_dict)

        # Normalizing all the particle counts with the total number of particles, and then multiplying by 100 to get a
        # percentage
        for size in self.size_list:
            for beach_state in self.beach_state_list:
                timeseries_dict[size][beach_state] /= timeseries_dict[size]['total_divide']
                timeseries_dict[size][beach_state] *= 100.

        # Creating the figure
        ax = vUtils.base_figure(fig_size=self.figure_size, ax_range=self.ax_range, y_label=self.y_label,
                                x_label=self.x_label, ax_label_size=self.ax_label_size,
                                ax_ticklabel_size=self.ax_ticklabel_size, shape=self.figure_shape,
                                plot_num=self.number_of_plots, legend_axis=True, log_yscale=False, x_time_axis=True,
                                width_ratios=[1, 0.3], all_x_labels=True)

        # Adding the subfigure titles
        for index, beach_state in enumerate(self.beach_state_list):
            ax[index].set_title(subfigure_title(index, beach_state), fontsize=self.ax_label_size)

        # Now, adding in the actual data
        for index_size, size in enumerate(self.size_list):
            for index_beach, beach_state in enumerate(self.beach_state_list):
                line_color = vUtils.discrete_color_from_cmap(index_size, subdivisions=len(self.size_list))
                ax[index_beach].plot(time_list, timeseries_dict[size][beach_state], linestyle='-',
                                     color=line_color, label=size_label(size))
        # Creating a legend
        size_number = self.size_list.__len__()
        size_colors = [plt.plot([], [], c=vUtils.discrete_color_from_cmap(index_size, subdivisions=size_number),
                                label=size_label(size), linestyle='-')[0] for index_size, size in enumerate(self.size_list)]
        ax[-1].legend(handles= size_colors, fontsize=self.ax_label_size, loc='upper right')
        ax[-1].axis('off')

        file_name = self.output_direc + 'SizeTransport_beach_state_timeseries.jpg'
        plt.savefig(file_name, bbox_inches='tight')


def subfigure_title(index, beach_state):
    return '({}) {}'.format(string.ascii_lowercase[index], beach_state)


def size_label(size):
    return r'r = {:.3f} mm'.format(size * 1e3)


def create_time_list(data_dict):
    time_list = []
    for time in data_dict['time']:
        time_list.append(datetime(settings.STARTYEAR, 1, 1, 12) + timedelta(seconds=time))
    return time_list


