import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import string
from datetime import datetime, timedelta


class FragmentationKaandorpPartial_timeseries:
    def __init__(self, scenario, figure_direc, shore_time, lambda_frag, rho, simulation_length, weight):
        # Simulation parameters
        self.scenario = scenario
        self.shore_time = shore_time
        self.lambda_frag = lambda_frag
        self.rho = rho
        self.simulation_length = simulation_length
        self.class_num = settings.SIZE_CLASS_NUMBER
        self.weight = weight
        # Data parameters
        self.output_direc = figure_direc + 'timeseries/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/FragmentationKaandorpPartial/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'timeseries'
        self.beach_state_list = ['beach', 'afloat', 'total']
        # Figure parameters
        self.figure_size = (20, 14)
        self.figure_shape = (self.beach_state_list.__len__(), 1)
        self.ax_label_size = 16
        self.ax_ticklabel_size = 14
        self.legend_size = 14
        self.y_label = 'Particle Number'
        self.x_label = 'Time'
        self.xmin, self.xmax = datetime(settings.STARTYEAR, 1, 1), \
                               datetime(settings.STARTYEAR + self.simulation_length, 1, 1)
        self.ymin, self.ymax = 1e0, 1e5
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.number_of_plots = self.beach_state_list.__len__()

    def plot(self):
        # Loading data
        timeseries_dict = {}
        for size_class in range(self.class_num):
            timeseries_dict[size_class] = {}
            data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=self.scenario, prefix=self.prefix,
                                                                      data_direc=self.data_direc,
                                                                      shore_time=self.shore_time,
                                                                      lambda_frag=self.lambda_frag, rho=self.rho,
                                                                      postprocess=True)
            for beach_state in self.beach_state_list:
                timeseries_dict[size_class][beach_state] = data_dict[beach_state][size_class][self.weight]

        # creating a time axis
        time_list = []
        for time in data_dict['time']:
            time_list.append(datetime(settings.STARTYEAR, 1, 1, 12) + timedelta(seconds=time))

        # Creating figure
        ax = vUtils.base_figure(fig_size=self.fig_size, ax_range=self.ax_range, y_label=self.y_label,
                                x_label=self.x_label, ax_label_size=self.ax_label_size,
                                ax_ticklabel_size=self.ax_ticklabel_size, shape=self.figure_shape,
                                plot_num=self.number_of_plots, legend_axis=True, log_yscale=True, x_time_axis=True,
                                width_ratios=[1, 0.3], all_x_labels=True)

        # Setting the subfigure titles
        for ax_index in range(self.number_of_plots):
            ax[ax_index].set_title(subfigure_title(ax_index, self.beach_state_list[ax_index]),
                                   fontsize=self.ax_label_size)

        # Creating a legend
        size_colors = [plt.plot([], [], c=vUtils.discrete_color_from_cmap(sc, subdivisions=self.class_num),
                                label=size_label(sc), linestyle='-')[0] for sc in range(self.class_num)]
        ax[-1].legend(handles=size_colors, fontsize=self.legend_size, loc='upper right')

        # Plotting the various fractions
        for size_class in range(self.class_num):
            for beach_index, beach_state in enumerate(self.beach_state_list):
                color_size = vUtils.discrete_color_from_cmap(size_class, subdivisions=self.class_num)
                ax[beach_index].plot(time_list, timeseries_dict[size_class][beach_state], linestyle='-',
                                     c=color_size)

        # Saving the figure
        plt.savefig(file_name(self.output_direc, self.shore_time, self.lambda_frag), bbox_inches='tight')


def file_name(output_direc, shore_time, lambda_frag):
    str_format = shore_time, lambda_frag
    return output_direc + 'FragmentationKaandorpPartial_beach_state_timeseries_ST={}_lamf={}.png'.format(*str_format)


def subfigure_title(index, beach_state):
    """
    setting the title of the subfigure
    :param index:
    :return:
    """
    alphabet = string.ascii_lowercase
    return '({}) {}'.format(alphabet[index], beach_state)


def size_label(size_class):
    particle_size = settings.INIT_SIZE * settings.P_FRAG ** size_class
    return 'Size class {}, d = {:.2f} mm'.format(size_class, particle_size * 1e3)