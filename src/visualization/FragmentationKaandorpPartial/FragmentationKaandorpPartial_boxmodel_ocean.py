import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import string
from Analysis.FragmentationKaandorp_boxmodel import FragmentationKaandorp_box_model
from copy import deepcopy
import numpy as np


class FragmentationKaandorpPartial_boxmodel_ocean:
    def __init__(self, figure_direc, shore_time=20, lambda_frag=388, rho=920, ocean_frag=True, sim_length=10,
                 size_class_number=settings.SIZE_CLASS_NUMBER):
        # Data parameters
        self.output_direc = figure_direc + 'size_distribution/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/FragmentationKaandorpPartial/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'size_distribution'
        # Simulation parameters
        self.shore_time = shore_time
        self.lambda_frag = lambda_frag
        self.rho = rho
        self.class_num = settings.SIZE_CLASS_NUMBER
        self.sim_length = sim_length
        self.ocean_frag = ocean_frag
        self.reservoir_list = ['ocean', 'coastal', 'beach']
        self.lambda_fO_list = self.lambda_frag * np.array([1000000])
        self.size_class_number = size_class_number
        # Figure parameters
        self.fig_size = (16, 10)
        self.fig_shape = (3, 2)
        self.x_label = 'Size (mm)'
        self.y_label = r'Particle Number (n)'
        self.twiny_label = r'Particle Mass (g)'
        self.ax_ticklabel_size = 12
        self.ax_label_size = 14
        self.legend_size = 12
        self.xmin, self.xmax = 1e-3, 2e2
        self.ymin, self.ymax = 1e2, 1e6
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.twin_ymin, self.twin_ymax = 1e-4, 1e2
        self.twin_ax_range = self.xmax, self.xmin, self.twin_ymax, self.twin_ymin
        self.number_of_plots = self.fig_shape[0] * self.fig_shape[1]
        self.cmap = 'viridis'

    def plot(self):
        # Getting the sizes of the size classes, and we convert from meters to mm
        size_classes = utils.size_range(size_class_number=self.class_num, units='mm')

        # Loading the box model data, first the base model without ocean fragmentation
        base_box_data = FragmentationKaandorp_box_model(sim_length=self.sim_length, lambda_f=388,
                                                        size_classes=self.size_class_number).load_box_model()
        time_index = base_box_data['mass']['final_index']
        base_mass, base_number = base_box_data['mass'][time_index], base_box_data['number'][time_index]
        # Next the ocean fragmentation cases
        lambda_fO_mass, lambda_fO_number = {}, {}
        for lambda_fO in self.lambda_fO_list:
            box_model_data = FragmentationKaandorp_box_model(sim_length=self.sim_length, lambda_f=388,
                                                             size_classes=self.size_class_number, ocean_frag=True,
                                                             lambda_f_ocean=lambda_fO).load_box_model()
            lambda_fO_mass[lambda_fO], lambda_fO_number[lambda_fO] = box_model_data['mass'][time_index], \
                                                                     box_model_data['number'][time_index]

        # Creating the figure
        ax, twin_ax = vUtils.base_figure(fig_size=self.fig_size, ax_range=self.ax_range, x_label=self.x_label,
                                         y_label=self.y_label, ax_ticklabel_size=self.ax_ticklabel_size,
                                         ax_label_size=self.ax_label_size, shape=self.fig_shape,
                                         plot_num=self.number_of_plots, log_yscale=True, log_xscale=True,
                                         all_x_labels=True, all_y_labels=True, add_twinx=True,
                                         twinx_y_label=self.twiny_label, twinx_ax_range=self.twin_ax_range,
                                         log_twinxscale=True, legend_axis=True, width_ratios=[1, 1, 0.5])

        # Labelling the subfigures
        for index_ax in range(self.number_of_plots):
            ax[index_ax].set_title(self.subfigure_title(index_ax), fontsize=self.ax_label_size)

        # Plotting the model distributions from the box model
        for index_reservoir, reservoir in enumerate(self.reservoir_list):
            ax[2 * index_reservoir].plot(size_classes, base_number[reservoir], linestyle='-', c='k')
            twin_ax[2 * index_reservoir + 1].plot(size_classes, base_mass[reservoir], linestyle='-', c='k')
            for index_fO, lambda_fO in enumerate(self.lambda_fO_list):
                c = vUtils.discrete_color_from_cmap(index=index_fO, subdivisions=self.lambda_fO_list.size,
                                                    cmap=self.cmap)
                ax[2 * index_reservoir].plot(size_classes, lambda_fO_number[lambda_fO][reservoir], linestyle='-', c=c)
                twin_ax[2 * index_reservoir + 1].plot(size_classes, lambda_fO_mass[lambda_fO][reservoir], linestyle='-', c=c)

        # Adding a legend
        line_base = [plt.plot([], [], c='k', label='base', linestyle='-')[0]]
        line_colors = [plt.plot([], [], c=vUtils.discrete_color_from_cmap(index_fO, subdivisions=self.lambda_fO_list.size, cmap=self.cmap),
                                label=label(fO), linestyle='-')[0] for index_fO, fO in enumerate(self.lambda_fO_list)]
        ax[-1].legend(handles=line_base + line_colors, fontsize=self.legend_size, loc='upper right')

        # Saving the figure
        str_format = self.shore_time, self.lambda_frag, self.size_class_number
        file_name = self.output_direc + 'boxmodel_ocean_frag-ST={}-lambda_f={}_size_classes={}.png'.format(*str_format)
        plt.savefig(file_name, bbox_inches='tight', dpi=300)


    def subfigure_title(self, index):
        alphabet = string.ascii_lowercase
        subtitle_list = ['Counts', 'Mass']
        return '({}) {} - {}'.format(alphabet[index], self.reservoir_list[index // 2], subtitle_list[index % 2])


def label(lambda_fO):
    return r"$\lambda_{f,O}$ = " + '{:.3d} days'.format(lambda_fO)