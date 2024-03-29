import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import string
from Analysis.FragmentationKaandorp_boxmodel import FragmentationKaandorp_box_model
from copy import deepcopy
import numpy as np


class FragmentationKaandorpPartial_boxmodel_ocean:
    def __init__(self, figure_direc, shore_time=26, lambda_frag_list=[388, 35000], rho=920, ocean_frag=True,
                 sim_length=10, size_class_number=10):
        # Data parameters
        self.output_direc = figure_direc + 'size_distribution/'
        self.data_direc = settings.DATA_OUTPUT_DIREC + 'size_distribution/FragmentationKaandorpPartial/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'size_distribution'
        # Simulation parameters
        self.shore_time = shore_time
        self.lambda_frag_list = lambda_frag_list
        self.rho = rho
        self.sim_length = sim_length
        self.ocean_frag = ocean_frag
        self.reservoir_list = ['ocean', 'coastal', 'beach']
        self.lambda_fO_list = np.array([10000, 1000, 100, 10, 5, 1])
        self.size_class_number = size_class_number
        # Figure parameters
        self.fig_size = (17, 10)
        self.fig_shape = (3, 2)
        self.x_label = 'Size (mm)'
        self.y_label = r'Normalized Particle Number (n mm$^{-1}$)'
        self.twiny_label = r'Normalized Particle Mass (g mm$^{-1}$)'
        self.ax_ticklabel_size = 12
        self.ax_label_size = 14
        self.legend_size = 12
        self.xmin, self.xmax = 1e-3, 1e1
        # self.ymin, self.ymax = {388: (1e2, 1e11), 35000: (1e3, 1e8)}[self.lambda_frag]
        self.ymin, self.ymax = 1e2, 1e12
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        # self.twin_ymin, self.twin_ymax = {388: (1e2, 1e5), 35000: (1e-1, 1e6)}[self.lambda_frag]
        self.twin_ymin, self.twin_ymax = 1e-2, 2e5
        self.twin_ax_range = self.xmax, self.xmin, self.twin_ymax, self.twin_ymin
        self.number_of_plots = self.fig_shape[0] * self.fig_shape[1]
        self.cmap = 'viridis'
        self.line_style = {388: '-', 35000: '--'}
        self.line_width = 3

    def plot(self):
        # Getting the sizes of the size classes, and we convert from meters to mm
        size_classes = utils.size_range(size_class_number=self.size_class_number, units='mm')

        # Loading the box model data, first the base model without ocean fragmentation
        base_mass, base_number = {}, {}
        for lambda_frag in self.lambda_frag_list:
            base_box_data = FragmentationKaandorp_box_model(sim_length=self.sim_length, lambda_f=lambda_frag,
                                                            size_classes=self.size_class_number,
                                                            steady_state=True).load_box_model()
            base_mass[lambda_frag], base_number[lambda_frag] = base_box_data['mass'], base_box_data['number']

        # Next the ocean fragmentation cases
        lambda_fO_mass, lambda_fO_number = {}, {}
        for lambda_frag in self.lambda_frag_list:
            lambda_fO_mass[lambda_frag], lambda_fO_number[lambda_frag] = {}, {}
            for lambda_fO in self.lambda_fO_list:
                box_model_data = FragmentationKaandorp_box_model(sim_length=self.sim_length, lambda_f=lambda_frag,
                                                                 size_classes=self.size_class_number, ocean_frag=True,
                                                                 lambda_f_ocean=lambda_frag * lambda_fO,
                                                                 steady_state=True).load_box_model()
                lambda_fO_mass[lambda_frag][lambda_fO] = box_model_data['mass']
                lambda_fO_number[lambda_frag][lambda_fO] = box_model_data['number']

        # Creating the figure
        ax, twin_ax = vUtils.base_figure(fig_size=self.fig_size, ax_range=self.ax_range, x_label=self.x_label,
                                         y_label=self.y_label, ax_ticklabel_size=self.ax_ticklabel_size,
                                         ax_label_size=self.ax_label_size, shape=self.fig_shape,
                                         plot_num=self.number_of_plots, log_yscale=True, log_xscale=True,
                                         all_x_labels=True, all_y_labels=False, add_twinx=True,
                                         twinx_y_label=self.twiny_label, twinx_ax_range=self.twin_ax_range,
                                         log_twinxscale=True, legend_axis=True, width_ratios=[1, 1, 0.5])

        # Labelling the subfigures
        for index_ax in range(self.number_of_plots):
            ax[index_ax].set_title(self.subfigure_title(index_ax), fontsize=self.ax_label_size)

        # Plotting the model distributions from the box model
        for lambda_frag in self.lambda_frag_list:
            line = self.line_style[lambda_frag]
            for index_reservoir, reservoir in enumerate(self.reservoir_list):
                ax[2 * index_reservoir].plot(size_classes, base_number[lambda_frag][reservoir], linestyle=line, c='k',
                                             linewidth=self.line_width)
                twin_ax[2 * index_reservoir + 1].plot(size_classes, base_mass[lambda_frag][reservoir],
                                                      linestyle=line, c='k', linewidth=self.line_width)
                # for index_fO, lambda_fO in enumerate(self.lambda_fO_list):
                #     c = vUtils.discrete_color_from_cmap(index=index_fO, subdivisions=self.lambda_fO_list.size,
                #                                         cmap=self.cmap)
                #     ax[2 * index_reservoir].plot(size_classes, lambda_fO_number[lambda_frag][lambda_fO][reservoir],
                #                                  linestyle=line, c=c, linewidth=self.line_width)
                #     twin_ax[2 * index_reservoir + 1].plot(size_classes, lambda_fO_mass[lambda_frag][lambda_fO][reservoir],
                #                                           linestyle=line, c=c, linewidth=self.line_width)

        # Adding a legend
        line_base_1 = [plt.plot([], [], c='k', label=r"$\lambda_{f, B}$ = 388 days", linestyle='-',
                                linewidth=self.line_width)[0]]
        line_base_2 = [plt.plot([], [], c='k', label=r"$\lambda_{f, B}$ = 35000 days", linestyle='--',
                                linewidth=self.line_width)[0]]
        line_colors = [plt.plot([], [], c=vUtils.discrete_color_from_cmap(index_fO, subdivisions=self.lambda_fO_list.size, cmap=self.cmap),
                                label=label(fO), linestyle='-', linewidth=self.line_width)[0] for index_fO, fO in enumerate(self.lambda_fO_list)]
        ax[-1].legend(handles=line_base_1 + line_base_2 + line_colors, fontsize=self.legend_size, loc='upper right')

        # Saving the figure
        str_format = self.shore_time, self.size_class_number
        file_name = self.output_direc + 'boxmodel_ocean_frag-ST={}_size_classes={}.png'.format(*str_format)
        plt.savefig(file_name, bbox_inches='tight', dpi=300)

    @staticmethod
    def subfigure_title(index):
        alphabet = string.ascii_lowercase
        subtitle_list = ['Counts', 'Mass']
        reservoir_list = ['Open water', 'Coastal', 'Beach']
        return '({}) {} - {}'.format(alphabet[index], reservoir_list[index // 2], subtitle_list[index % 2])


def label(lambda_fO):
    return r"$\lambda_{f,O}$ = " + '{} '.format(lambda_fO) + r'$\times$ $\lambda_{f, B}$'