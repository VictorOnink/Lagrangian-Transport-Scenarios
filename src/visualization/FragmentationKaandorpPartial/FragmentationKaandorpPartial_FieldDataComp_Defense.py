import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import string
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd


class FragmentationKaandorpPartial_FieldDataComp_Defense:
    def __init__(self, figure_direc, scenario, shore_time, lambda_frag_list, rho, sink=True, with_input=True,
                 input_list=['LebretonDivision']):
        # Data parameters
        self.output_direc = figure_direc + 'size_distribution/'
        self.data_direc = settings.DATA_OUTPUT_DIREC + 'size_distribution/FragmentationKaandorpPartial/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'size_distribution'
        self.beach_state_list = ['adrift_open_surf', 'adrift_10km_surf', 'beach']
        # Simulation parameters
        self.scenario = scenario
        self.shore_time = shore_time
        self.lambda_frag_list = lambda_frag_list
        self.rho = rho
        self.class_num = settings.SIZE_CLASS_NUMBER
        self.sink = sink
        if self.sink:
            self.count, self.mass = 'particle_number_sink', 'particle_mass_sink'
        else:
            self.count, self.mass = 'particle_number', 'particle_mass'
        self.input_list = input_list
        self.input_line_style = {'LebretonDivision': 'dotted', 'LebretonKaandorpInit': '-'}
        # Figure parameters
        self.fig_size = (14, 14)
        self.fig_shape = (2, 2)
        self.x_label = 'Size (mm)'
        self.y_label = r'Normalized Particle Number (n mm$^{-1}$)'
        self.twiny_label = r'Normalized Particle Mass (g mm$^{-1}$)'
        self.ax_ticklabel_size = 12
        self.ax_label_size = 14
        self.legend_size = 12
        self.xmin, self.xmax = 1e-1, 2e2
        self.ymin, self.ymax = 1e-1, 1e3
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.twin_ymin, self.twin_ymax = 1e-3, 1e1
        self.twin_ax_range = self.xmax, self.xmin, self.twin_ymax, self.twin_ymin
        self.number_of_plots = 3
        self.field_marker = 'X'
        self.field_line = '--'
        self.line_width = 3
        self.line_width_observations = 2
        self.marker_size = 8
        self.legend_loc = 'upper right'
        self.with_input = with_input

    def plot(self):
        # Getting the sizes of the size classes, and we convert from meters to mm
        size_classes = utils.size_range(size_class_number=self.class_num, units='mm')

        # Loading the data
        data_dict = {}
        for input_scenario in self.input_list:
            data_dict[input_scenario] = {}
            for lambda_frag in self.lambda_frag_list:
                data = vUtils.FragmentationKaandorpPartial_load_data(scenario=self.scenario, prefix=self.prefix,
                                                                     data_direc=self.data_direc, shore_time=self.shore_time,
                                                                     lambda_frag=lambda_frag, rho=self.rho,
                                                                     input=input_scenario, postprocess=True)
                data_dict[input_scenario][lambda_frag] = {}
                for beach_state in self.beach_state_list:
                    data_dict[input_scenario][lambda_frag][beach_state] = {}
                    data_dict[input_scenario][lambda_frag][beach_state][self.count] = data[beach_state][self.count]
                    data_dict[input_scenario][lambda_frag][beach_state][self.mass] = data[beach_state][self.mass]
        time_index = data['final_index'] // 2
        field_dict = utils.load_obj(vUtils.FragmentationKaandorpPartial_fielddata_filename())

        # Creating the figure
        ax = vUtils.base_figure(fig_size=self.fig_size, ax_range=self.ax_range, x_label=self.x_label,
                                y_label=self.y_label, ax_ticklabel_size=self.ax_ticklabel_size,
                                ax_label_size=self.ax_label_size, shape=self.fig_shape,
                                plot_num=self.number_of_plots, log_yscale=True, log_xscale=True,
                                all_x_labels=True, all_y_labels=False)

        # Labelling the subfigures
        for index_ax in range(self.number_of_plots):
            ax[index_ax].set_title(self.subfigure_title(index_ax),
                                   fontsize=self.ax_label_size)

        # Plotting the model distributions
        for ax_index, sub_ax in enumerate(ax):
            for input_scenario in self.input_list:
                for lambda_index, lambda_frag in enumerate(self.lambda_frag_list):
                    c = vUtils.discrete_color_from_cmap(index=lambda_index, subdivisions=self.lambda_frag_list.__len__())
                    linestyle = self.input_line_style[input_scenario]
                    bin_norm_data = data_dict[input_scenario][lambda_frag][self.beach_state_list[ax_index // 2]][self.count][time_index] / size_classes
                    norm_factor = bin_norm_data[0]
                    sub_ax.plot(size_classes, bin_norm_data / norm_factor, linestyle=linestyle, color=c,
                                linewidth=self.line_width)
            # Add shading to indicate the region where sampling by trawl net is more difficult
            if ax_index < 4:
                sub_ax.fill_betweenx(np.arange(1e-4, 1e4), self.xmin, 0.33, color='grey', alpha=0.2)

        # for lambda_index, lambda_frag in enumerate(self.lambda_frag_list):
        #     c = vUtils.discrete_color_from_cmap(index=lambda_index, subdivisions=self.lambda_frag_list.__len__())
        #     twin_ax[1].plot([], [], color=c, label=self.label(lambda_frag), linestyle='-', linewidth=self.line_width)
        # twin_ax[1].legend(fontsize=self.legend_size, loc=self.legend_loc)

        for lambda_index, lambda_frag in enumerate(self.lambda_frag_list):
            c = vUtils.discrete_color_from_cmap(index=lambda_index, subdivisions=self.lambda_frag_list.__len__())
            ax[0].plot([], [], color=c, label=self.label(lambda_frag), linestyle='-', linewidth=self.line_width)
        ax[0].legend(fontsize=self.legend_size, loc=self.legend_loc)

        # Adding input data if self.with_input == True
        if self.with_input:
            self.add_input_distribution(ax)

        # Saving the figure
        str_format = self.shore_time, self.rho, self.sink
        file_name = self.output_direc + 'SizeSpectrumFieldData-ST={}-rho={}-sink={}_DEFENSE.png'.format(*str_format)
        plt.savefig(file_name, bbox_inches='tight', dpi=400)
        plt.close()

    def add_input_distribution(self, ax):
        # The field data from Zeri et al. (2021) https://doi.org/10.3390/su13105328
        size_bins = np.arange(0, 5.2, 0.1)
        lengths = pd.read_excel(settings.DATA_INPUT_DIREC + "Zeri_Sustainability_MP_sizes.xls",
                                sheet_name='Combined')['LENGTH']
        field_data, _ = np.histogram(lengths, size_bins)
        size_bins = (size_bins[1:] + size_bins[:-1]) / 2
        # Creating the interpolation function from the bin sizes and the uncorrected field data
        interpolation_function = interp1d(size_bins, field_data)
        # Calculating the particle number at the sizes equivalent to the size class
        number_inter = interpolation_function(utils.size_range(units='mm', size_class_number=settings.SIZE_CLASS_NUMBER))
        number_norm = number_inter / utils.size_range(units='mm', size_class_number=settings.SIZE_CLASS_NUMBER)
        number_norm = number_norm / number_norm[0]
        # Plotting the normalized count and mass inputs
        ax[2].plot(utils.size_range(units='mm', size_class_number=settings.SIZE_CLASS_NUMBER),
                   number_norm, marker=self.field_marker, linestyle=self.field_line, color='tab:cyan',
                   label=r'Zeri et al. (2021) Input', linewidth=self.line_width, markersize=self.marker_size + 2)
        ax[2].legend(fontsize=self.legend_size, loc=self.legend_loc)

    def subfigure_title(self, index):
        alphabet = string.ascii_lowercase
        title_dict = {'adrift_open_surf': 'Adrift - Open water', 'adrift_10km_surf': 'Adrift - Coastal', 'beach': 'Beach'}
        return '({}) {}'.format(alphabet[index], title_dict[self.beach_state_list[index]])

    @staticmethod
    def label(lambda_frag):
        return r'$\lambda_f=$' + '{} days'.format(lambda_frag)
