import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import string
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd



class FragmentationKaandorpPartial_FieldDataComp:
    def __init__(self, figure_direc, scenario, shore_time, lambda_frag_list, rho, sink=True, with_input=True,
                 input_list=['LebretonDivision']):
        # Data parameters
        self.output_direc = figure_direc + 'size_distribution/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/FragmentationKaandorpPartial/'
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
        self.input_line_style = {'LebretonDivision': '-', 'LebretonKaandorpInit': 'dotted'}
        # Figure parameters
        self.fig_size = (14, 14)
        self.fig_shape = (self.beach_state_list.__len__(), 2)
        self.x_label = 'Size (mm)'
        self.y_label = r'Normalized Particle Number (n mm$^{-1}$)'
        self.twiny_label = r'Normalized Particle Mass (g mm$^{-1}$)'
        self.ax_ticklabel_size = 12
        self.ax_label_size = 14
        self.legend_size = 10
        self.xmin, self.xmax = 1e-1, 2e2
        self.ymin, self.ymax = 1e-4, 1e4
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.twin_ymin, self.twin_ymax = 1e-5, 1e1
        self.twin_ax_range = self.xmax, self.xmin, self.twin_ymax, self.twin_ymin
        self.number_of_plots = self.fig_shape[0] * self.fig_shape[1]
        self.field_marker = 'x'
        self.field_line = '--'
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
        ax, twin_ax = vUtils.base_figure(fig_size=self.fig_size, ax_range=self.ax_range, x_label=self.x_label,
                                         y_label=self.y_label, ax_ticklabel_size=self.ax_ticklabel_size,
                                         ax_label_size=self.ax_label_size, shape=self.fig_shape,
                                         plot_num=self.number_of_plots, log_yscale=True, log_xscale=True,
                                         all_x_labels=True, all_y_labels=False, add_twinx=True,
                                         twinx_y_label=self.twiny_label, twinx_ax_range=self.twin_ax_range,
                                         log_twinxscale=True)

        # Labelling the subfigures
        for index_ax in range(self.number_of_plots):
            ax[index_ax].set_title(subfigure_title(index_ax, self.beach_state_list[index_ax // 2]),
                                   fontsize=self.ax_label_size)

        # Plotting the model distributions
        for ax_index, sub_ax in enumerate(ax):
            for input_scenario in self.input_list:
                for lambda_index, lambda_frag in enumerate(self.lambda_frag_list):
                    c = vUtils.discrete_color_from_cmap(index=lambda_index, subdivisions=self.lambda_frag_list.__len__())
                    linestyle = self.input_line_style[input_scenario]
                    if ax_index % 2 == 0:
                        bin_norm_data = data_dict[input_scenario][lambda_frag][self.beach_state_list[ax_index // 2]][self.count][time_index] / size_classes
                        norm_factor = bin_norm_data[0]
                        sub_ax.plot(size_classes, bin_norm_data / norm_factor, linestyle=linestyle, color=c)
                    else:
                        bin_norm_data = data_dict[input_scenario][lambda_frag][self.beach_state_list[ax_index // 2]][self.mass][time_index] / size_classes
                        norm_factor = bin_norm_data[0]
                        twin_ax[ax_index].plot(size_classes, bin_norm_data / norm_factor, linestyle=linestyle, color=c)

        # Field data - open ocean
        norm_factor = field_dict['Cozar']['pdf_counts'][14]
        ax[0].plot(field_dict['Cozar']['bin_midpoint'], field_dict['Cozar']['pdf_counts'] / norm_factor,
                   marker=self.field_marker, linestyle=self.field_line, color='tab:red', label='Cozar et al. (2015)')
        ax[0].legend(fontsize=self.legend_size, loc=self.legend_loc)

        # Field data - coastal waters
        norm_factor = field_dict['RuizOrejon']['pdf_counts'][6]
        ax[2].plot(field_dict['RuizOrejon']['bin_midpoint'], field_dict['RuizOrejon']['pdf_counts'] / norm_factor,
                   marker=self.field_marker, linestyle=self.field_line, color='tab:red',
                   label=r'Ruiz-Orej$\`o$n et al. (2018)')
        ax[2].legend(fontsize=self.legend_size, loc=self.legend_loc)

        norm_factor = field_dict['Gundogdu2017']['pdf_counts'][14]
        ax[2].plot(field_dict['Gundogdu2017']['bin_midpoint'], field_dict['Gundogdu2017']['pdf_counts'] / norm_factor,
                   marker=self.field_marker, linestyle=self.field_line, color='tab:blue',
                   label=r'Gundogdu & Cem (2017)')
        ax[2].legend(fontsize=self.legend_size, loc=self.legend_loc)

        norm_factor = field_dict['DeHaanSanchez']['pdf_counts'][14]
        ax[2].plot(field_dict['DeHaanSanchez']['bin_midpoint'], field_dict['DeHaanSanchez']['pdf_counts'] / norm_factor,
                   marker=self.field_marker, linestyle=self.field_line, color='tab:orange',
                   label=r'de Haan et al. (submitted)')
        ax[2].legend(fontsize=self.legend_size, loc=self.legend_loc)

        # Field data - beach, microplastic counts
        norm_factor = field_dict['Fok']['pdf_counts'][5]
        ax[4].plot(field_dict['Fok']['bin_midpoint'], field_dict['Fok']['pdf_counts'] / norm_factor,
                   marker=self.field_marker, linestyle=self.field_line, color='tab:red', label='Fok et al. (2017)')
        norm_factor = field_dict['Constant1']['pdf_counts'][-2]
        ax[4].plot(field_dict['Constant1']['bin_midpoint'], field_dict['Constant1']['pdf_counts'] / norm_factor,
                   marker=self.field_marker, linestyle=self.field_line, color='tab:blue',
                   label='Constant et al. (2019), site 1')
        norm_factor = field_dict['Constant2']['pdf_counts'][-2]
        ax[4].plot(field_dict['Constant2']['bin_midpoint'], field_dict['Constant2']['pdf_counts'] / norm_factor,
                   marker=self.field_marker, linestyle=self.field_line, color='tab:orange',
                   label='Constant et al. (2019), site 2')
        norm_factor = field_dict['Merlino2020']['pdf_counts'][-1]
        ax[4].plot(field_dict['Merlino2020']['bin_midpoint'], field_dict['Merlino2020']['pdf_counts'] / norm_factor,
                   marker=self.field_marker, linestyle=self.field_line, color='tab:green',
                   label='Merlino et al. (2020)')
        ax[4].legend(fontsize=self.legend_size, loc=self.legend_loc)

        # Field data - beach, microplastic mass
        norm_factor = field_dict['Fok']['pdf_mass'][5]
        twin_ax[5].plot(field_dict['Fok']['bin_midpoint'], field_dict['Fok']['pdf_mass'] / norm_factor,
                        marker=self.field_marker, linestyle=self.field_line, color='tab:red', label='Fok et al. (2017)')
        twin_ax[5].legend(fontsize=self.legend_size, loc=self.legend_loc)

        # Adding input data if self.with_input == True
        if self.with_input:
            self.add_input_distribution(ax)

        # Adding a legend for the model line colors in the top right panel
        input_names = ["k = 0 input", r'Zeri et al. (2021)']
        for index_input, input_scenario in enumerate(self.input_list):
            twin_ax[1].plot([], [], color='k', label=input_names[index_input],
                            linestyle=self.input_line_style[input_scenario])
        for lambda_index, lambda_frag in enumerate(self.lambda_frag_list):
            c = vUtils.discrete_color_from_cmap(index=lambda_index, subdivisions=self.lambda_frag_list.__len__())
            twin_ax[1].plot([], [], color=c, label=label(lambda_frag), linestyle='-')
        twin_ax[1].legend(fontsize=self.legend_size, loc=self.legend_loc)

        # Saving the figure
        str_format = self.shore_time, self.rho, self.sink
        file_name = self.output_direc + 'SizeSpectrumFieldData-ST={}-rho={}-sink={}.png'.format(*str_format)
        plt.savefig(file_name, bbox_inches='tight', dpi=400)
        plt.close()

    def add_input_distribution(self, ax):
        # The field data from Zeri et al. (2021) https://doi.org/10.3390/su13105328
        size_bins = np.arange(0, 5.2, 0.1)
        lengths = pd.read_excel(settings.DATA_INPUT_DIR_SERVERS[settings.SERVER] + "Zeri_Sustainability_MP_sizes.xls",
                                sheet_name='Combined')['LENGTH']
        field_data, _ = np.histogram(lengths, size_bins)
        size_bins = (size_bins[1:] + size_bins[:-1]) / 2
        # Creating the interpolation function from the bin sizes and the uncorrected field data
        interpolation_function = interp1d(size_bins, field_data)
        # Calculating the particle number at the sizes equivalent to the size class
        number_inter = interpolation_function(utils.size_range(units='mm', size_class_number=settings.SIZE_CLASS_NUMBER))
        number_norm = number_inter / number_inter[0]
        # Converting the particle number to the particle mass
        mass_inter = np.zeros(number_inter.shape, dtype=float)
        for size_class in range(mass_inter.size):
            mass_inter[size_class] = number_inter[size_class] / (2 ** (settings.DN * size_class))
        # Normalize by the particle mass in size class k = 0
        mass_inter /= mass_inter[0]
        # Plotting the normalized count and mass inputs
        ax[2].plot(utils.size_range(units='mm', size_class_number=settings.SIZE_CLASS_NUMBER),
                   number_norm, marker=self.field_marker, linestyle=self.field_line, color='tab:cyan',
                   label=r'Zeri et al. (2021) Input')
        ax[2].legend(fontsize=self.legend_size, loc=self.legend_loc)
        norm_factor = utils.size_range(units='mm', size_class_number=settings.SIZE_CLASS_NUMBER)
        ax[3].plot(utils.size_range(units='mm', size_class_number=settings.SIZE_CLASS_NUMBER),
                   mass_inter / norm_factor, marker=self.field_marker, linestyle=self.field_line, color='tab:cyan',
                   label=r'Zeri et al. (2021) Input')
        ax[3].legend(fontsize=self.legend_size, loc=self.legend_loc)



def subfigure_title(index, beach_state):
    alphabet = string.ascii_lowercase
    title_dict = {'adrift_open_surf': 'Adrift - open ocean', 'adrift_10km_surf': 'Adrift - coastal', 'beach': 'Beach'}
    return '({}) {}'.format(alphabet[index], title_dict[beach_state])


def label(lambda_frag):
    return r'$\lambda_f=$' + '{} days'.format(lambda_frag)
