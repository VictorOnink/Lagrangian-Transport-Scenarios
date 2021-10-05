import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import string
import numpy as np


class FragmentationKaandorpPartial_FieldDataComp:
    def __init__(self, figure_direc, scenario, shore_time, lambda_frag_list, rho, sink):
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
        if sink:
            self.count, self.mass = 'particle_number', 'particle_mass'
        else:
            self.count, self.mass = 'particle_number_sink', 'particle_mass_sink'
        # Figure parameters
        self.fig_size = (14, 14)
        self.fig_shape = (self.beach_state_list.__len__(), 2)
        self.x_label = 'Size (mm)'
        self.y_label = r'Normalized Particle Number (n mm$^{-1}$)'
        self.twiny_label = r'Normalized Particle Mass (g mm$^{-1}$)'
        self.ax_ticklabel_size = 12
        self.ax_label_size = 14
        self.legend_size = 12
        self.xmin, self.xmax = 1e-1, 2e2
        self.ymin, self.ymax = 1e-4, 2e5
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.number_of_plots = self.fig_shape[0] * self.fig_shape[1]
        self.field_marker = 'x'
        self.field_line = '--'

    def plot(self):
        # Getting the sizes of the size classes, and we convert from meters to mm
        size_classes = utils.size_range(size_class_number=self.class_num, units='mm')

        # Loading the data
        data_dict = {}
        for lambda_frag in self.lambda_frag_list:
            data = vUtils.FragmentationKaandorpPartial_load_data(scenario=self.scenario, prefix=self.prefix,
                                                                 data_direc=self.data_direc, shore_time=self.shore_time,
                                                                 lambda_frag=lambda_frag, rho=self.rho,
                                                                 postprocess=True)
            data_dict[lambda_frag] = {}
            for beach_state in self.beach_state_list:
                data_dict[lambda_frag][beach_state] = {}
                data_dict[lambda_frag][beach_state][self.count] = data[beach_state][self.count]
                data_dict[lambda_frag][beach_state][self.mass] = data[beach_state][self.mass]
                print(data[beach_state][self.count][data['final_index']])
                print(data[beach_state][self.mass][data['final_index']])
        time_index = data['final_index']
        field_dict = utils.load_obj(vUtils.FragmentationKaandorpPartial_fielddata_filename())

        # Creating the figure
        ax, twin_ax = vUtils.base_figure(fig_size=self.fig_size, ax_range=self.ax_range, x_label=self.x_label,
                                         y_label=self.y_label, ax_ticklabel_size=self.ax_ticklabel_size,
                                         ax_label_size=self.ax_label_size, shape=self.fig_shape,
                                         plot_num=self.number_of_plots, log_yscale=True, log_xscale=True,
                                         all_x_labels=True, all_y_labels=False, add_twinx=True,
                                         twinx_y_label=self.twiny_label, twinx_ax_range=self.ax_range,
                                         log_twinxscale=True)

        # Labelling the subfigures
        for index_ax in range(self.number_of_plots):
            ax[index_ax].set_title(subfigure_title(index_ax, self.beach_state_list[index_ax // 2]),
                                   fontsize=self.ax_label_size)

        # Plotting the model distributions
        for ax_index, sub_ax in enumerate(ax):
            for lambda_index, lambda_frag in enumerate(self.lambda_frag_list):
                c = vUtils.discrete_color_from_cmap(index=lambda_index, subdivisions=self.lambda_frag_list.__len__())
                if ax_index % 2 == 0:
                    norm_factor = data_dict[lambda_frag][self.beach_state_list[ax_index // 2]][self.count][time_index][0]
                    sub_ax.plot(size_classes, data_dict[lambda_frag][self.beach_state_list[ax_index // 2]][self.count][time_index] / norm_factor,
                                linestyle='-', color=c, label=label(lambda_frag))
                else:
                    norm_factor = data_dict[lambda_frag][self.beach_state_list[ax_index // 2]][self.mass][time_index][0]
                    twin_ax[ax_index].plot(size_classes,
                                           data_dict[lambda_frag][self.beach_state_list[ax_index // 2]][self.mass][time_index] / norm_factor,
                                           linestyle='-', color=c, label=label(lambda_frag))
                print(norm_factor)

        # Field data - open ocean
        norm_factor = field_dict['Cozar']['pdf_counts'][14]
        ax[0].plot(field_dict['Cozar']['bin_midpoint'], field_dict['Cozar']['pdf_counts'] / norm_factor,
                   marker=self.field_marker, linestyle=self.field_line, color='tab:red', label='Cozar et al. (2015)')

        # Field data - coastal waters
        norm_factor = field_dict['RuizOrejon']['pdf_counts'][6]
        ax[2].plot(field_dict['RuizOrejon']['bin_midpoint'], field_dict['RuizOrejon']['pdf_counts'] / norm_factor,
                   marker=self.field_marker, linestyle=self.field_line, color='tab:red',
                   label=r'Ruiz-Orej$\`o$n et al. (2018)')

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
        # Field data - beach, microplastic mass
        norm_factor = field_dict['Fok']['pdf_mass'][5]
        twin_ax[5].plot(field_dict['Fok']['bin_midpoint'], field_dict['Fok']['pdf_mass'] / norm_factor,
                        marker=self.field_marker, linestyle=self.field_line, color='tab:red', label='Fok et al. (2017)')

        # Adding legends
        for sub_ax in ax:
            sub_ax.legend(fontsize=self.legend_size, loc='upper right')

        # Saving the figure
        str_format = self.shore_time, self.rho
        file_name = self.output_direc + 'SizeSpectrumFieldData-ST={}-rho={}.png'.format(*str_format)
        plt.savefig(file_name, bbox_inches='tight')


def subfigure_title(index, beach_state):
    alphabet = string.ascii_lowercase
    title_dict = {'adrift_open_surf': 'Adrift - open ocean', 'adrift_10km_surf': 'Adrift - coastal', 'beach': 'Beach'}
    return '({}) {}'.format(alphabet[index], title_dict[beach_state])


def label(lambda_frag):
    return r'$\lambda_f=$' + '{} days'.format(lambda_frag)
