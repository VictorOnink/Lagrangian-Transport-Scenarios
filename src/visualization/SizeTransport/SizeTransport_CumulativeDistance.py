import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import numpy as np
from advection_scenarios import advection_files
from copy import deepcopy


class SizeTransport_CumulativeDistance:
    def __init__(self, scenario, figure_direc, size_list, simulation_years=2, rho=920, tau=0.0):
        # Simulation parameters
        self.scenario = scenario
        self.rho = rho
        self.size_list = size_list
        self.tau = tau
        self.simulation_years = simulation_years
        # Data parameters
        self.output_direc = figure_direc + 'mix/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'statistics/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'basic_statistics'
        self.advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                                 advection_scenario='CMEMS_MEDITERRANEAN',
                                                                 repeat_dt=None)
        # Figure parameters
        self.figure_size = (20, 8)
        self.variable_list = ['z', 'distance_vertical', 'distance_horizontal']
        self.number_of_plots = self.variable_list.__len__()
        self.figure_shape = (1, self.variable_list.__len__() + 1)
        self.ax_label_size = 14
        self.ax_ticklabel_size = 14
        self.legend_size = 12
        self.max_depth = np.nanmax(self.advection_scenario.file_names['DEPTH'])
        self.variable_domain = {'z': np.arange(0, self.max_depth, 0.1),
                                'distance_vertical': np.logspace(0, 8, num=1000),
                                'distance_horizontal': np.logspace(0, 8, num=1000)}

    def plot(self):
        # Loading the data
        variable_dict = dict.fromkeys(self.variable_list)
        for key in variable_dict.keys():
            variable_dict[key] = np.zeros(shape=self.variable_domain[key].shape, dtype=float)
        size_dict = dict.fromkeys(self.size_list)
        for keys in size_dict.keys():
            size_dict[keys] = deepcopy(variable_dict)
        for index_size, size in enumerate(self.size_list):
            data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix, data_direc=self.data_direc,
                                                       size=size, rho=self.rho, tau=self.tau, restart=self.simulation_years)
            for variable in self.variable_list:
                size_dict[size][variable] += cumulative_fraction(data_dict, variable, self.variable_domain)

        # Create the figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1],
                              width_ratios=[1] * self.number_of_plots + [0.35])
        ax = []
        for column in range(gs.ncols):
            ax_sub = fig.add_subplot(gs[0, column])
            if column == 0:
                ax_sub.set_ylabel(r'Fraction of Total (%)', fontsize=self.ax_label_size)
            ax_sub.tick_params(axis='both', labelsize=self.ax_ticklabel_size)
            ax_sub.set_ylim([0, 100])
            ax.append(ax_sub)
        # Plotting the max depth
        ax[0].set_xlabel('Max depth (m)', fontsize=self.ax_label_size)
        ax[0].set_xscale('log')
        ax[0].set_xlim([1e0, 1e3])
        # the cumulative vertical distance
        ax[1].set_xlabel('Cumulative vertical distance (m)', fontsize=self.ax_label_size)
        ax[1].set_xscale('log')
        ax[1].set_xlim([1e1, 5e6])
        # the cumulative horizontal distance
        ax[2].set_xlabel('Cumulative horizontal distance (km)', fontsize=self.ax_label_size)
        ax[2].set_xscale('log')
        ax[2].set_xlim([1e1, 2e5])

        # Plotting the data
        for index_size, size in enumerate(self.size_list):
            for index_var, variable in enumerate(self.variable_list):
                line_color = vUtils.discrete_color_from_cmap(index_size, subdivisions=self.size_list.__len__())
                ax[index_var].plot(self.variable_domain[variable], size_dict[size][variable], linestyle='-',
                                   color=line_color, label=size_label(size))
        # Creating a legend
        size_colors = [plt.plot([], [], c=vUtils.discrete_color_from_cmap(index_size, subdivisions=self.size_list.__len__()),
                                label=size_label(size), linestyle='-')[0] for index_size, size in enumerate(self.size_list)]
        ax_legend = fig.add_subplot(gs[:, 3])
        ax_legend.legend(handles=size_colors, fontsize=self.legend_size, loc='upper right')
        ax_legend.axis('off')

        # Saving the figure
        file_name = self.output_direc + 'Vertical_Horizontal_Distance_cumulative.png'
        plt.savefig(file_name, bbox_inches='tight')


def size_label(size):
    return r'r = {:.3f} mm'.format(size * 1e3)


def cumulative_fraction(data_dict, variable, variable_domain):
    var_data = data_dict[variable]['total']['max']
    particle_number = var_data.size
    output_array = np.zeros(variable_domain[variable].size, dtype=float)
    for step in range(output_array.size):
        output_array[step] += np.nansum(var_data < variable_domain[variable][step]) / particle_number * 100.
    return output_array
