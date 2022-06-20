import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import string
from datetime import datetime, timedelta
import numpy as np


class SizeTransport_reservoirs_Defense:
    def __init__(self, scenario, figure_direc, rho_list=[920], single_plot=False,
                 fixed_resus=False, resus_time=7):
        utils.print_statement('Creating the SizeTransport reservoirs figure', to_print=True)
        # Simulation parameters
        self.scenario = scenario
        self.rho_list = rho_list
        self.size_list = np.array([5000, 2500, 1250, 625, 313, 156, 78, 39, 20, 10, 5, 2]) * settings.SIZE_FACTOR
        self.tau = 0.0
        self.fixed_resus = fixed_resus
        self.resus_time = resus_time
        # Data parameters
        self.output_direc = figure_direc + 'timeseries/'
        self.data_direc = settings.DATA_OUTPUT_DIREC + 'timeseries/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'timeseries'
        # Figure parameters
        self.figure_size = (14, 8)
        self.ax_label_size = 18
        self.ax_ticklabel_size = 16
        self.legend_size = 16
        self.markersize = 140
        self.markersize_legend = 10
        self.xmax, self.xmin = 1e1, 1e-3
        self.ymax, self.ymin = 100, 0
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.beach_state_list = ['beach', 'adrift']
        self.coastal_list = ['total', 'coastal_zone', 'coastal_10km', 'coastal_20km']
        self.y_label = 'Percentage of total particles per size class (%)'
        self.x_label = 'Size (mm)'
        self.rho_marker_dict = {30: 'X', 920: 'o', 980: 's', 1020: 'D'}
        self.state_color = {'beach': 'r', 'adrift': 'b'}
        self.single_plot = single_plot
        self.figure_size = {True: (12, 10), False: (16, 8)}[self.single_plot]
        self.number_of_plots = {True: 1, False: 2}[self.single_plot]
        self.width_ratio = {1: [1, 0.3], 2: [1, 1, 0.4]}[self.number_of_plots]
        self.figure_shape = (1, self.number_of_plots)

    def plot(self):
        # Loading the data
        timeseries_dict = dict.fromkeys(self.rho_list)
        for rho in self.rho_list:
            timeseries_dict[rho] = {}
            for size in self.size_list:
                timeseries_dict[rho][size] = {}
                for fixed_resus in [True, False]:
                    timeseries_dict[rho][size][fixed_resus] = {}
                    if fixed_resus:
                        data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix,
                                                                   data_direc=self.data_direc, size=size, rho=rho,
                                                                   tau=self.tau, fixed_resus=fixed_resus,
                                                                   resus_time=self.resus_time)
                    else:
                        data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix,
                                                                   data_direc=self.data_direc, size=size, rho=rho,
                                                                   tau=self.tau, fixed_resus=fixed_resus,
                                                                   resus_time=utils.get_resuspension_timescale(L=size, rho_p=rho))

                    for beach_state in self.beach_state_list:
                        if beach_state in ['beach']:
                            timeseries_dict[rho][size][fixed_resus][beach_state] = data_dict[beach_state]
                        elif beach_state in ['adrift']:
                            timeseries_dict[rho][size][fixed_resus][beach_state] = {}
                            for coastal in self.coastal_list:
                                timeseries_dict[rho][size][fixed_resus][beach_state][coastal] = data_dict[beach_state][coastal]
                    timeseries_dict[rho][size][fixed_resus]['total_divide'] = data_dict['total'][-1]

        # Normalizing all the particle counts with the total number of particles, and then multiplying by 100 to get a
        # percentage
        for rho in self.rho_list:
            for size in self.size_list:
                for fixed_resus in [False, True]:
                    for beach_state in self.beach_state_list:
                        if beach_state in ['beach']:
                            timeseries_dict[rho][size][fixed_resus][beach_state] /= timeseries_dict[rho][size][fixed_resus]['total_divide']
                            timeseries_dict[rho][size][fixed_resus][beach_state] *= 100.
                            timeseries_dict[rho][size][fixed_resus][beach_state] = np.nanmean(timeseries_dict[rho][size][fixed_resus][beach_state])
                        elif beach_state in ['adrift']:
                            for c in self.coastal_list:
                                timeseries_dict[rho][size][fixed_resus][beach_state][c] /= timeseries_dict[rho][size][fixed_resus]['total_divide']
                                timeseries_dict[rho][size][fixed_resus][beach_state][c] *= 100.
                                timeseries_dict[rho][size][fixed_resus][beach_state][c] = np.nanmean(timeseries_dict[rho][size][fixed_resus][beach_state][c])

        # Creating the figure
        ax = vUtils.base_figure(fig_size=self.figure_size, ax_range=self.ax_range, y_label=self.y_label,
                                x_label=self.x_label, ax_label_size=self.ax_label_size,
                                ax_ticklabel_size=self.ax_ticklabel_size, shape=self.figure_shape,
                                plot_num=self.number_of_plots, legend_axis=True, log_yscale=False, log_xscale=True,
                                width_ratios=self.width_ratio, all_x_labels=True)

        # If we have multiple subplots, labelling the subplots
        if self.single_plot:
            title = {True: r'$\lambda_R =$' + '{}'.format(self.resus_time) + ' days',
                     False: 'Size-dependent resuspension'}[self.fixed_resus]
            ax[0].set_title(title, weight='bold', fontsize=self.ax_label_size)
        else:
            title_list = ['(a) Size-dependent resuspension',
                          r' (b) $\lambda_R =$' + '{}'.format(self.resus_time) + ' days']
            for index in range(self.number_of_plots):
                ax[index].set_title(title_list[index], fontsize=self.ax_label_size)

        # Now, adding in the actual data
        if self.single_plot:
            for rho in self.rho_list:
                for index_size, size in enumerate(self.size_list):
                    for beach_state in self.beach_state_list:
                        if beach_state in ['beach']:
                            ax[0].scatter(size * 1e3, timeseries_dict[rho][size][self.fixed_resus][beach_state],
                                          marker=self.rho_marker_dict[rho], edgecolors='b',
                                          facecolors='b', s=self.markersize)
                        elif beach_state in ['adrift']:
                            coastal_zone = 'coastal_10km'
                            ax[0].scatter(size * 1e3, timeseries_dict[rho][size][self.fixed_resus][beach_state][coastal_zone],
                                          marker=self.rho_marker_dict[rho], edgecolors='g',
                                          facecolors='g', s=self.markersize)
                            # To get the open ocean, we subtract the coastal fraction from the total fraction
                            open_ocean = timeseries_dict[rho][size][self.fixed_resus][beach_state]['total'] - timeseries_dict[rho][size][self.fixed_resus][beach_state][coastal_zone]
                            ax[0].scatter(size * 1e3, open_ocean,
                                          marker=self.rho_marker_dict[rho], edgecolors='r',
                                          facecolors='r', s=self.markersize)
        else:
            for rho in self.rho_list:
                for index_size, size in enumerate(self.size_list):
                    for beach_state in self.beach_state_list:
                        for index_fix, fixed_resus in enumerate([False, True]):
                            if beach_state in ['beach']:
                                ax[index_fix].scatter(size * 1e3, timeseries_dict[rho][size][fixed_resus][beach_state],
                                                      marker=self.rho_marker_dict[rho], edgecolors='b',
                                                      facecolors='none', s=self.markersize)
                            elif beach_state in ['adrift']:
                                coastal_zone = 'coastal_10km'
                                ax[index_fix].scatter(size * 1e3, timeseries_dict[rho][size][fixed_resus][beach_state][coastal_zone],
                                                      marker=self.rho_marker_dict[rho], edgecolors='g',
                                                      facecolors='none', s=self.markersize)
                                # To get the open ocean, we subtract the coastal fraction from the total fraction
                                open_ocean = timeseries_dict[rho][size][fixed_resus][beach_state]['total'] - timeseries_dict[rho][size][fixed_resus][beach_state][coastal_zone]
                                ax[index_fix].scatter(size * 1e3, open_ocean,
                                                      marker=self.rho_marker_dict[rho], edgecolors='r',
                                                      facecolors='none', s=self.markersize)

        # Creating a legend
        # rho_lines = [plt.plot([], [], self.rho_marker_dict[rho], c='k', label=r'$\rho=$' + str(rho) + r' kg m$^{-3}$',
        #                       )[0] for rho in self.rho_list]
        beach = [plt.plot([], [], 'o', c='b', label='Beach', markersize=self.markersize_legend)[0]]
        coastal = [plt.plot([], [], 'o', c='g', label='Coastal', markersize=self.markersize_legend)[0]]
        open = [plt.plot([], [], 'o', c='r', label='Open water', markersize=self.markersize_legend)[0]]
        ax[-1].legend(handles=beach + coastal + open, fontsize=self.legend_size, loc='upper right')
        ax[-1].axis('off')

        file_name = self.file_name()
        plt.savefig(file_name, bbox_inches='tight', dpi=400)
        plt.close('all')

    def file_name(self, file_type='.png'):
        name = 'SizeTransport_reservoirs_DEFENSE_'
        if self.rho_list.__len__() != 4:
            name += 'rho_{}_'.format(self.rho_list[0])
        if self.single_plot:
            if self.fixed_resus:
                name += 'fixed_resus_{}'.format(self.resus_time)
            else:
                name += 'size_dependent'
        else:
            name += '{}'.format(self.resus_time)
        return self.output_direc + name + file_type


def beach_label(beach_state):
    if beach_state == 'beach':
        return 'Beach'
    elif beach_state == 'adrift':
        return 'Adrift'


