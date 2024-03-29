import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string
import cmocean.cm as cmo


class SizeTransport_MaxDistance:
    def __init__(self, scenario, figure_direc, rho, tau=0, fixed_resus=False, resus_time=50, subselection=False):
        # Simulation parameters
        self.scenario = scenario
        self.rho = rho
        self.subselection = subselection
        self.size_list = np.array({False: [5000, 2500, 1250, 625, 313, 156, 78, 39, 20, 10, 5, 2],
                                   True: [5000, 2]}[self.subselection]) * settings.SIZE_FACTOR
        self.tau = tau
        self.fixed_resus = fixed_resus
        self.resus_time = resus_time
        # Data parameters
        self.output_direc = figure_direc + 'max_distance/'
        self.data_direc = settings.DATA_OUTPUT_DIREC + 'max_distance/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'maximum_distance'
        # Figure parameters
        self.figure_size = {False: (20, 10), True: (20, 8)}[self.subselection]
        self.ax_label_size = 16
        self.ax_ticklabel_size = 12
        self.number_of_plots = self.size_list.__len__()
        self.figure_shape = {False: (4, 3), True: (1, 2)}[self.subselection]
        self.width_ratios = {False: [1, 1, 1, 0.1], True: [1, 1, 0.07]}[self.subselection]
        self.adv_file_dict = advection_files.AdvectionFiles(advection_scenario='CMEMS_MEDITERRANEAN').file_names
        self.spatial_domain = np.nanmin(self.adv_file_dict['LON']), np.nanmax(self.adv_file_dict['LON']), \
                              np.nanmin(self.adv_file_dict['LAT']), np.nanmax(self.adv_file_dict['LAT'])
        self.cmap = cmo.haline_r

    def plot(self):
        # Loading the data
        data_dict = {}
        for size in self.size_list:
            size_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix,
                                                       data_direc=self.data_direc, fixed_resus=self.fixed_resus,
                                                       size=size, rho=self.rho, tau=self.tau,
                                                       resus_time=self.resus_time)
            data_dict[size] = size_dict['median_max_distance']
        LON, LAT = size_dict['site_lon'], size_dict['site_lat']

        # Creating the base figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1] + 1,
                              width_ratios=self.width_ratios)

        ax_list = []
        for rows in range(self.figure_shape[0]):
            for columns in range(self.figure_shape[1]):
                ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                           domain=self.spatial_domain, label_size=self.ax_ticklabel_size,
                                                           lat_grid_step=5, lon_grid_step=10, resolution='10m',
                                                           auto_aspect=False))

        # Setting the colormap, and adding a colorbar
        norm = colors.LogNorm(vmin=1, vmax=400)
        cbar_label, extend = r"Median maximum distance from shore (km)", None
        cmap = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        cax = fig.add_subplot(gs[:, -1])
        cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend=extend)
        cbar.set_label(cbar_label, fontsize=self.ax_label_size)
        cbar.ax.tick_params(which='major', labelsize=self.ax_ticklabel_size, length=14, width=2)
        cbar.ax.tick_params(which='minor', labelsize=self.ax_ticklabel_size, length=7, width=2)

        # Adding subfigure titles
        for index, ax in enumerate(ax_list):
            ax.set_title(self.subfigure_title(index), fontsize=self.ax_label_size)

        # Plotting the maximum distances
        for index, size in enumerate(self.size_list):
            ax_list[index].scatter(LON, LAT, cmap=self.cmap, c=data_dict[size], norm=norm, zorder=50, s=8)

        # Saving the figure
        plt.savefig(self.plot_save_name(), bbox_inches='tight', dpi=400)
        plt.close('all')

    def plot_save_name(self, file_type='.png'):
        name = self.output_direc + 'Max_distance_rho={}'.format(self.rho)
        if self.fixed_resus:
            name += '_fixed_resus_{}'.format(self.resus_time)
        if self.subselection:
            name += '_subselect'
        return name + file_type

    def subfigure_title(self, index):
        return '({}) d = {:.3f} mm'.format(string.ascii_lowercase[index], self.size_list[index] * 1e3)
