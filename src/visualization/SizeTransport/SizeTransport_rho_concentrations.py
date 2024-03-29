import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string
import cmocean.cm as cmo


class SizeTransport_rho_concentrations:
    def __init__(self, scenario, figure_direc, size, beach_state, time_selection, rho_list, depth_level, tau=0):
        # Simulation parameters
        self.scenario = scenario
        self.rho_list = rho_list
        self.time_selection = time_selection
        self.beach_state = beach_state
        self.size = size
        self.tau = tau
        self.depth_level = depth_level
        # Data parameters
        self.output_direc = figure_direc + 'concentrations/'
        self.data_direc = settings.DATA_OUTPUT_DIREC + 'concentrations/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'horizontal_concentration'
        # Figure parameters
        self.figure_size = (20, 10)
        self.figure_shape = (2, 2)
        self.ax_label_size = 14
        self.ax_ticklabel_size = 14
        self.number_of_plots = self.rho_list.__len__()
        self.adv_file_dict = advection_files.AdvectionFiles(advection_scenario='CMEMS_MEDITERRANEAN').file_names
        self.spatial_domain = np.nanmin(self.adv_file_dict['LON']),  np.nanmax(self.adv_file_dict['LON']), \
                              np.nanmin(self.adv_file_dict['LAT']), np.nanmax(self.adv_file_dict['LAT'])
        self.cmap = cmo.thermal

    def plot(self):
        # Loading the data
        concentration_dict = {}
        if self.time_selection == 'average':
            key_concentration = "overall_concentration"
        else:
            key_concentration = utils.analysis_simulation_year_key(self.time_selection)
        for index, rho in enumerate(self.rho_list):
            data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix,
                                                       data_direc=self.data_direc,
                                                       size=self.size, rho=rho, tau=self.tau)
            concentration_dict[index] = data_dict[key_concentration][self.beach_state][self.depth_level]
        Lon, Lat = np.meshgrid(data_dict['lon'], data_dict['lat'])

        # Normalizing the concentration by the lowest non-zero concentration over all the sizes
        normalization_factor = 1e10
        for size in concentration_dict.keys():
            concentration = concentration_dict[size]
            min_non_zero = np.nanmin(concentration[concentration > 0])
            if min_non_zero < normalization_factor:
                normalization_factor = min_non_zero
        for size in concentration_dict.keys():
            concentration_dict[size] /= normalization_factor

        # Setting zero values to nan
        for size in concentration_dict.keys():
            concentration_dict[size][concentration_dict[size] == 0] = np.nan

        # Creating the base figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1] + 1, width_ratios=[1, 1, 0.1])

        ax_list = []
        for rows in range(self.figure_shape[0]):
            for columns in range(self.figure_shape[1]):
                ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                           domain=self.spatial_domain, label_size=self.ax_label_size,
                                                           lat_grid_step=5, lon_grid_step=10, resolution='10m'))

        # Setting the colormap, and adding a colorbar
        norm = set_normalization(self.beach_state)
        cbar_label, extend = r"Relative Concentration ($C/C_{min}$)", 'max'
        cmap = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        cax = fig.add_subplot(gs[:, -1])
        cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend=extend)
        cbar.set_label(cbar_label, fontsize=self.ax_label_size)
        cbar.ax.tick_params(which='major', labelsize=self.ax_ticklabel_size, length=14, width=2)
        cbar.ax.tick_params(which='minor', labelsize=self.ax_ticklabel_size, length=7, width=2)

        # Adding subfigure titles
        for index, ax in enumerate(ax_list):
            ax.set_title(subfigure_title(index, self.size, self.rho_list), weight='bold', fontsize=self.ax_label_size)

        # The actual plotting of the figures
        for index in range(self.rho_list.__len__()):
            if self.beach_state in ['adrift']:
                ax_list[index].pcolormesh(Lon, Lat, concentration_dict[index], norm=norm, cmap=self.cmap,
                                          zorder=200)
            else:
                ax_list[index].scatter(Lon.flatten(), Lat.flatten(), c=concentration_dict[index], norm=norm,
                                       cmap=self.cmap,
                                       zorder=200)

        # Saving the figure
        file_name = self.plot_save_name()
        plt.savefig(file_name, bbox_inches='tight')
        plt.close('all')

    def plot_save_name(self, flowdata='CMEMS_MEDITERRANEAN', startyear=2010):
        selection_dict = {'average': 'TotalAverage', 0: 'year_0', 1: 'year_1', 2: 'year_2'}
        str_format = flowdata, self.beach_state, selection_dict[self.time_selection], self.size * 1E6, startyear, self.depth_level
        return self.output_direc + 'rho/SizeTransport_RHO_{}_{}_{}_size_{:.3f}_y_{}_{}.png'.format(*str_format)


def set_normalization(beach_state):
    """
    Setting the normalization that we use for the colormap
    :param beach_state: adrift, beach or seabed
    :return:
    """
    if beach_state == 'adrift':
        vmin, vmax = 1, 1e4
    elif beach_state == 'beach':
        vmin, vmax = 1, 1e4
    return colors.LogNorm(vmin=vmin, vmax=vmax)


def subfigure_title(index, size, rho_list):
    return '({}) r = {:.3f} mm, '.format(string.ascii_lowercase[index], size * 1e3) + r'$\rho=$' + \
           '{} kg m'.format(rho_list[index]) + r'$^{-3}$'




