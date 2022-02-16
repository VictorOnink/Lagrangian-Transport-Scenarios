import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string
import cmocean.cm as cmo


class SizeTransport_concentration_seasons:
    def __init__(self, scenario, figure_direc, size, time_selection, depth_level, beach_state='adrift', rho=920, tau=0):
        # Simulation parameters
        self.scenario = scenario
        self.rho = rho
        self.time_selection = time_selection
        self.beach_state = beach_state
        self.size = size
        self.tau = tau
        self.depth_level = depth_level
        # Data parameters
        self.output_direc = figure_direc + 'concentrations/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'horizontal_concentration_monthly'
        # Figure parameters
        self.figure_size = (20, 10)
        self.figure_shape = (2, 2)
        self.ax_label_size = 14
        self.ax_ticklabel_size = 14
        self.number_of_plots = 4
        self.adv_file_dict = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                            advection_scenario='CMEMS_MEDITERRANEAN',
                                                            repeat_dt=None).file_names
        self.spatial_domain = np.nanmin(self.adv_file_dict['LON']), np.nanmax(self.adv_file_dict['LON']), \
                              np.nanmin(self.adv_file_dict['LAT']), np.nanmax(self.adv_file_dict['LAT'])
        self.cmap = cmo.thermal

    def plot(self):
        # Loading the data
        concentration_dict = {}
        key_concentration = utils.analysis_simulation_year_key(self.time_selection)
        for month in range(1, 13):
            data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix,
                                                       data_direc=self.data_direc,
                                                       size=self.size, rho=self.rho, tau=self.tau)
            concentration_dict[month] = data_dict[key_concentration][month][self.beach_state][self.depth_level]
        Lon, Lat = np.meshgrid(data_dict['lon'], data_dict['lat'])

        season_dict = {}
        for season, month_list in enumerate([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]):
            season_dict[season] = np.zeros(shape=(data_dict['lat'].size, data_dict['lon'].size))
            print('season_dict[season] {}'.format(season_dict[season].shape))
            print('concentration_dict[month] {}'.format(concentration_dict[month].shape))
            for month in month_list:
                season_dict[season] += concentration_dict[month]
            season_dict[season] /= month_list.__len__()

        # Normalizing the concentration by the lowest non-zero concentration over all the sizes
        normalization_factor = 1e10
        for season in season_dict.keys():
            min_non_zero = np.nanmin(season_dict[season][season_dict[season] > 0])
            if min_non_zero < normalization_factor:
                normalization_factor = min_non_zero
        for season in season_dict.keys():
            season_dict[season] /= normalization_factor

        # Setting zero values to nan
        for season in season_dict.keys():
            season_dict[season][season_dict[season] == 0] = np.nan

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
        norm = colors.LogNorm(vmin=1, vmax=1e4)
        cbar_label, extend = r"Relative Concentration ($C/C_{min}$)", 'max'
        cmap = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        cax = fig.add_subplot(gs[:, -1])
        cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend=extend)
        cbar.set_label(cbar_label, fontsize=self.ax_label_size)
        cbar.ax.tick_params(which='major', labelsize=self.ax_ticklabel_size, length=14, width=2)
        cbar.ax.tick_params(which='minor', labelsize=self.ax_ticklabel_size, length=7, width=2)

        # Adding subfigure titles
        season_title = ['Winter (JFM)', 'Spring (AMJ)', 'Summer (JAS)', 'Fall (OND)']
        for index, ax in enumerate(ax_list):
            ax.set_title(season_title[index], weight='bold', fontsize=self.ax_label_size)

        # The actual plotting of the figures
        for index in range(self.rho_list.__len__()):
            ax_list[index].pcolormesh(Lon, Lat, season_dict[index], norm=norm, cmap=self.cmap, zorder=200)

        # Saving the figure
        file_name = self.plot_save_name()
        plt.savefig(file_name, bbox_inches='tight')
        plt.close('all')

    def plot_save_name(self, flowdata='CMEMS_MEDITERRANEAN', startyear=2010):
        selection_dict = {0: 'year_0', 1: 'year_1', 2: 'year_2'}
        str_format = flowdata, self.beach_state, selection_dict[self.time_selection], self.size * 1E6, startyear, \
                     self.depth_level
        return self.output_direc + 'rho/SizeTransport_SEASON_{}_{}_{}_size_{:.3f}_y_{}_{}.png'.format(*str_format)
