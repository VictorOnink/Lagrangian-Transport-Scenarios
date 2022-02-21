import settings
import utils
from netCDF4 import Dataset
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean.cm as cmo
import cartopy.crs as ccrs
import string
from Analysis.CMEMS_mediterranean_mean_MLD import CMEMS_mediterranean_mean_MLD


class General_season_average:
    def __init__(self, scenario, figure_direc, variable):
        # Scenario specific variables
        self.scenario = scenario
        self.figure_direc = figure_direc
        self.variable = variable
        # Data variables
        self.output_direc = figure_direc + 'General/'
        utils.check_direc_exist(self.output_direc)
        self.file_dict = self.scenario.file_dict
        # Figure variables
        self.figure_size = (16, 12)
        self.figure_shape = (2, 2)
        self.ax_label_size = 16
        self.ax_ticklabel_size = 14
        self.spatial_domain = np.nanmin(self.file_dict['LON']), np.nanmax(self.file_dict['LON']), \
                              np.nanmin(self.file_dict['LAT']), np.nanmax(self.file_dict['LAT'])
        self.season_list = {'wind': ['DJF', 'MAM', 'JJA', 'SON'], 'MLD': ['JFM', 'AMJ', 'JAS', 'OND']}[self.variable]
        self.cmap = {'wind': cmo.speed, 'MLD': cmo.haline_r}[self.variable]
        self.step = 3

    def plot(self):
        # Loading the data
        variable_dict = {}
        for season in self.season_list:
            if self.variable == 'MLD':
                MLD_data = CMEMS_mediterranean_mean_MLD().calculate_seasonal_mean()
                variable_dict[season] = (MLD_data[2010][season] + MLD_data[2011][season] + MLD_data[2012][season]) / 3
            elif self.variable == 'wind':
                dataset = Dataset(file_name_variable(variable=self.variable, season=season))
                variable_dict[season] = {}
                u10, v10 = np.nanmean(dataset.variables['u10'][:], axis=0), np.nanmean(dataset.variables['v10'][:], axis=0)
                variable_dict[season]['magnitude'] = np.sqrt(np.square(u10) + np.square(v10))
                variable_dict[season]['u10'] = np.divide(u10, variable_dict[season]['magnitude'])
                variable_dict[season]['v10'] = np.divide(v10, variable_dict[season]['magnitude'])

        # Creating the figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1] + 1, width_ratios=[1, 1, 0.1])

        ax_list = []
        for rows in range(self.figure_shape[0]):
            for columns in range(self.figure_shape[1]):
                ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                           domain=self.spatial_domain, add_gridlines=False,
                                                           resolution='10m', land_zorder=90, ocean_zorder=0))
        # Adding subplot titles
        for ax_index, ax in enumerate(ax_list):
            ax.set_title(subfigure_title(ax_index, self.season_list[ax_index]), fontsize=self.ax_label_size)

        # Setting the colormap normalization
        if self.variable == 'MLD':
            norm, extend = colors.LogNorm(vmin=1e0, vmax=1e2), 'max'
            cbar_label = 'Depth (m)'
        elif self.variable == 'wind':
            norm, extend = colors.LogNorm(vmin=1e-2, vmax=1e0), 'both'
            cbar_label = r"Wind speed (m s$^{-1}$)"
        cmap = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        cax = fig.add_subplot(gs[:, -1])
        cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend=extend)
        cbar.set_label(cbar_label, fontsize=self.ax_label_size)
        cbar.ax.tick_params(which='major', labelsize=self.ax_ticklabel_size, length=14, width=2)
        cbar.ax.tick_params(which='minor', labelsize=self.ax_ticklabel_size, length=7, width=2)

        # Plotting the actual mean variables
        if self.variable == 'MLD':
            Lon, Lat = np.meshgrid(dataset.variables['lon'][:], dataset.variables['lat'][:])
        elif self.variable == 'wind':
            Lon, Lat = np.meshgrid(dataset.variables['longitude'][:], dataset.variables['latitude'][:])
        for ax_index, ax in enumerate(ax_list):
            if self.variable == 'MLD':
                ax.pcolormesh(Lon, Lat, variable_dict[self.season_list[ax_index]], cmap=self.cmap, zorder=10,
                              norm=norm)
            if self.variable == 'wind':
                ax.pcolormesh(Lon, Lat, variable_dict[self.season_list[ax_index]]['magnitude'], cmap=self.cmap,
                              zorder=10, norm=norm)
                ax.quiver(Lon[::self.step, ::self.step], Lat[::self.step, ::self.step],
                          variable_dict[self.season_list[ax_index]]['u10'][::self.step, ::self.step],
                          variable_dict[self.season_list[ax_index]]['v10'][::self.step, ::self.step], zorder=15,
                          scale=7, scale_units='inches')

        # Saving the figure
        str_format = self.variable, settings.ADVECTION_DATA
        file_name = self.output_direc + 'Seasonal_average_{}_{}_2010-12.png'.format(*str_format)
        plt.savefig(file_name, bbox_inches='tight')


def subfigure_title(index, season):
    alphabet = string.ascii_lowercase
    return '({}) {}'.format(alphabet[index], season)


def file_name_variable(variable, season):
    if variable == 'MLD':
        return settings.DATA_DIR_SERVERS[settings.SERVER] + 'CMEMS_MED/mean_{}_AMXL.nc'.format(season)
    elif variable == 'wind':
        return settings.DATA_DIR_SERVERS[settings.SERVER] + 'Wind/{}_2010_2012_wind_mean.nc'.format(season)