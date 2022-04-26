import settings
import utils
from netCDF4 import Dataset
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean.cm as cmo
import cartopy.crs as ccrs


class General_average_wind_speed:
    def __init__(self, scenario, figure_direc):
        # Scenario specific variables
        self.scenario = scenario
        self.figure_direc = figure_direc
        # Data variables
        self.output_direc = figure_direc + 'General/'
        utils.check_direc_exist(self.output_direc)
        self.data_direc = settings.DATA_DIREC + 'Wind/'
        self.lon_grid, self.lat_grid = self.scenario.file_dict['LON'], self.scenario.file_dict['LAT']
        # Figure variables
        self.figure_size = (10, 8)
        self.figure_shape = (1, 1)
        self.ax_label_size = 14
        self.ax_ticklabel_size = 12
        self.spatial_domain = np.nanmin(self.lon_grid), np.nanmax(self.lon_grid), \
                              np.nanmin(self.lat_grid), np.nanmax(self.lat_grid)
        self.cmap = cmo.speed

    def plot(self):
        # Loading the data
        wind_file = self.data_direc + 'ERA5-wind10m-average.nc'
        dataset = Dataset(wind_file)
        u10, v10 = dataset.variables['u10'][0, :, :], dataset.variables['v10'][0, :, :]
        Lon, Lat = np.meshgrid(dataset.variables['longitude'][:], dataset.variables['latitude'][:])


        # Getting the wind speed
        wind_magnitude = np.sqrt(np.square(u10) + np.square(v10))

        # Creating the figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1] + 1, width_ratios=[1, 0.1])

        ax = []
        for rows in range(self.figure_shape[0]):
            for columns in range(self.figure_shape[1]):
                ax.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                      domain=self.spatial_domain, lat_grid_step=5, lon_grid_step=10,
                                                      resolution='10m', land_zorder=90))
        # Setting the normalization of the colormap
        normalization = colors.LogNorm(vmin=1e-1, vmax=1e1)

        # The actual plotting
        depth_plot = plt.pcolormesh(Lon, Lat, wind_magnitude, norm=normalization, cmap=self.cmap, zorder=50,
                                    transform=ccrs.PlateCarree())

        cax = fig.add_subplot(gs[:, -1])
        cbar = plt.colorbar(depth_plot, cax=cax, orientation='vertical', extend='both')
        cbar.set_label(r'Mean Wind Speed (m s$^{-1}$)', fontsize=self.ax_label_size)
        cbar.ax.tick_params(which='major', labelsize=self.ax_ticklabel_size, length=14, width=2)
        cbar.ax.tick_params(which='minor', labelsize=self.ax_ticklabel_size, length=7, width=2)

        file_name = self.output_direc + 'AverageWindSpeed_2010-2015.png'
        plt.savefig(file_name, bbox_inches='tight')






