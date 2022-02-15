import settings
import utils
from netCDF4 import Dataset
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean.cm as cmo


class General_bathymetry:
    def __init__(self, scenario, figure_direc):
        # Scenario specific variables
        self.scenario = scenario
        self.figure_direc = figure_direc
        # Data variables
        self.output_direc = figure_direc + 'General/'
        utils.check_direc_exist(self.output_direc)
        # Figure variables
        self.figure_size = (10, 8)
        self.figure_shape = (1, 1)
        self.ax_label_size = 14
        self.ax_ticklabel_size = 12
        self.cmap = cmo.speed

    def plot(self):
        # Loading the data
        dataset = Dataset(self.scenario.file_dict['BATH_filenames'])
        bath_dict = {}
        for variable in ['LON', 'LAT', 'DEPTH']:  # Note, DEPTH is actually the bathymetry
            bath_dict[variable] = dataset.variables[self.scenario.file_dict['BATH_variables'][variable]][:]

        # Setting the zero depths to nan values
        bath_dict['DEPTH'][bath_dict['DEPTH'] == 0] = np.nan

        # Getting just the regions where the depth is between 10 and 11 meters
        bath_dict['DEPTH'][bath_dict['DEPTH'] < 1] = np.nan
        bath_dict['DEPTH'][bath_dict['DEPTH'] > 100] = np.nan

        depth_list = bath_dict['DEPTH'].flatten()
        selection = ~np.isnan(depth_list)
        Lon, Lat = np.meshgrid(bath_dict['LON'], bath_dict['LAT'])
        lon_list = Lon.flatten()[selection]
        lat_list = Lat.flatten()[selection]

        # Getting the spatial domain for the figure
        spatial_domain = np.nanmin(bath_dict['LON']), np.nanmax(bath_dict['LON']), \
                         np.nanmin(bath_dict['LAT']), np.nanmax(bath_dict['LAT'])

        # Creating the figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1] + 1, width_ratios=[1, 0.1])

        ax = []
        for rows in range(self.figure_shape[0]):
            for columns in range(self.figure_shape[1]):
                ax.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                      domain=spatial_domain, lat_grid_step=5, lon_grid_step=10,
                                                      resolution='10m'))
        # Setting the normalization of the colormap
        normalization = colors.LogNorm(vmin=1e1, vmax=5e3)

        # The actual plotting
        # Lon, Lat = np.meshgrid(bath_dict['LON'], bath_dict['LAT'])
        # depth_plot = plt.pcolormesh(Lon, Lat, bath_dict['DEPTH'], norm=normalization, cmap=self.cmap, zorder=1000)

        ax[0].plot(lon_list, lat_list, '.r', zorder=10000)

        # cax = fig.add_subplot(gs[:, -1])
        # cbar = plt.colorbar(depth_plot, cax=cax, orientation='vertical', extend='both')
        # cbar.set_label('Depth (m)', fontsize=self.ax_label_size)
        # cbar.ax.tick_params(which='major', labelsize=self.ax_ticklabel_size, length=14, width=2)
        # cbar.ax.tick_params(which='minor', labelsize=self.ax_ticklabel_size, length=7, width=2)

        file_name = self.output_direc + 'Bathymetry.png'
        plt.savefig(file_name, bbox_inches='tight')
