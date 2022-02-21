import settings
import utils
import xarray as xr
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean.cm as cmo
from matplotlib.ticker import AutoMinorLocator


class General_mean_tidal_Kz:
    def __init__(self, scenario, figure_direc):
        # Scenario specific variables
        self.scenario = scenario
        self.lon_grid, self.lat_grid = self.scenario.file_dict['LON'], self.scenario.file_dict['LAT']
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
        # Load the data
        Kz_data = self.load_data_file()

        # Calculate the vertical mean Kz value
        Kz_mean = np.nanmean(Kz_data['TIDAL_Kz'], axis=(1, 2))

        # Creating the figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1])

        ax = fig.add_subplot(gs[0, 0])
        ax.set_yscale('log')
        # ax.set_ylim([3000, 1])
        ax.set_xscale('log')
        # ax.set_xlim([1e-7, 1e-3])

        ax.set_ylabel('Depth (m)')
        ax.set_xlabel(r'Tidal $K_z$ (m$^2$ s$^{-1}$)')

        ax.plot(Kz_data['depth'], Kz_mean)

        ax.set_aspect('auto', adjustable=None)
        file_name = self.output_direc + 'Tidal_Kz_mean.png'
        plt.savefig(file_name, bbox_inches='tight')

    def load_data_file(self):
        # Loading the Kz file
        Kz_data = xr.load_dataset(utils.get_input_directory(server=settings.SERVER) + 'CMEMS_MEDITERRANEAN_Kz_TIDAL.nc')
        data_dict = {'depth': Kz_data.depth.values, 'LON': Kz_data.lon.values, 'LAT': Kz_data.lat.values,
                     'TIDAL_Kz': Kz_data.TIDAL_Kz.values[0, :, :, :]}
        # Setting all values not in the Mediterranean to np.nan
        Lat, _, Lon = np.meshgrid(data_dict['LAT'], data_dict['depth'], data_dict['LON'])
        data_dict['TIDAL_Kz'][Lat < self.lat_grid.min()] = np.nan
        data_dict['TIDAL_Kz'][Lat > self.lat_grid.max()] = np.nan
        data_dict['TIDAL_Kz'][Lon < self.lon_grid.min()] = np.nan
        data_dict['TIDAL_Kz'][Lon > self.lon_grid.max()] = np.nan
        return data_dict