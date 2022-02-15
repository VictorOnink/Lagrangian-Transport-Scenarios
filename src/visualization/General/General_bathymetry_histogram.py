import settings
import utils
from netCDF4 import Dataset
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean.cm as cmo


class General_bathymetry_histogram:
    def __init__(self, scenario, figure_direc, depth_selection='all'):
        # Scenario specific variables
        self.scenario = scenario
        self.figure_direc = figure_direc
        self.depth_selection = depth_selection
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
        depth = dataset.variables[self.scenario.file_dict['BATH_variables']['DEPTH']][:]
        depth[depth <= 0] = np.nan

        # Filter out the cells in the nearshore and offshore
        distance2shore = Dataset(self.scenario.file_dict['DISTANCE_filename']).variables['distance'][:]
        print('max {} min {}'.format(np.nanmax(distance2shore), np.nanmin(distance2shore)))
        if self.depth_selection in ['offshore']:
            depth[distance2shore < 50] = np.nan
            print('In the offshore case, we have {} cells'.format(np.sum(~np.isnan(depth))))
        elif self.depth_selection in ['nearshore']:
            depth[distance2shore > 50] = np.nan
            print('In the nearshore case, we have {} cells'.format(np.sum(~np.isnan(depth))))
        elif self.depth_selection in ['coastal']:
            depth[distance2shore > 10] = np.nan
            print('In the coastal case, we have {} cells'.format(np.sum(~np.isnan(depth))))
        else:
            print('We have {}'.format(self.depth_selection))

        # Flattening the array, and removing all depths == 0
        depth = depth.flatten()
        depth = depth[~np.isnan(depth)]

        # Set the depth bins, and then calculate a histogram of the depths
        depth_bins = np.logspace(0, 3.5)
        depth_bins_mid = (depth_bins[:-1] + depth_bins[1:]) / 2
        histogram_depths, _ = np.histogram(depth, bins=depth_bins)

        # Normalize all the depth bins by the total number of cells
        histogram_depths = np.divide(histogram_depths, np.nansum(histogram_depths))
        histogram_depths *= 100

        # Creating the figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1])

        ax = fig.add_subplot(gs[0, 0])
        ax.set_yscale('symlog')
        ax.set_ylim([-3000, -1])

        ax.set_ylabel('Depth (m)')
        ax.set_xlabel('Percentage of all ocean cells')

        ax.plot(histogram_depths, -1 * depth_bins_mid)

        ax.set_aspect('auto', adjustable=None)
        file_name = self.output_direc + 'Bathymetry_histogram_{}.png'.format(self.depth_selection)
        plt.savefig(file_name, bbox_inches='tight')
