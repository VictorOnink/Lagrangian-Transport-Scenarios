import settings
import utils
from netCDF4 import Dataset
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean.cm as cmo


class General_bathymetry_histogram:
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
        depth = dataset.variables[self.scenario.file_dict['BATH_variables']['DEPTH']][:]

        # Flattening the array, and removing all depths == 0
        depth = depth.flatten()
        depth = depth[depth > 0]

        # Set the depth bins, and then calculate a histogram of the depths
        depth_bins = np.logspace(0, 3.5)
        depth_bins_mid = (depth_bins[:-1] + depth_bins[1:]) / 2
        histogram_depths, _ = np.histogram(depth, bins=depth_bins)

        # Normalize all the depth bins by the total number of cells
        histogram_depths = np.divide(histogram_depths, np.sum(histogram_depths))
        histogram_depths *= 100

        # Creating the figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1])

        ax = fig.add_subplot(gs[0, 0])
        ax.set_xlim([0, 20])
        ax.set_yscale('symlog')
        ax.set_ylim([-3000, -1])

        ax.plot(histogram_depths, -1 * depth_bins_mid)

        ax.set_aspect('auto', adjustable=None)
        file_name = self.output_direc + 'Bathymetry_histogram.png'
        plt.savefig(file_name, bbox_inches='tight')
