import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import cmocean.cm as cmo


class SizeTransport_peak_depth:
    def __init__(self, scenario, figure_direc, size, time_selection, rho=920):
        # Figure Parameters
        self.figure_size = (20, 20)
        self.figure_shape = (2, 2)
        self.ax_label_size = 18
        self.ax_ticklabel_size = 16
        self.number_of_plots = 4
        self.adv_file_dict = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                            advection_scenario='CMEMS_MEDITERRANEAN',
                                                            repeat_dt=None).file_names
        self.spatial_domain = np.nanmin(self.adv_file_dict['LON']),  np.nanmax(self.adv_file_dict['LON']), \
                              np.nanmin(self.adv_file_dict['LAT']), np.nanmax(self.adv_file_dict['LAT'])
        self.cmap = cmo.speed
        # Data parameters
        self.output_direc = figure_direc + 'vertical_profile/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'spatial_vertical_profile'
        # Simulation parameters
        self.scenario = scenario
        self.size = size
        self.time_selection = time_selection
        self.year = 2010 + self.time_selection
        self.rho = rho
        self.tau = 0.0

    def plot(self):
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
        norm = colors.Normalize(vmin=1, vmax=100)
        cbar_label, extend = r"Relative Concentration ($C/C_{min}$)", 'max'
        cmap = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        cax = fig.add_subplot(gs[:, -1])
        cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend=extend)
        cbar.set_label(cbar_label, fontsize=self.ax_label_size)
        cbar.ax.tick_params(which='major', labelsize=self.ax_ticklabel_size, length=14, width=2)
        cbar.ax.tick_params(which='minor', labelsize=self.ax_ticklabel_size, length=7, width=2)

        # Saving the figure
        file_name = self.plot_save_name()
        plt.savefig(file_name, bbox_inches='tight')
        plt.close('all')

    def plot_save_name(self, file_type='.png'):
        str_format = self.size, self.rho, self.time_selection
        return self.output_direc + 'VerticalPeak_size={:.2E}_rho={}_year={}'.format(*str_format) + file_type
