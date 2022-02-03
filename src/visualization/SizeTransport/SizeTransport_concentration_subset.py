import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string
import cmocean.cm as cmo


class SizeTransport_concentration_subset:
    def __init__(self, scenario, figure_direc, time_selection,  size_list, rho=920, tau=0,
                 fixed_resus=False, resus_time=50):
        # Simulation parameters
        self.scenario = scenario
        self.rho = rho
        self.time_selection = time_selection
        self.beach_state = 'adrift'
        self.size_list = size_list
        self.tau = tau
        self.depth_list = ['column', 'surface_1m']
        self.fixed_resus = fixed_resus
        self.resus_time = resus_time
        # Data parameters
        self.output_direc = figure_direc + 'concentrations/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'horizontal_concentration'
        # Figure parameters
        self.figure_size = (20, 20)
        self.figure_shape = (self.size_list.__len__(), 2)
        self.ax_label_size = 18
        self.ax_ticklabel_size = 16
        self.number_of_plots = self.size_list.__len__() * 2
        self.adv_file_dict = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                            advection_scenario='CMEMS_MEDITERRANEAN',
                                                            repeat_dt=None).file_names
        self.spatial_domain = np.nanmin(self.adv_file_dict['LON']),  np.nanmax(self.adv_file_dict['LON']), \
                              np.nanmin(self.adv_file_dict['LAT']), np.nanmax(self.adv_file_dict['LAT'])
        self.cmap = cmo.thermal

    def plot(self):
        # Loading the data
        concentration_dict = {'column': {}, 'surface_1m': {}}
        key_concentration = utils.analysis_simulation_year_key(self.time_selection)
        for index, size in enumerate(self.size_list):
            for depth_level in concentration_dict.keys():
                data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix,
                                                           data_direc=self.data_direc, fixed_resus=self.fixed_resus,
                                                           size=size, rho=self.rho, tau=self.tau,
                                                           resus_time=self.resus_time)
                if self.beach_state in ['adrift']:
                    concentration_array = data_dict[key_concentration]['adrift'][depth_level]
                concentration_dict[depth_level][index] = concentration_array
        Lon, Lat = np.meshgrid(data_dict['lon'], data_dict['lat'])

        # Normalizing the concentration by the lowest non-zero concentration over all the sizes
        normalization_factor = 1e10
        for depth_level in concentration_dict.keys():
            for size in concentration_dict[depth_level].keys():
                concentration = concentration_dict[depth_level][size]
                min_non_zero = np.nanmin(concentration[concentration > 0])
                if min_non_zero < normalization_factor:
                    normalization_factor = min_non_zero
        for depth_level in concentration_dict.keys():
            for size in concentration_dict[depth_level].keys():
                concentration_dict[depth_level][size] /= normalization_factor

        # Setting zero values to nan
        for beach_state in concentration_dict.keys():
            for size in concentration_dict[beach_state].keys():
                concentration_dict[beach_state][size][concentration_dict[beach_state][size] == 0] = np.nan

        # Creating the base figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1] + 1,
                              width_ratios=[1, 1, 0.1])

        ax_list = []
        for rows in range(self.figure_shape[0]):
            for columns in range(self.figure_shape[1]):
                ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                           domain=self.spatial_domain,
                                                           label_size=self.ax_label_size,
                                                           lat_grid_step=5, lon_grid_step=10, resolution='10m'))

        # Setting the colormap, and adding a colorbar
        norm = self.set_normalization()
        cbar_label, extend = r"Relative Concentration ($C/C_{min}$)", 'max'
        cmap = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        cax = fig.add_subplot(gs[:, -1])
        cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend=extend)
        cbar.set_label(cbar_label, fontsize=self.ax_label_size)
        cbar.ax.tick_params(which='major', labelsize=self.ax_ticklabel_size, length=14, width=2)
        cbar.ax.tick_params(which='minor', labelsize=self.ax_ticklabel_size, length=7, width=2)

        # Adding subfigure titles
        for index, ax in enumerate(ax_list):
            ax.set_title(self.subfigure_title(index), weight='bold', fontsize=self.ax_label_size)

        # The actual plotting of the figures
        for index, size in enumerate(self.size_list):
            if self.beach_state in ['adrift']:
                ax_list[index].pcolormesh(Lon, Lat, concentration_dict['adrift'][index], norm=norm, cmap=self.cmap,
                                          zorder=200)
            else:
                ax_list[index].scatter(Lon.flatten(), Lat.flatten(), c=concentration_dict['beach'][index].flatten(),
                                       norm=norm, cmap=self.cmap, zorder=200)

        # Saving the figure
        print(self.plot_save_name())
        plt.savefig(self.plot_save_name(), bbox_inches='tight')
        plt.close('all')

    def plot_save_name(self, file_type='.png'):
        year = {0: 'year_0', 1: 'year_1', 2: 'year_2'}[self.time_selection]
        str_format = self.rho, self.rho, year
        name = self.output_direc + 'rho_{}/SizeTransport_rho={}_subset_year={}'.format(*str_format)
        if self.fixed_resus:
            name += '_fixed_resus_{}'.format(self.resus_time)
        return name + file_type

    @staticmethod
    def set_normalization():
        """
        Setting the normalization that we use for the colormap
        :param beach_state: adrift, beach or seabed
        :return:
        """
        return colors.LogNorm(vmin=1e0, vmax=1e4)

    def subfigure_title(self, index):
        return '({}) r = {:.3f} mm'.format(string.ascii_lowercase[index], self.size_list[index // 2] * 1e3)





