import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string
import cmocean.cm as cmo


class FragmentationKaandorpPartial_Concentration:
    def __init__(self, scenario, figure_direc, rho, shore_time, beach_state, simulation_year, lambda_frag, mass,
                 sink=True, input='LebretonDivision'):
        # Simulation parameters
        self.scenario = scenario
        self.shore_time = shore_time
        self.rho = rho
        self.beach_state = beach_state
        self.simulation_year = simulation_year
        self.lambda_frag = lambda_frag
        self.mass = mass
        self.weight = {False: 'particle_number', True: 'particle_mass'}[self.mass]
        self.sink = sink
        if self.sink:
            self.weight += '_sink'
        self.input = input
        # Data parameters
        self.output_direc = figure_direc + 'concentrations/'
        self.data_direc = settings.DATA_OUTPUT_DIREC + 'concentrations/FragmentationKaandorpPartial/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'horizontal_concentration'
        self.key_concentration = utils.analysis_simulation_year_key(self.simulation_year)
        # Figure parameters
        self.figure_size = (20, 10)
        self.figure_shape = (3, 2)
        self.ax_label_size = 14
        self.ax_ticklabel_size = 14
        self.number_of_plots = self.figure_shape[0] * self.figure_shape[1]
        self.adv_file_dict = advection_files.AdvectionFiles(advection_scenario='CMEMS_MEDITERRANEAN').file_names
        self.spatial_domain = np.nanmin(self.adv_file_dict['LON']),  np.nanmax(self.adv_file_dict['LON']), \
                              np.nanmin(self.adv_file_dict['LAT']), np.nanmax(self.adv_file_dict['LAT'])
        self.cmap = cmo.thermal

    def plot(self):
        # Loading the data
        data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=self.scenario, prefix=self.prefix,
                                                                  data_direc=self.data_direc, shore_time=self.shore_time,
                                                                  lambda_frag=self.lambda_frag, rho=self.rho,
                                                                  postprocess=True, input=self.input)
        concentration_dict = data_dict[self.key_concentration][self.beach_state][self.weight]
        Lon, Lat = np.meshgrid(data_dict['lon'], data_dict['lat'])

        # Normalizing the concentrations by the lowest non-zero concentration over all the sizes
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

        # Creating the figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1] + 1,
                              width_ratios=[1] * self.figure_shape[1] + [0.1])

        ax_list = []
        for rows in range(self.figure_shape[0]):
            for columns in range(self.figure_shape[1]):
                ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                           domain=self.spatial_domain, lat_grid_step=5,
                                                           lon_grid_step=10, resolution='10m', auto_aspect=False))
        # Setting the colormap and creating the colorbar
        norm = set_normalization(self.beach_state, self.simulation_year, self.weight)
        cbar_label = "Relative {}".format({True: 'Mass', False: 'Count'}[self.mass]) + r" Concentration ($C/C_{min}$)"
        extend = 'both'
        cmap = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        cax = fig.add_subplot(gs[:, -1])
        cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend=extend)
        cbar.set_label(cbar_label, fontsize=self.ax_label_size)
        cbar.ax.tick_params(which='major', labelsize=self.ax_ticklabel_size, length=14, width=2)
        cbar.ax.tick_params(which='minor', labelsize=self.ax_ticklabel_size, length=7, width=2)

        # Adding subfigure titles
        for index, ax in enumerate(ax_list):
            ax.set_title(subfigure_title(index), fontsize=self.ax_label_size)

        # The actual plotting of the figures
        for index, ax in enumerate(ax_list):
            if self.beach_state in ['adrift']:
                ax.pcolormesh(Lon, Lat, concentration_dict[index], norm=norm, cmap=self.cmap, zorder=200)
            else:
                ax.scatter(Lon.flatten(), Lat.flatten(), c=concentration_dict[index], norm=norm, cmap=self.cmap,
                           zorder=200)

        plt.savefig(self.file_name(), bbox_inches='tight')
        plt.close()

    def file_name(self):
        str_format = self.input, self.beach_state, self.simulation_year, self.lambda_frag, self.shore_time, self.rho, \
                     self.weight, self.sink
        prefix = self.output_direc + 'Concentrations_'
        return prefix + '{}_{}_year={}_lambda_f={}_st={}_rho={}_{}_sink={}.png'.format(*str_format)


def subfigure_title(index):
    """
    setting the title of the subfigure
    :param index:
    :param size:
    :param rho:
    :return:
    """
    particle_size = utils.size_range(single_size_class=index, units='mm')
    title = '({}) k = {}, d = {:.3f} mm'.format(string.ascii_lowercase[index], index, particle_size)
    return title


def set_normalization(beach_state, simulation_year, weight):
    """
    Setting the normalization that we use for the colormap
    :param beach_state: adrift, beach or seabed
    :param difference: is it absolute concentrations or the difference relative to reference size
    :return:
    """
    if beach_state == 'adrift':
        if simulation_year == 0:
            if weight in ['particle_number', 'particle_number_sink']:
                vmin, vmax = 1e1, 1e5
            elif weight in ['particle_mass', 'particle_mass_sink']:
                vmin, vmax = 1e2, 1e9
        elif simulation_year == 1:
            if weight in ['particle_number', 'particle_number_sink']:
                vmin, vmax = 1e1, 1e5
            elif weight in ['particle_mass', 'particle_mass_sink']:
                vmin, vmax = 1e2, 1e9
        elif simulation_year == 2:
            if weight in ['particle_number', 'particle_number_sink']:
                vmin, vmax = 1e1, 1e5
            elif weight in ['particle_mass', 'particle_mass_sink']:
                vmin, vmax = 1e2, 1e9
    elif beach_state == 'beach':
        if simulation_year == 0:
            vmin, vmax = 1e6, 1e15
        elif simulation_year == 1:
            vmin, vmax = 1e15, 1e24
    return colors.LogNorm(vmin=vmin, vmax=vmax)


