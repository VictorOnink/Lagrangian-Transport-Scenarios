import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string
import cmocean.cm as cmo


class FragmentationKaandorpPartial_lonlat_averages:
    def __init__(self, scenario, figure_direc, lambda_frag, time_selection, beach_state, mass, rho=920, shore_time=20):
        # Simulation parameters
        self.scenario = scenario
        self.rho = rho
        self.time_selection = time_selection
        self.size = settings.INIT_SIZE
        self.lambda_frag = lambda_frag
        self.beach_state = beach_state
        self.beach_state_list = ['beach', 'adrift']
        self.dimension_list = ["lon_counts", "lat_counts"]
        self.weight = {False: 'particle_number_sink', True: 'particle_mass_sink'}[mass]
        self.total_list = ['total_lon_{}'.format(self.weight), 'total_lat_{}'.format(self.weight)]
        self.size_class_list = np.arange(0, settings.SIZE_CLASS_NUMBER)
        self.key_concentration = utils.analysis_simulation_year_key(self.time_selection)
        self.step = 0.5
        self.shore_time = shore_time
        # Data parameters
        self.output_direc = figure_direc + 'concentrations/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/FragmentationKaandorpPartial/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'lonlat_concentration'
        # Figure parameters
        self.figure_size = (20, 10)
        self.figure_shape = (2, 2)
        self.ax_label_size = 14
        self.ax_ticklabel_size = 14
        self.adv_file_dict = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                            advection_scenario='CMEMS_MEDITERRANEAN',
                                                            repeat_dt=None).file_names
        self.spatial_domain = np.nanmin(self.adv_file_dict['LON']), np.nanmax(self.adv_file_dict['LON']), \
                              np.nanmin(self.adv_file_dict['LAT']), np.nanmax(self.adv_file_dict['LAT'])
        self.cmap = 'viridis_r'
        self.line_types = {'beach': '-', 'adrift': '-'}

    def plot(self):
        # Loading the data
        data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=self.scenario, prefix=self.prefix,
                                                                  data_direc=self.data_direc,
                                                                  shore_time=self.shore_time,
                                                                  lambda_frag=self.lambda_frag, rho=self.rho,
                                                                  postprocess=True)
        concentration_dict, lon, lat = self.histogram_reduction(data_dict=data_dict)

        # Normalizing the concentrations by the total number of particles in the simulation
        for beach_state in self.beach_state_list:
            for size_class in self.size_class_list:
                concentration_dict[beach_state][size_class]["lon_counts"] /= concentration_dict[self.total_list[0]]
                concentration_dict[beach_state][size_class]["lat_counts"] /= concentration_dict[self.total_list[1]]

        # Creating the map
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1], width_ratios=[1, 0.5],
                              height_ratios=[1, 0.5], wspace=0.1, hspace=0.1)

        ax_map = vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=0, column=0, add_gridlabels=False,
                                             domain=self.spatial_domain, label_size=self.ax_label_size,
                                             lat_grid_step=2, lon_grid_step=5, resolution='10m', land_color='grey',
                                             border_color='white',
                                             x_grid_locator=np.arange(start=-5, stop=40, step=5))
        ax_map.set_title("{} - {}".format(self.beach_state, self.weight), fontsize=self.ax_label_size + 2,
                         fontweight='bold')

        # Creating the axis for the longitudes
        ax_lon = fig.add_subplot(gs[1, 0])
        ax_lon.set_xlim((self.spatial_domain[0], self.spatial_domain[1]))
        ax_lon.grid(which='major', axis='x', linestyle='-')
        ax_lon.set_xlabel(r'Longitude ($^{\circ}$)', fontsize=self.ax_label_size)
        ax_lon.set_ylabel(r'Particle Fraction', fontsize=self.ax_label_size)
        ax_lon.tick_params(axis='both', labelsize=self.ax_ticklabel_size)
        ax_lon.set_yscale('log')
        ax_lon.set_ylim((1e-5, 1e0))

        # Creating the axis for the latitudes
        ax_lat = fig.add_subplot(gs[0, 1])
        ax_lat.set_ylim((self.spatial_domain[2], self.spatial_domain[3]))
        ax_lat.grid(which='major', axis='y', linestyle='-')
        ax_lat.yaxis.set_label_position("right")
        ax_lat.yaxis.tick_right()
        ax_lat.xaxis.set_label_position("top")
        ax_lat.xaxis.tick_top()
        ax_lat.set_xlabel(r'Particle Fraction', fontsize=self.ax_label_size)
        ax_lat.set_ylabel(r'Latitude ($^{\circ}$)', fontsize=self.ax_label_size)
        ax_lat.tick_params(axis='both', labelsize=self.ax_ticklabel_size)
        ax_lat.set_xscale('log')
        ax_lat.set_xlim((1e-5, 1e0))

        # Adding a legend
        ax_legend = fig.add_subplot(gs[1, 1])
        size_colors = [plt.plot([], [], c=vUtils.discrete_color_from_cmap(index_size, subdivisions=settings.SIZE_CLASS_NUMBER),
                                label=size_label(size_class), linestyle='-')[0] for index_size, size_class in enumerate(self.size_class_list)]
        ax_legend.legend(handles=size_colors, fontsize=self.ax_label_size, loc='upper left', ncol=1)
        ax_legend.axis('off')

        # Plotting the longitudes
        for index_size, size_class in enumerate(self.size_class_list):
            ax_lon.plot(lon, concentration_dict[self.beach_state][size_class]["lon_counts"],
                        linestyle=self.line_types[self.beach_state],
                        c=vUtils.discrete_color_from_cmap(index_size, subdivisions=settings.SIZE_CLASS_NUMBER))

        # Plotting the latitudes
        for index_size, size_class in enumerate(self.size_class_list):
            ax_lat.plot(concentration_dict[self.beach_state][size_class]["lat_counts"], lat,
                        linestyle=self.line_types[self.beach_state],
                        c=vUtils.discrete_color_from_cmap(index_size, subdivisions=settings.SIZE_CLASS_NUMBER))

        # Saving the figure
        str_format = self.time_selection, self.beach_state, self.rho, self.lambda_frag, self.weight
        fig_name = self.output_direc + "lon_lat_year={}_beach_state={}_rho={}_lambda_f={}_weight={}.png".format(*str_format)
        plt.savefig(fig_name, bbox_inches='tight')

    def histogram_reduction(self, data_dict):
        # The new bin edges
        bins_lon = np.arange(self.spatial_domain[0], self.spatial_domain[1] + self.step, step=self.step)
        bins_lat = np.arange(self.spatial_domain[2], self.spatial_domain[3] + self.step, step=self.step)
        # The output dictionary with the new concentrations
        output_dict = {}
        # Selecting just the specified year of the data
        coordinate = {"lon_counts": data_dict['lon'], "lat_counts": data_dict['lat']}
        bin = {"lon_counts": bins_lon, "lat_counts": bins_lat}
        year_data = data_dict[self.key_concentration]

        for beach_state in self.beach_state_list:
            output_dict[beach_state] = {}
            for size_class in self.size_class_list:
                output_dict[beach_state][size_class] = {}
                for lonlat in self.dimension_list:
                    output_dict[beach_state][size_class][lonlat], _ = np.histogram(a=coordinate[lonlat],
                                                                                   bins=bin[lonlat],
                                                                                   weights=year_data[beach_state][self.weight][size_class][lonlat])
        for key in self.total_list:
            output_dict[key] = year_data[key]
        # Calculate the bin edge midpoints
        bin_mid_lon = 0.5 * bins_lon[1:] + 0.5 * bins_lon[:-1]
        bin_mid_lat = 0.5 * bins_lat[1:] + 0.5 * bins_lat[:-1]

        return output_dict, bin_mid_lon, bin_mid_lat


def size_label(size_class):
    particle_size = settings.INIT_SIZE * settings.P_FRAG ** size_class
    return 'Size class {}, d = {:.2f} mm'.format(size_class, particle_size * 1e3)


