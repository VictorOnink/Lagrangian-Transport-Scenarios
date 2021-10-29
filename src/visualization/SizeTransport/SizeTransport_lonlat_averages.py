import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string
import cmocean.cm as cmo


class SizeTransport_lonlat_averages:
    def __init__(self, scenario, figure_direc, size_list, time_selection, beach_state, rho=920, tau=0):
        # Simulation parameters
        self.scenario = scenario
        self.rho = rho
        self.time_selection = time_selection
        self.size_list = size_list
        self.tau = tau
        self.beach_state = beach_state
        self.beach_state_list = ['beach', 'adrift']
        self.dimension_list = ["lon_counts", "lat_counts"]
        self.key_concentration = utils.analysis_simulation_year_key(self.time_selection)
        self.step = 1
        # Data parameters
        self.output_direc = figure_direc + 'concentrations/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'lonlat_concentration'
        # Figure parameters
        self.figure_size = (20, 10)
        self.figure_shape = (2, 2)
        self.ax_label_size = 14
        self.ax_ticklabel_size = 14
        self.number_of_plots = self.size_list.__len__()
        self.adv_file_dict = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                            advection_scenario='CMEMS_MEDITERRANEAN',
                                                            repeat_dt=None).file_names
        self.spatial_domain = np.nanmin(self.adv_file_dict['LON']), np.nanmax(self.adv_file_dict['LON']), \
                              np.nanmin(self.adv_file_dict['LAT']), np.nanmax(self.adv_file_dict['LAT'])
        self.cmap = 'viridis_r'
        self.line_types = {'beach': '-', 'adrift': '-'}

    def plot(self):
        # Loading the data
        concentration_dict = {}
        key_concentration = utils.analysis_simulation_year_key(self.time_selection)
        for index, size in enumerate(self.size_list):
            data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix,
                                                       data_direc=self.data_direc,
                                                       size=size, rho=self.rho, tau=self.tau)
            concentration_dict[index], lon, lat = self.histogram_reduction(data_dict=data_dict)

        # Normalizing the concentrations by the total number of particles in the simulation
        for size in concentration_dict.keys():
            for beach_state in self.beach_state_list:
                for lonlat in self.dimension_list:
                    concentration_dict[size][beach_state][lonlat] /= np.nansum(
                        concentration_dict[size][beach_state][lonlat])

        # Setting zero values to NAN
        for size in concentration_dict.keys():
            for beach_state in self.beach_state_list:
                for lonlat in self.dimension_list:
                    zero = concentration_dict[size][beach_state][lonlat] == 0
                    concentration_dict[size][beach_state][lonlat][zero] = np.nan

        # Creating the map
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1], width_ratios=[1, 0.5],
                              height_ratios=[1, 0.5], wspace=0.1, hspace=0.1)

        ax_map = vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=0, column=0, add_gridlabels=False,
                                             domain=self.spatial_domain, label_size=self.ax_label_size,
                                             lat_grid_step=2, lon_grid_step=5, resolution='10m', land_color='grey',
                                             border_color='white',
                                             x_grid_locator=np.arange(start=-5, stop=40, step=5))
        ax_map.set_title(self.beach_state, fontsize=self.ax_label_size + 2, fontweight='bold')

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
        # beach_lines = [plt.plot([], [], c='k', label=beach_state, linestyle=self.line_types[beach_state])[0]
        #                for beach_state in self.beach_state_list]

        size_colors = [plt.plot([], [], c=vUtils.discrete_color_from_cmap(index_size, subdivisions=self.size_list.__len__()),
                                label=size_label(size), linestyle='-')[0] for index_size, size in enumerate(self.size_list)]
        ax_legend.legend(handles=size_colors, fontsize=self.ax_label_size, loc='upper left', ncol=2)
        ax_legend.axis('off')

        # Plotting the longitudes
        for index_size, size in enumerate(self.size_list):
            ax_lon.plot(lon, concentration_dict[index_size][self.beach_state]["lon_counts"],
                        linestyle=self.line_types[self.beach_state],
                        c=vUtils.discrete_color_from_cmap(index_size, subdivisions=self.size_list.__len__()))

        # Plotting the latitudes
        for index_size, size in enumerate(self.size_list):
            ax_lat.plot(concentration_dict[index_size][self.beach_state]["lat_counts"], lat,
                        linestyle=self.line_types[self.beach_state],
                        c=vUtils.discrete_color_from_cmap(index_size, subdivisions=self.size_list.__len__()))

        # Saving the figure
        str_format = self.time_selection, self.beach_state, self.rho
        fig_name = self.output_direc + "lon_lat_year={}_beach_state={}_rho={}.png".format(*str_format)
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
        data_dict = data_dict[self.key_concentration]
        for size in data_dict.keys():
            output_dict[size] = {}
            for beach_state in self.beach_state_list:
                output_dict[size][beach_state] = {}
                for lonlat in self.dimension_list:
                    output_dict[size][beach_state][lonlat], _ = np.histogram(a=coordinate[lonlat], bins=bin[lonlat],
                                                                             weights=data_dict[size][beach_state][lonlat])
        # Calculate the bin edge midpointw
        bin_mid_lon = 0.5 * bins_lon[1:] + 0.5 * bins_lon[:-1]
        bin_mid_lat = 0.5 * bins_lat[1:] + 0.5 * bins_lat[:-1]

        return output_dict, bin_mid_lon, bin_mid_lat


def size_label(size):
    return r'r = {:.3f} mm'.format(size * 1e3)



