import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
from datetime import datetime, timedelta
import matplotlib.animation as animation
import cmocean.cm as cmo
import string


class FragmentationKaandorpPartial_Animation:
    def __init__(self, scenario, figure_direc, shore_time, rho, simulation_years, ocean_frag=False):
        # Figure parameters
        self.figure_size = (20, 10)
        self.ax_label_size = 18
        self.tick_label_size = 16
        self.grid_shape = (2, 3)
        self.adv_file_dict = advection_files.AdvectionFiles(advection_scenario='CMEMS_MEDITERRANEAN').file_names
        self.spatial_domain = np.nanmin(self.adv_file_dict['LON']),  np.nanmax(self.adv_file_dict['LON']), \
                              np.nanmin(self.adv_file_dict['LAT']), np.nanmax(self.adv_file_dict['LAT'])
        self.cmap = cmo.haline_r
        self.fps = 10
        # Data parameters
        self.output_direc = figure_direc + 'animations/'
        self.data_direc = settings.DATA_OUTPUT_DIREC + 'timeslices/FragmentationKaandorpPartial/'
        utils.check_direc_exist(self.output_direc)
        self.ocean_frag = ocean_frag
        # Simulation parameters
        self.scenario = scenario
        self.shore_time = shore_time
        self.lambda_frag = 388
        self.rho = rho
        self.simulation_years = simulation_years

    def animate(self):
        # Creating the base figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.grid_shape[0], ncols=self.grid_shape[1] + 1, width_ratios=[1, 1, 1, 0.1])
        ax_list = []
        for rows in range(self.grid_shape[0]):
            for columns in range(self.grid_shape[1]):
                ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                           domain=self.spatial_domain, label_size=self.tick_label_size,
                                                           lat_grid_step=5, lon_grid_step=10, resolution='10m',
                                                           ocean_color='black', border_color='white'))
        # Setting the colormap for the particle depth
        norm = colors.Normalize(vmin=0.0, vmax=100.0)
        cmap = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        cax = fig.add_subplot(gs[:, -1])
        cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend='max')
        cbar.set_label(r"Depth (m)", fontsize=self.ax_label_size)
        cbar.ax.tick_params(which='major', labelsize=self.tick_label_size, length=14, width=2)
        cbar.ax.tick_params(which='minor', labelsize=self.tick_label_size, length=7, width=2)

        # Setting the time range for which we want to create the simulation
        current_time, end_time = datetime(2010, 1, 1, 0), datetime(2010 + self.simulation_years, 1, 1, 0)
        time_step, time_list = timedelta(hours=24), []
        while current_time < end_time:
            time_list.append(current_time)
            current_time += time_step
        frame_number = len(time_list)

        # Setting a text box for the simulation date
        props = dict(boxstyle='round', facecolor='white', alpha=1)
        text = ax_list[3].text(0.02, 0.02, 'initial', horizontalalignment='left', verticalalignment='bottom',
                               transform=ax_list[3].transAxes, bbox=props, fontsize=self.ax_label_size, zorder=200)

        # Setting the initial values of the x and y coordinates, which will later be filled by lon and lat
        plot_list = []
        for ax_index, ax in enumerate(ax_list):
            plot_list.append(ax.scatter(0, 0, c=0, s=7, alpha=1, zorder=1000, cmap=self.cmap, norm=norm))
            ax.set_title(subfigure_title(ax_index), fontsize=self.ax_label_size)

        # Initialization function
        def init():
            for plot in plot_list:
                plot.set_offsets(np.c_[[], []])
            text.set_text('initial 2')
            return plot_list

        # Animation function
        def animate_function(frame_index):
            utils.print_statement("we are at index {} of {}".format(frame_index, frame_number), to_print=True)
            date = time_list[frame_index].strftime("%Y-%m-%d-%H-%M-%S")
            # Loading the data
            prefix = 'timeslices_{}'.format(date)
            data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=self.scenario, prefix=prefix,
                                                                      data_direc=self.data_direc,
                                                                      lambda_frag=self.lambda_frag,
                                                                      rho=self.rho, shore_time=self.shore_time,
                                                                      ocean_frag=self.ocean_frag, postprocess=True)
            lon, lat, depth = data_dict['lon'], data_dict['lat'], data_dict['z'].astype(int)
            size_class = data_dict['size_class']
            # Looping through the axes corresponding to the different size classes
            for ax_index, size in enumerate(ax_list):
                size_selection = size_class == ax_index
                lon_select, lat_select, depth_select = lon[size_selection], lat[size_selection], depth[size_selection]
                # Updating the plot on each axis with the data
                plot_list[ax_index].set_offsets(np.c_[lon_select, lat_select])
                plot_list[ax_index].set_array(depth_select)
            text.set_text(time_list[frame_index].strftime("%Y-%m-%d"))
            return plot_list

        # Calling the animator
        animator = animation.FuncAnimation(plt.gcf(), animate_function, init_func=init, frames=frame_number,
                                           interval=100, blit=True)

        # Saving the animation
        animator.save(filename=animation_save_name(output_direc=self.output_direc, shore_time=self.shore_time,
                                                   ocean_frag=self.ocean_frag),
                      fps=self.fps, extra_args=['-vcodec', 'libx264'])


def animation_save_name(output_direc, shore_time, ocean_frag, flowdata='CMEMS_MEDITERRANEAN', file_type='.mov'):
    return output_direc + 'FragmentationKaandorpPartial_OFRAG_{}_{}_st={}'.format(ocean_frag, flowdata, shore_time) + \
           file_type


def subfigure_title(index):
    size = settings.INIT_SIZE * 2 ** (-1 * index) * 1e3
    return '({}) size class {}, d = {:.3f} mm'.format(string.ascii_letters[index], index, size)