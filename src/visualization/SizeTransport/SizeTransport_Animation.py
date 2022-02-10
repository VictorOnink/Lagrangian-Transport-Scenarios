import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
from datetime import datetime, timedelta
import string
import matplotlib.animation as animation
import cmocean.cm as cmo


class SizeTransport_Animation:
    def __init__(self, scenario, figure_direc, size_list, simulation_years, rho=920, tau=0.0):
        # Figure parameters
        self.figure_size = (20, 10)
        self.ax_label_size = 16
        self.tick_label_size = 16
        self.grid_shape = (2, 2)
        self.adv_file_dict = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                            advection_scenario='CMEMS_MEDITERRANEAN',
                                                            repeat_dt=None).file_names
        self.spatial_domain = np.nanmin(self.adv_file_dict['LON']),  np.nanmax(self.adv_file_dict['LON']), \
                              np.nanmin(self.adv_file_dict['LAT']), np.nanmax(self.adv_file_dict['LAT'])
        self.cmap = cmo.haline_r
        self.fps = 10
        # Data parameters
        self.output_direc = figure_direc + 'animations/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeslices/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        # Simulation parameters
        self.scenario = scenario
        self.size_list = size_list
        self.lambda_frag = 388
        self.rho = rho
        self.tau = tau
        self.simulation_years = simulation_years

    def animate(self):
        # Creating the base figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.grid_shape[0], ncols=self.grid_shape[1] + 1, width_ratios=[1, 1, 0.1])
        ax_list = []
        for rows in range(self.grid_shape[0]):
            for columns in range(self.grid_shape[1]):
                ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                           domain=self.spatial_domain, label_size=self.tick_label_size,
                                                           lat_grid_step=5, lon_grid_step=10, resolution='10m',
                                                           ocean_color='black', border_color='white',
                                                           add_gridlines=False))
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
        frame_number = time_list.__len__()

        # Setting a text box for the simulation date
        props = dict(boxstyle='round', facecolor='white', alpha=1)
        text = ax_list[2].text(0.02, 0.02, 'initial', horizontalalignment='left', verticalalignment='bottom',
                               transform=ax_list[3].transAxes, bbox=props, fontsize=self.ax_label_size, zorder=200)

        # Setting the initial values of the x and y coordinates, which will later be filled by lon and lat
        plot_list = []
        for ax_index, ax in enumerate(ax_list):
            plot_list.append(ax.scatter(0, 0, c=0, s=7, alpha=1, zorder=1000, cmap=self.cmap, norm=norm))
            ax.set_title(subfigure_title(ax_index, self.size_list[ax_index], self.rho), fontsize=self.ax_label_size)

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
            for index, size in enumerate(self.size_list):
                # Loading the dictionary with the data
                prefix = 'timeslices_{}'.format(date)
                data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=prefix, data_direc=self.data_direc,
                                                           size=size, rho=self.rho, tau=self.tau)
                lon, lat, depth = data_dict['lon'], data_dict['lat'], data_dict['z'].astype(int)
                # Updating the plot on each axis with the data
                plot_list[index].set_offsets(np.c_[lon, lat])
                plot_list[index].set_array(depth)
            text.set_text(time_list[frame_index].strftime("%Y-%m-%d"))
            return plot_list

        # Calling the animator
        animator = animation.FuncAnimation(plt.gcf(), animate_function, init_func=init, frames=frame_number,
                                           interval=100, blit=True)

        # Saving the animation
        animator.save(filename=animation_save_name(output_direc=self.output_direc, rho=self.rho),
                      fps=self.fps, extra_args=['-vcodec', 'libx264'])


def animation_save_name(output_direc, rho, flowdata='CMEMS_MEDITERRANEAN', startyear=2010):
    return output_direc + 'SizeTransport_{}_rho_{}_y_{}.mov'.format(flowdata, rho, startyear)


def subfigure_title(index, size, rho):
    alphabet = string.ascii_lowercase
    return '({}) r = {:.3f} mm, '.format(alphabet[index], size * 1e3) + r'$\rho$ = ' + '{} kg m'.format(rho) + r'$^{-3}$'


def animation_save_name(output_direc, rho, flowdata='CMEMS_MEDITERRANEAN', startyear=2010):
    return output_direc + 'SizeTransport_{}_rho_{}_y_{}.mov'.format(flowdata, rho, startyear)