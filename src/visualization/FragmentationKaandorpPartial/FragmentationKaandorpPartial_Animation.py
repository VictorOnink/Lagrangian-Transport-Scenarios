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


def FragmentationKaandorpPartial_Animation(scenario, figure_direc, shore_time, lambda_frag,
                                           figsize=(20, 10), ax_label_size=18, tick_label_size=16):
    """
    Here we want to make an animation of the
    :return:
    """
    # Getting the size of the domain that we want to plot for
    advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                        advection_scenario='CMEMS_MEDITERRANEAN',
                                                        repeat_dt=None)
    adv_file_dict = advection_scenario.file_names

    spatial_domain = np.nanmin(adv_file_dict['LON']),  np.nanmax(adv_file_dict['LON']), \
                     np.nanmin(adv_file_dict['LAT']), np.nanmax(adv_file_dict['LAT'])

    # Setting the folder within which we have the output, and where we have the saved timeslices
    output_direc = figure_direc + 'animations/'
    utils.check_direc_exist(output_direc)
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeslices/FragmentationKaandorpPartial/'

    # Creating the base figure
    gridspec_shape = (1, 1)
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=gridspec_shape[0], ncols=gridspec_shape[1] + 1, width_ratios=[1, 0.1])

    ax_list = []
    for rows in range(gridspec_shape[0]):
        for columns in range(gridspec_shape[1]):
            ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                       domain=spatial_domain, label_size=tick_label_size,
                                                       lat_grid_step=5, lon_grid_step=10, resolution='10m'))

    # Setting the colormap, that we will use for coloring the scatter plot according to the particle depth. Then, adding
    # a colorbar.
    norm = colors.Normalize(vmin=0.0, vmax=100.0)
    cmap = plt.cm.ScalarMappable(cmap=cmo.haline_r, norm=norm)
    cax = fig.add_subplot(gs[:, -1])
    cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend='max')
    cbar.set_label(r"Depth (m)", fontsize=ax_label_size)
    cbar.ax.tick_params(which='major', labelsize=tick_label_size, length=14, width=2)
    cbar.ax.tick_params(which='minor', labelsize=tick_label_size, length=7, width=2)

    # Setting the time range for which we want to create the simulation
    current_time = datetime(2010, 1, 1, 0)
    end_time = datetime(2010, 12, 31, 12)
    time_step = timedelta(hours=12)
    time_list = []
    while current_time < end_time:
        time_list.append(current_time)
        current_time += time_step
    frame_number = len(time_list)

    # Setting a text box to give the date
    props = dict(boxstyle='round', facecolor='white', alpha=1)
    text = ax_list[0].text(0.02, 0.02, 'initial', horizontalalignment='left', verticalalignment='bottom',
                           transform=ax_list[0].transAxes, bbox=props, fontsize=ax_label_size, zorder=200)

    # Now, the actual animation part
    # Setting the initial values of the x and y, which will later be filled by lon and lat
    plot_list = []
    for ax in ax_list:
        plot_list.append(ax.scatter(0, 0, c=0, s=10, alpha=1, zorder=1000, cmap=cmo.haline_r, norm=norm))

    # Initializing the plots on each axis
    def init():
        for plot in plot_list:
            plot.set_offsets(np.c_[[], []])
        text.set_text('initial 2')
        return plot_list

    def animate(frame_index):
        utils.print_statement("we are at index {} of {}".format(frame_index, frame_number), to_print=True)
        date = time_list[frame_index].strftime("%Y-%m-%d-%H-%M-%S")
        for ax_index, size in enumerate(ax_list):
            # Loading the dictionary with the data
            prefix = 'timeslices_{}'.format(date)
            data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                                                      lambda_frag=lambda_frag, rho=920, shore_time=shore_time)
            lon, lat, depth = data_dict['lon'], data_dict['lat'], data_dict['z'].astype(int)
            # Updating the plot on each axis with the data
            plot_list[ax_index].set_offsets(np.c_[lon, lat])
            plot_list[ax_index].set_array(depth)
        text.set_text(time_list[frame_index].strftime("%Y-%m-%d"))
        return plot_list

    # Calling the animator
    animator = animation.FuncAnimation(plt.gcf(), animate, init_func=init,
                                       frames=frame_number, interval=100, blit=True)

    # Saving the animation
    animator.save(filename=animation_save_name(output_direc=output_direc, shore_time=shore_time,
                                               lambda_frag=lambda_frag),
                  fps=10, extra_args=['-vcodec', 'libx264'])


def animation_save_name(output_direc, shore_time, lambda_frag, flowdata='CMEMS_MEDITERRANEAN', file_type='.mov'):
    return output_direc + 'FragmentationKaandorpPartial_st={}_lam_f={}'.format(flowdata, shore_time, lambda_frag) + \
           file_type