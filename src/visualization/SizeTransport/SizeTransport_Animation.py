import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
from advection_scenarios import advection_files
import numpy as np
from datetime import datetime, timedelta
import string
import matplotlib.animation as animation


def SizeTransport_Animation(scenario, figure_direc, figsize=(20, 10), fontsize=14):
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
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeslices/{}/'.format('SizeTransport')

    # Creating the base figure
    gridspec_shape = (2, 3)
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=gridspec_shape[0], ncols=gridspec_shape[1])

    ax_list = []
    for rows in range(gridspec_shape[0]):
        for columns in range(gridspec_shape[1]):
            ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                       domain=spatial_domain,
                                                       lat_grid_step=5, lon_grid_step=10, resolution='10m'))
    # Defining the particle sizes and densities that we want to plot, and adding subfigure titles to the corresponding
    # subfigures
    size_list = np.array([500, 100, 50, 10, 5, 1]) * 1e-5
    rho_list = np.ones(size_list.shape, dtype=int) * 920
    for index, ax in enumerate(ax_list):
        ax.set_title(subfigure_title(index, size_list[index], rho_list[index]), weight='bold', fontsize=fontsize)

    # Setting the time range for which we want to create the simulation
    current_time = datetime(2010, 1, 1, 0)
    end_time = datetime(2013, 1, 1, 0)
    time_step = timedelta(hours=12)
    time_list = []
    while current_time < end_time:
        time_list.append(current_time)
        current_time += time_step
    frame_number = len(time_list)

    # Now, the actual animation part
    # Setting the initial values of the x and y, which will later be filled by lon and lat
    plot_list = []
    for ax in ax_list:
        plot_list.append(ax.scatter(0, 0, size=4, alpha=1, zorder=1000)[0])

    # Initializing the plots on each axis
    def init():
        for plot in plot_list:
            plot.set_data([], [])
        return plot_list

    def animate(frame_index):
        for index, size in enumerate(size_list):
            # Loading the dictionary with the data
            prefix = 'timeslices_{}'.format(time_list[frame_index].strftime("%Y-%m-%d-%H-%M-%S"))
            data_dict = vUtils.SizeTransport_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                                       size=size, rho=rho_list[index])
            lon, lat = data_dict['lon'], data_dict['lat']
            # Updating the plot on each axis with the data
            plot_list[index].set_data(lon, lat)
        return plot_list

    # Calling the animator
    animator = animation.FuncAnimation(plt.gcf(), animate, init_func=init,
                                       frames=2, interval=100, blit=True)

    # Saving the animation
    animator.save(filename=output_direc + 'blah.mov', fps=2, extra_args=['-vcodec', 'libx264'])


def subfigure_title(index, size, rho):
    """
    setting the title of the subfigure
    :param index:
    :param size:
    :param rho:
    :return:
    """
    alphabet = string.ascii_lowercase
    return '({}) r = {} mm, '.format(alphabet[index], size * 1e3) + r'$\rho$ = ' + '{} kg m'.format(rho) + r'$^{-3}$'