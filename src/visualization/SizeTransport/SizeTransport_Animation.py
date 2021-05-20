import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from advection_scenarios import advection_files
import numpy as np
from datetime import datetime, timedelta
import string


def SizeTransport_Animation(figure_direc, figsize=(20, 10), fontsize=14):
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

    # Setting the folder within which we have the output
    output_direc = figure_direc + 'animations/'

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
    rho_list = np.ones(size_list.shape) * 920
    for index, ax in enumerate(ax_list):
        ax.set_title(subfigure_title(index, size_list[index], rho_list[index]), weight='bold', fontsize=fontsize)

    # Setting the time range for which we want to create the simulation
    start_time = datetime(2010, 1, 1, 0)
    end_time = datetime(2013, 1, 1, 0)
    time_step = timedelta(hours=12)

    # Setting the output name of the animation, and saving the output
    output_name = output_direc + 'test.jpg'
    plt.savefig(output_name, bbox_inches='tight')


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