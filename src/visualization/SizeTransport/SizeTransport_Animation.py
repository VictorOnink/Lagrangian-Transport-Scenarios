import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from advection_scenarios import advection_files
import numpy as np


def SizeTransport_Animation(figure_direc, figsize=(14, 10)):
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
    gridspec_shape = (2, 2)
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=gridspec_shape[0], ncols=gridspec_shape[1])

    ax_list = []
    for rows in range(gridspec_shape[0]):
        for columns in range(gridspec_shape[1]):
            ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs[rows, columns], domain=spatial_domain,
                                                       lat_grid_step=5, lon_grid_step=10, resolution='10m'))

    # Setting the output name of the animation, and saving the output
    output_name = output_direc + 'test.jpg'
    plt.savefig(output_name)
