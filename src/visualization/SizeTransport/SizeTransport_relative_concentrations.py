import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string


def SizeTransport_relative_concentrations(scenario, figure_direc, size_list, rho_list, selection, figsize=(20, 10),
                                          fontsize=14):
    # Setting the folder within which we have the output, and where we have the saved timeslices
    output_direc = figure_direc + 'concentrations/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format('SizeTransport')

    # Loading in the data
    prefix = 'horizontal_concentration'
    concentration_dict = {}
    selection_dict = {'average': "overall_concentration", 0: 0, 1: 1, 2: 2}
    for index, size in enumerate(size_list):
        data_dict = vUtils.SizeTransport_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                                   size=size, rho=rho_list[index])
        concentration_dict[size] = data_dict[selection_dict[selection]]
    lon, lat = data_dict['lon'], data_dict['lat']
    Lon, Lat = np.meshgrid(lon, lat)

    # Normalizing the concentration by the lowest non-zero concentration over all the sizes
    normalization_factor = 1e10
    for size in concentration_dict.keys():
        concentration = concentration_dict[size]
        min_non_zero = np.nanmin(concentration[concentration > 0])
        if min_non_zero < normalization_factor:
            normalization_factor = min_non_zero
    for size in concentration_dict.keys():
        concentration_dict[size] /= normalization_factor
        concentration_dict[size][concentration_dict[size] == 0] = np.nan

    # Getting the size of the domain that we want to plot for
    advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                        advection_scenario='CMEMS_MEDITERRANEAN',
                                                        repeat_dt=None)
    adv_file_dict = advection_scenario.file_names

    spatial_domain = np.nanmin(adv_file_dict['LON']), np.nanmax(adv_file_dict['LON']), \
                     np.nanmin(adv_file_dict['LAT']), np.nanmax(adv_file_dict['LAT'])

    # Creating the base figure
    gridspec_shape = (2, 3)
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=gridspec_shape[0], ncols=gridspec_shape[1] + 1, width_ratios=[1, 1, 1, 0.1])

    ax_list = []
    for rows in range(gridspec_shape[0]):
        for columns in range(gridspec_shape[1]):
            ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                       domain=spatial_domain,
                                                       lat_grid_step=5, lon_grid_step=10, resolution='10m'))

    # Setting the colormap, that we will use for coloring the scatter plot according to the particle depth. Then, adding
    # a colorbar.
    norm = colors.LogNorm(vmin=1, vmax=1e5)
    cmap_name = 'inferno_r'
    cmap = plt.cm.ScalarMappable(cmap=cmap_name, norm=norm)
    cax = fig.add_subplot(gs[:, -1])
    cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend='max')
    cbar.set_label(r"Relative Concentration ($C/C_{min}$)", fontsize=fontsize)
    cbar.ax.tick_params(which='major', labelsize=fontsize - 2, length=14, width=2)
    cbar.ax.tick_params(which='minor', labelsize=fontsize - 2, length=7, width=2)

    # Defining the particle sizes and densities that we want to plot, and adding subfigure titles to the corresponding
    # subfigures
    for index, ax in enumerate(ax_list):
        ax.set_title(subfigure_title(index, size_list[index], rho_list[index]), weight='bold', fontsize=fontsize)

    # The actual plotting of the figures
    for index, size in enumerate(size_list):
        ax_list[index].pcolormesh(Lon, Lat, concentration[size], norm=norm, cmap=cmap_name)

    # Saving the figure
    file_name = animation_save_name(output_direc, np.nanmean(rho_list), selection)
    plt.savefig(file_name, bbox_inches='tight')


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


def animation_save_name(output_direc, rho, selection, flowdata='CMEMS_MEDITERRANEAN', startyear=2010):
    selection_dict = {'average': 'TotalAverage', 0: 'year_0', 1: 'year_1', 2: 'year_2'}
    return output_direc + 'SizeTransport_{}_{}_rho_{}_y_{}.mov'.format(flowdata, selection_dict[selection],
                                                                       rho, startyear)