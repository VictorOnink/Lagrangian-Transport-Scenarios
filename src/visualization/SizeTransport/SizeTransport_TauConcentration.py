import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string
import cmocean.cm as cmo


def SizeTransport_TauConcentration(scenario, figure_direc, size_selection, tau_list, beach_state, rho=920,
                                   time_selection=0, figsize=(20, 10), fontsize=14):
    # Setting the folder within which we have the output, and where we have the saved timeslices
    output_direc = figure_direc + 'concentrations/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format('SizeTransport')

    # Loading in the data
    prefix = 'horizontal_concentration'
    concentration_dict = {}
    if time_selection == 'average':
        key_concentration = "overall_concentration"
    else:
        key_concentration = utils.analysis_simulation_year_key(time_selection)
    for index, tau in enumerate(tau_list):
        data_dict = vUtils.SizeTransport_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                                   size=size_selection, rho=rho, tau=tau)
        concentration_dict[index] = data_dict[key_concentration][beach_state]
    lon, lat = data_dict['lon'], data_dict['lat']
    Lon, Lat = np.meshgrid(lon, lat)

    # Normalizing the concentration by the lowest non-zero concentration over all the sizes
    normalization_factor = 1e10
    for size in concentration_dict.keys():
        concentration = concentration_dict[size]
        if np.sum(concentration > 0) > 0:
            min_non_zero = np.nanmin(concentration[concentration > 0])
            if min_non_zero < normalization_factor:
                normalization_factor = min_non_zero
    for size in concentration_dict.keys():
        concentration_dict[size] /= normalization_factor

    # Setting zero values to nan
    for size in concentration_dict.keys():
        concentration_dict[size][concentration_dict[size] == 0] = np.nan

    # Getting the size of the domain that we want to plot for
    advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                        advection_scenario='CMEMS_MEDITERRANEAN',
                                                        repeat_dt=None)
    adv_file_dict = advection_scenario.file_names

    spatial_domain = np.nanmin(adv_file_dict['LON']), np.nanmax(adv_file_dict['LON']), \
                     np.nanmin(adv_file_dict['LAT']), np.nanmax(adv_file_dict['LAT'])

    # Creating the base figure
    gridspec_shape = (2, 2)
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=gridspec_shape[0], ncols=gridspec_shape[1] + 1, width_ratios=[1, 1, 0.1])

    ax_list = []
    for rows in range(gridspec_shape[0]):
        for columns in range(gridspec_shape[1]):
            ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                       domain=spatial_domain,
                                                       lat_grid_step=5, lon_grid_step=10, resolution='10m'))

    # Setting the colormap, that we will use for coloring the scatter plot according to the particle depth. Then, adding
    # a colorbar.
    norm = set_normalization(beach_state)
    cmap_name = cmo.solar
    cbar_label, extend = r"Relative Concentration ($C/C_{min}$)", 'max'
    cmap = plt.cm.ScalarMappable(cmap=cmap_name, norm=norm)
    cax = fig.add_subplot(gs[:, -1])
    cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend=extend)
    cbar.set_label(cbar_label, fontsize=fontsize)
    cbar.ax.tick_params(which='major', labelsize=fontsize - 2, length=14, width=2)
    cbar.ax.tick_params(which='minor', labelsize=fontsize - 2, length=7, width=2)

    # adding subfigure titles to the corresponding subfigures
    for index, ax in enumerate(ax_list):
        ax.set_title(subfigure_title(index, size_selection, tau_list), weight='bold', fontsize=fontsize)

    # The actual plotting of the figures
    for index, tau in enumerate(tau_list):
        if beach_state in ['adrift']:
            ax_list[index].pcolormesh(Lon, Lat, concentration_dict[index], norm=norm, cmap=cmap_name, zorder=200)
        else:
            ax_list[index].scatter(Lon.flatten(), Lat.flatten(), c=concentration_dict[index], norm=norm, cmap=cmap_name,
                                   zorder=200)
    # Saving the figure
    file_name = plot_save_name(output_direc, rho, time_selection, beach_state, size_selection)
    plt.savefig(file_name, bbox_inches='tight')


def plot_save_name(output_direc, rho, time_selection, beach_state, size,
                   flowdata='CMEMS_MEDITERRANEAN', startyear=2010):
    selection_dict = {'average': 'TotalAverage', 0: 'year_0', 1: 'year_1', 2: 'year_2'}
    return output_direc + 'SizeTransport_TAUCOMP_{}_{}_{:.3f}_{}_rho_{}_y_{}.png'.format(flowdata, beach_state, size * 1e3,
                                                                                     selection_dict[time_selection],
                                                                                     rho, startyear)


def set_normalization(beach_state):
    """
    Setting the normalization that we use for the colormap
    :param beach_state: adrift, beach or seabed
    :return:
    """
    if beach_state == 'adrift':
        vmin, vmax = 1, 1e4
    elif beach_state == 'beach':
        vmin, vmax = 1, 1e5
    elif beach_state == 'seabed':
        vmin, vmax = 1, 1e6
    return colors.LogNorm(vmin=vmin, vmax=vmax)


def subfigure_title(index, size, tau_list):
    """
    setting the title of the subfigure
    :param index:
    :param size:
    :param rho:
    :return:
    """
    title = '({}) r = {:.3f} mm'.format(string.ascii_lowercase[index], size * 1e3)
    title += r', $\tau_{crit} = $' + '{} Pa'.format(tau_list[index])
    return title
