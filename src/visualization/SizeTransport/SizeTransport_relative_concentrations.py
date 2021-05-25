import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string


def SizeTransport_relative_concentrations(scenario, figure_direc, size_list, rho_list, beach_state, time_selection,
                                          difference=False, reference_size=None, figsize=(20, 10), fontsize=14):
    """

    :param scenario: loading the scenario class to get the grid data
    :param figure_direc: which directory to save the figure in
    :param size_list: the particle sizes that we are plotting
    :param rho_list: the particle densities we are plotting
    :param beach_state: which particles do we want to consider for figures: 'beach', 'afloat' or 'seabed'
    :param time_selection: if 'average', then it uses the average concentration calculated over the entire simulation,
                           otherwise it takes the year of the simulation given as input
    :param difference: if false, we plot just the normalized concentrations, otherwise we plot the difference relative
                       to a given reference_size particle
    :param reference_size: if not specified, we just take the smallest particle size in size_list for the difference
                           calculation
    :param figsize: size of the figure output
    :param fontsize: font size of the figure text/labels
    :return:
    """
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
    for index, size in enumerate(size_list):
        data_dict = vUtils.SizeTransport_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                                   size=size, rho=rho_list[index])
        concentration_dict[size] = data_dict[key_concentration][beach_state]
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

    # Calculating differences relative to the reference particle size
    if difference:
        if reference_size is None:
            reference_size = np.nanmin(size_list)
        for keys in concentration_dict.keys():
            concentration_dict[keys] -= concentration_dict[reference_size]

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
    norm = set_normalization(beach_state, difference)
    cmap_name = {False: 'inferno_r', True: 'bwr'}[difference]
    if not difference:
        cbar_label, extend = r"Relative Concentration ($C/C_{min}$)", 'max'
    elif difference:
        cbar_label, extend = r"Concentration Difference relative to {} mm".format(reference_size * 1e3), 'both'
    cmap = plt.cm.ScalarMappable(cmap=cmap_name, norm=norm)
    cax = fig.add_subplot(gs[:, -1])
    cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend=extend)
    cbar.set_label(cbar_label, fontsize=fontsize)
    cbar.ax.tick_params(which='major', labelsize=fontsize - 2, length=14, width=2)
    cbar.ax.tick_params(which='minor', labelsize=fontsize - 2, length=7, width=2)

    # Defining the particle sizes and densities that we want to plot, and adding subfigure titles to the corresponding
    # subfigures
    for index, ax in enumerate(ax_list):
        ax.set_title(subfigure_title(index, size_list[index], rho_list[index]), weight='bold', fontsize=fontsize)

    # The actual plotting of the figures
    for index, size in enumerate(size_list):
        if beach_state in ['afloat']:
            ax_list[index].pcolormesh(Lon, Lat, concentration_dict[size], norm=norm, cmap=cmap_name)
        else:
            ax_list[index].scatter(Lon.flatten(), Lat.flatten(), c=concentration_dict[size], norm=norm, cmap=cmap_name)

    # Saving the figure
    file_name = plot_save_name(output_direc, np.nanmean(rho_list), time_selection, difference, beach_state,
                               reference_size)
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


def plot_save_name(output_direc, rho, time_selection, difference, beach_state, reference_size,
                   flowdata='CMEMS_MEDITERRANEAN', startyear=2010):
    selection_dict = {'average': 'TotalAverage', 0: 'year_0', 1: 'year_1', 2: 'year_2'}
    difference_dict = {True: 'Difference', False: 'absolute'}
    if not difference:
        return output_direc + 'SizeTransport_{}_{}_{}_{}_rho_{}_y_{}.png'.format(flowdata, beach_state,
                                                                                  difference_dict[difference],
                                                                                  selection_dict[time_selection], rho,
                                                                                  startyear)
    else:
        return output_direc + 'SizeTransport_{}_{}_{}_refsize_{}_{}_rho_{}_y_{}.png'.format(flowdata, beach_state,
                                                                                 difference_dict[difference],
                                                                                            reference_size,
                                                                                 selection_dict[time_selection], rho,
                                                                                 startyear)


def set_normalization(beach_state, difference):
    """
    Setting the normalization that we use for the colormap
    :param beach_state: afloat, beach or seabed
    :param difference: is it absolute concentrations or the difference relative to reference size
    :return:
    """
    if not difference:
        if beach_state == 'afloat':
            vmin, vmax = 1, 1e4
        elif beach_state == 'beach':
            vmin, vmax = 1, 1e5
        elif beach_state == 'seabed':
            vmin, vmax = 1, 1e5
        return colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        linthresh,linscale, vmin, vmax = 1e2, 1, -1e4, 1e4
        return colors.SymLogNorm(linthresh=linthresh, linscale=linscale, vmin=vmin, vmax=vmax)
