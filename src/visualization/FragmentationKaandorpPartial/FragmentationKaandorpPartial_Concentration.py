import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string


def FragmentationKaandorpPartial_Concentration(scenario, figure_direc, rho, shore_time, beach_state, simulation_year,
                                               lambda_frag, figsize=(20, 10), ax_ticklabel_size=12, ax_label_size=14):
    # Setting the folder within which we have the output, and where we have the saved timeslices
    output_direc = figure_direc + 'concentrations/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format('SizeTransport')

    # Loading in the data
    prefix = 'horizontal_concentration'
    if simulation_year == 'average':
        key_concentration = "overall_concentration"
    else:
        key_concentration = utils.analysis_simulation_year_key(simulation_year)
    data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=scenario, prefix=prefix,
                                                              data_direc=data_direc, shore_time=shore_time,
                                                              lambda_frag=lambda_frag, rho=rho, postprocess=True)
    concentration_dict = data_dict[key_concentration][beach_state]
    Lon, Lat = np.meshgrid(data_dict['lon'], data_dict['lat'])

    # Normalizing the concentration by the lowest non-zero concentration over all the sizes
    normalization_factor = 1e10
    for size in concentration_dict.keys():
        concentration = concentration_dict[size]
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
    norm = set_normalization(beach_state)
    cmap_name = 'inferno_r'
    cbar_label, extend = r"Relative Concentration ($C/C_{min}$)", 'max'
    cmap = plt.cm.ScalarMappable(cmap=cmap_name, norm=norm)
    cax = fig.add_subplot(gs[:, -1])
    cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend=extend)
    cbar.set_label(cbar_label, fontsize=ax_label_size)
    cbar.ax.tick_params(which='major', labelsize=ax_ticklabel_size, length=14, width=2)
    cbar.ax.tick_params(which='minor', labelsize=ax_ticklabel_size, length=7, width=2)

    # Defining the particle sizes and densities that we want to plot, and adding subfigure titles to the corresponding
    # subfigures
    for index, ax in enumerate(ax_list):
        ax.set_title(subfigure_title(index), weight='bold', fontsize=ax_label_size)

    # The actual plotting of the figures
    for index, ax in enumerate(ax_list):
        if beach_state in ['adrift']:
            ax.pcolormesh(Lon, Lat, concentration_dict[index], norm=norm, cmap=cmap_name, zorder=200)
        else:
            ax.scatter(Lon.flatten(), Lat.flatten(), c=concentration_dict[index], norm=norm, cmap=cmap_name, zorder=200)

    file_name = output_direc + 'Concentrations_{}_year={}_lambda_f={}_st={}_rho={}.png'.format(beach_state,
                                                                                               simulation_year,
                                                                                               lambda_frag, shore_time,
                                                                                               rho)
    plt.savefig(file_name, bbox_inches='tight')


def subfigure_title(index):
    """
    setting the title of the subfigure
    :param index:
    :param size:
    :param rho:
    :return:
    """
    particle_size = utils.size_range(size_class_number=index, units='mm')
    title = '({}) r = {:.3f} mm'.format(string.ascii_lowercase[index], particle_size)
    return title


def set_normalization(beach_state):
    """
    Setting the normalization that we use for the colormap
    :param beach_state: adrift, beach or seabed
    :param difference: is it absolute concentrations or the difference relative to reference size
    :return:
    """
    if beach_state == 'adrift':
        vmin, vmax = 1, 1e4
    elif beach_state == 'beach':
        vmin, vmax = 1, 1e5
    return colors.LogNorm(vmin=vmin, vmax=vmax)


