import settings
import utils
from netCDF4 import Dataset
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean.cm as cmo
import cartopy.crs as ccrs
import string


def General_season_average_MLD(scenario, figure_direc, figsize=(16, 12), fontsize=16):
    utils.print_statement('Creating a figure of the season average MLD', to_print=True)
    # Setting the folder within which we have the output
    output_direc = figure_direc + 'General/'
    utils.check_direc_exist(output_direc)

    # Getting the grid data
    utils.print_statement('Getting the grid data', to_print=True)
    file_dict = scenario.file_dict
    LON_GRID, LAT_GRID = file_dict['LON'], file_dict['LAT']
    spatial_domain = np.nanmin(LON_GRID), np.nanmax(LON_GRID), \
                     np.nanmin(LAT_GRID), np.nanmax(LAT_GRID)

    # Loading the MLD mean fields
    utils.print_statement('Loading the MLD fields', to_print=True)
    season_list = ['DJF', 'MAM', 'JJA', 'SON']
    # season_list = ['20100101', '20100401', '20100701', '20101001']
    MLD_dict = {}
    for season in season_list:
        file_name = settings.DATA_DIR_SERVERS[settings.SERVER] + 'CMEMS_MED/mean_{}_AMXL.nc'.format(season)
        # file_name = settings.DATA_DIR_SERVERS[settings.SERVER] + 'CMEMS_MED/{}_d-CMCC--AMXL-MFSe3r1-MED-b20200901_re-sv01.00.nc'.format(season)
        dataset = Dataset(file_name)
        MLD_dict[season] = np.nanmean(dataset.variables['mlotst'][:], axis=0)

    # Creating the figure
    utils.print_statement('Setting up the figure', to_print=True)
    gridspec_shape = (2, 2)
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=gridspec_shape[0], ncols=gridspec_shape[1] + 1, width_ratios=[1, 1, 0.1])

    ax_list = []
    for rows in range(gridspec_shape[0]):
        for columns in range(gridspec_shape[1]):
            ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                       domain=spatial_domain, add_gridlines=False, resolution='10m',
                                                       land_zorder=90, ocean_zorder=0))
    # Adding subplot titles
    for ax_index, ax in enumerate(ax_list):
        ax.set_title(subfigure_title(ax_index, season_list[ax_index]), fontsize=fontsize)

    # Setting the colormap, that we will use for coloring the scatter plot according to the particle depth. Then, adding
    # a colorbar.
    # norm = colors.Normalize(vmin=0.0, vmax=100.0)
    norm = colors.LogNorm(vmin=1e0, vmax=1e2)
    cmap_name =cmo.haline_r
    cmap = plt.cm.ScalarMappable(cmap=cmap_name, norm=norm)
    cax = fig.add_subplot(gs[:, -1])
    cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend='max')
    cbar.set_label(r"Depth (m)", fontsize=fontsize)
    cbar.ax.tick_params(which='major', labelsize=fontsize - 2, length=14, width=2)
    cbar.ax.tick_params(which='minor', labelsize=fontsize - 2, length=7, width=2)

    # Plotting the actual mean MLD
    Lon, Lat = np.meshgrid(dataset.variables['lon'][:], dataset.variables['lat'][:])
    for ax_index, ax in enumerate(ax_list):
        ax.pcolormesh(Lon, Lat, MLD_dict[season_list[ax_index]], cmap=cmap_name, zorder=10)

    # Saving the figure
    utils.print_statement('Saving the figure', to_print=True)
    file_name = output_direc + 'Seasonal_average_MLD_{}_2010-12.png'.format(settings.ADVECTION_DATA)
    plt.savefig(file_name, bbox_inches='tight')


def subfigure_title(index, season):
    alphabet = string.ascii_lowercase
    return '({}) {}'.format(alphabet[index], season)
