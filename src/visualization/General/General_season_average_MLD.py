import settings
import utils
from netCDF4 import Dataset
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean.cm as cmo
import cartopy.crs as ccrs


def General_season_average_MLD(scenario, figure_direc, figsize=(10, 8), fontsize=14):
    # Setting the folder within which we have the output
    output_direc = figure_direc + 'General/'
    utils.check_direc_exist(output_direc)

    # Getting the grid data
    file_dict = scenario.file_dict
    LON_GRID, LAT_GRID = file_dict['LON'], file_dict['LAT']
    spatial_domain = np.nanmin(LON_GRID), np.nanmax(LON_GRID), \
                     np.nanmin(LAT_GRID), np.nanmax(LAT_GRID)

    # Loading the MLD mean fields
    season_list = ['DJF', 'MAM', 'JJA', 'SON']
    MLD_dict = {}
    for season in season_list:
        file_name = settings.DATA_DIR_SERVERS[settings.server] + '/CMEMS_MED/mean_{}_2010_2012_AXML.nc'.format(season)
        dataset = Dataset(file_name)
        MLD_dict[season] = np.nanmean(dataset.variables['mlotst'][:], axis=0)
    MLD_dict['latlon'] = np.meshgrid(dataset.variables['lat'][:], dataset.variables['lon'][:])

    # Creating the figure
    # Creating the base figure
    gridspec_shape = (2, 2)
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=gridspec_shape[0], ncols=gridspec_shape[1] + 1, width_ratios=[1, 1, 0.2])

    ax_list = []
    for rows in range(gridspec_shape[0]):
        for columns in range(gridspec_shape[1]):
            ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                       domain=spatial_domain,
                                                       lat_grid_step=5, lon_grid_step=10, resolution='10m',
                                                       land_zorder=90))

    # Setting the colormap, that we will use for coloring the scatter plot according to the particle depth. Then, adding
    # a colorbar.
    norm = colors.Normalize(vmin=0.0, vmax=100.0)
    cmap = plt.cm.ScalarMappable(cmap=cmo.haline_r, norm=norm)
    cax = fig.add_subplot(gs[:, -1])
    cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend='max')
    cbar.set_label(r"Depth (m)", fontsize=fontsize)
    cbar.ax.tick_params(which='major', labelsize=fontsize - 2, length=14, width=2)
    cbar.ax.tick_params(which='minor', labelsize=fontsize - 2, length=7, width=2)

    file_name = output_direc + 'Seasonal_average_MLD_{}_2010-12.png'.format(settings.ADVECTION_DATA)
    plt.savefig(file_name, bbox_inches='tight')
