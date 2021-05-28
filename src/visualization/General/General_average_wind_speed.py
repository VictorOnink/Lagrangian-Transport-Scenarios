import settings
import utils
from netCDF4 import Dataset
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean
import cartopy.crs as ccrs


def General_average_wind_speed(scenario, figure_direc, figsize=(10, 8), fontsize=14):
    # Setting the folder within which we have the output
    output_direc = figure_direc + 'General/'
    utils.check_direc_exist(output_direc)

    # Getting the grid data
    file_dict = scenario.file_dict
    LON_GRID, LAT_GRID = file_dict['LON'], file_dict['LAT']
    spatial_domain = np.nanmin(LON_GRID), np.nanmax(LON_GRID), \
                     np.nanmin(LAT_GRID), np.nanmax(LAT_GRID)


    # Loading the mean wind data
    wind_file = settings.DATA_DIR_SERVERS[settings.SERVER] + 'Wind/ERA5-wind10m-average.nc'
    dataset = Dataset(wind_file)
    u10, v10 = dataset.variables['u10'][0, :, :], dataset.variables['v10'][0, :, :]
    lon, lat = dataset.variables['longitude'][:], dataset.variables['latitude'][:]

    # Getting the wind speed
    wind_magnitude = np.sqrt(np.square(u10) + np.square(v10))

    # Creating the base figure
    gridspec_shape = (1, 1)
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=gridspec_shape[0], ncols=gridspec_shape[1] + 1, width_ratios=[1, 0.1])

    ax_list = []
    for rows in range(gridspec_shape[0]):
        for columns in range(gridspec_shape[1]):
            ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                       domain=spatial_domain,
                                                       lat_grid_step=5, lon_grid_step=10, resolution='10m',
                                                       land_zorder=90))

    # Setting the normalization of the colormap
    normalization = colors.LogNorm(vmin=1e-2, vmax=1e1)

    # The actual plotting
    Lon, Lat = np.meshgrid(lon, lat)
    cmap = cmocean.cm.speed
    depth_plot = plt.pcolormesh(Lon, Lat, wind_magnitude, norm=normalization, cmap=cmap, zorder=50,
                                transform=ccrs.PlateCarree())

    cax = fig.add_subplot(gs[:, -1])
    cbar = plt.colorbar(depth_plot, cax=cax, orientation='vertical', extend='both')
    cbar.set_label(r'Mean Wind Speed (m s$^{-1}$)', fontsize=fontsize)
    cbar.ax.tick_params(which='major', labelsize=fontsize - 2, length=14, width=2)
    cbar.ax.tick_params(which='minor', labelsize=fontsize - 2, length=7, width=2)

    file_name = output_direc + 'AverageWindSpeed_2010-2015.png'
    plt.savefig(file_name, bbox_inches='tight')



