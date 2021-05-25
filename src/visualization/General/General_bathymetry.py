import settings
import utils
from netCDF4 import Dataset
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def General_bathymetry(scenario, figure_direc, figsize=(10, 8), fontsize=14):
    """
    This takes the bathymetry file from an advection scenario and plots it
    :param scenario:
    :param figure_direc:
    :return:
    """
    # Setting the folder within which we have the output
    output_direc = figure_direc + 'General/'
    utils.check_direc_exist(output_direc)

    # Getting the bathymetry data
    file_dict = scenario.file_dict
    dataset = Dataset(file_dict['BATH_filenames'])
    bath_dict = {}
    # Here, DEPTH is actually the bathymetry
    for variable in ['LON', 'LAT', 'DEPTH']:
        bath_dict[variable] = dataset.variables[file_dict['BATH_variables'][variable]][:]

    # Setting the zero depths to nan values
    bath_dict['DEPTH'][bath_dict['DEPTH'] == 0] = np.nan

    # Getting the spatial domain for the figure
    spatial_domain = np.nanmin(bath_dict['LON']), np.nanmax(bath_dict['LON']), \
                     np.nanmin(bath_dict['LAT']), np.nanmax(bath_dict['LAT'])

    # Creating the base figure
    gridspec_shape = (1, 1)
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=gridspec_shape[0], ncols=gridspec_shape[1] + 1, width_ratios=[1, 0.1])

    ax_list = []
    for rows in range(gridspec_shape[0]):
        for columns in range(gridspec_shape[1]):
            ax_list.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                       domain=spatial_domain,
                                                       lat_grid_step=5, lon_grid_step=10, resolution='10m'))

    # Setting the normalization of the colormap
    normalization = colors.LogNorm(vmin=1e0, vmax=1e3)

    # The actual plotting
    Lon, Lat = np.meshgrid(bath_dict['LON'], bath_dict['LAT'])
    cmap = 'binary'
    depth_plot = ax_list[0].meshgrid(Lon, Lat, c=bath_dict['DEPTH'], norm=normalization, cmap='binary')

    cax = fig.add_subplot(gs[:, -1])
    cbar = plt.colorbar(depth_plot, cax=cax, orientation='vertical', extend='max')
    cbar.set_label('Depth (m)', fontsize=fontsize)
    cbar.ax.tick_params(which='major', labelsize=fontsize - 2, length=14, width=2)
    cbar.ax.tick_params(which='minor', labelsize=fontsize - 2, length=7, width=2)

    file_name = output_direc + 'Bathymetry.png'
    plt.savefig(file_name, bbox_inches='tight')
