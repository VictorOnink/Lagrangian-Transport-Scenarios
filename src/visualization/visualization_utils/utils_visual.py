import settings
import utils
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cpf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import numpy as np


def SizeTransport_load_data(scenario, prefix, data_direc, size,  rho, advection_data='CMEMS_MEDITERRANEAN',
                            shore_time=20, start_year=2010, input='Lebreton'):
    """
    Loading the data we want for SizeTransport analysis output, which will generally just differ in terms of which
    particle size and density
    :param scenario:
    :param prefix:
    :param data_direc:
    :param size:
    :param rho:
    :param advection_data:
    :param shore_time:
    :param start_year:
    :param input:
    :return:
    """
    file_name = scenario._file_names(new=True, advection_data=advection_data, shore_time=shore_time, init_size=size,
                                     init_density=rho, start_year=start_year, input=input, run=settings.RUN,
                                     restart=settings.RESTART)
    return data_direc + utils._analysis_save_file_name(input_file=file_name, prefix=prefix)


def cartopy_standard_map(fig, gridspec, row, column, domain, resolution='50m', add_gridlines=True, add_gridlabels=True,
                         lat_grid_step=20, lon_grid_step=30, label_size=14):
    """
    A nice basic function that can be used to create standardized maps
    :param axis:
    :return:
    """
    axis = fig.add_subplot(gridspec[row, column], projection=ccrs.PlateCarree())
    # Setting the domain of the map
    lon_min, lon_max, lat_min, lat_max = domain
    axis.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Adding coastlines, borders, land and ocean shapefiles
    axis.coastlines(resolution=resolution)
    axis.add_feature(cpf.BORDERS.with_scale(resolution), edgecolor='black', zorder=100)
    axis.add_feature(cpf.LAND.with_scale(resolution), facecolor='gray')
    axis.add_feature(cpf.OCEAN.with_scale(resolution), facecolor='white')

    # Adding gridlines and axis labels
    if add_gridlines:
        grid = axis.gridlines(crs=ccrs.PlateCarree(),  # specify the projection being used
                              draw_labels=add_gridlabels,  # Add labels
                              linestyle='-',  # style
                              color='black'
                              )
        # Here we can choose along which axes we want to have the labels included
        grid.top_labels = False
        grid.bottom_labels = True
        grid.left_labels = True
        grid.right_labels = False
        # Formatting of the labels, since that they include N/S or E/W
        grid.xformatter = LONGITUDE_FORMATTER
        grid.yformatter = LATITUDE_FORMATTER

        # Determine where we want labels
        grid.xlocator = mticker.FixedLocator(np.arange(int(lon_min), int(lon_max) + lon_grid_step, lon_grid_step))
        grid.ylocator = mticker.FixedLocator(np.arange(int(lat_min), int(lat_max) + lat_grid_step, lat_grid_step))
        # Here we can change the appearances of the labels
        grid.xlabel_style = {'size': label_size, 'color': 'black', 'weight': 'normal'}
        grid.ylabel_style = {'size': label_size, 'color': 'black', 'weight': 'normal'}

    # Setting the axis
    axis.set_aspect('auto', adjustable=None)
    return axis

