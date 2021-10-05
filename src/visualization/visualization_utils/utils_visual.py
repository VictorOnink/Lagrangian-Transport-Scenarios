import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cpf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.dates as mdates
import utils


def SizeTransport_linestyle_SEABED_CRIT(tau):
    """
    Returning the linestyle for line-based figures depending on the
    :param tau:
    :return:
    """
    return {100.0: 'dotted', 0.14: 'dashed', 0.025: 'dashdot', 0.0: 'solid'}[float(tau)]


def cartopy_standard_map(fig, gridspec, row, column, domain, resolution='50m', add_gridlines=True, add_gridlabels=True,
                         lat_grid_step=20, lon_grid_step=30, label_size=14, land_zorder=1, ocean_zorder=1,
                         line_zorder=101, border_color='black', land_color='lightgray', ocean_color='white'):
    """
    A nice basic function that can be used to create standardized maps
    :param fig:
    :param gridspec:
    :param row:
    :param column:
    :param domain:
    :param resolution:
    :param add_gridlines:
    :param add_gridlabels:
    :param lat_grid_step:
    :param lon_grid_step:
    :param label_size:
    :param land_zorder:
    :param ocean_zorder:
    :param line_zorder:
    :return:
    """
    axis = fig.add_subplot(gridspec[row, column], projection=ccrs.PlateCarree())
    # Setting the domain of the map
    lon_min, lon_max, lat_min, lat_max = domain
    axis.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Adding coastlines, borders, land and ocean shapefiles
    axis.coastlines(resolution=resolution)
    axis.add_feature(cpf.BORDERS.with_scale(resolution), edgecolor=border_color, zorder=line_zorder)
    axis.add_feature(cpf.LAND.with_scale(resolution), facecolor=land_color, zorder=land_zorder)
    axis.add_feature(cpf.OCEAN.with_scale(resolution), facecolor=ocean_color, zorder=ocean_zorder)

    # Adding gridlines and axis labels
    if add_gridlines:
        grid = axis.gridlines(crs=ccrs.PlateCarree(),  # specify the projection being used
                              draw_labels=add_gridlabels,  # Add labels
                              linestyle='-',  # style
                              color='black',
                              zorder=line_zorder
                              )
        # Here we can choose along which axes we want to have the labels included
        grid.top_labels = False
        grid.right_labels = False
        if column == 0:
            grid.left_labels = True
        else:
            grid.left_labels = False
        if row == (gridspec.nrows - 1):
            grid.bottom_labels = True
        else:
            grid.bottom_labels = False
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


def discrete_color_from_cmap(index, subdivisions, cmap='viridis_r'):
    cmap_steps = plt.cm.get_cmap(cmap, subdivisions)
    return cmap_steps(index)


def base_figure(fig_size, ax_range, y_label, x_label, ax_label_size, ax_ticklabel_size,
                shape=(1, 1), plot_num=1, all_x_labels=False, legend_axis=False, log_yscale=False, log_xscale=False,
                x_time_axis=False, width_ratios=None, height_ratios=None, all_y_labels=True, add_twinx=False,
                twinx_ax_range=None, twinx_y_label=None, twin_x_all_columns=False, log_twinxscale=False):
    """
    Function creating the base figure that we use as a foundation for almost all figures
    :param log_twinxscale:
    :param twinx_y_label:
    :param twinx_ax_range:
    :param twin_x_all_columns:
    :param add_twinx:
    :param log_yscale: if True, the y axis has a log scale
    :param log_xscale: if True, the x axis has a log scale
    :param x_time_axis: if True, the x axis is a time axis
    :param height_ratios: ratios of the height of the subfigures (must be a list with the same size as the number of
                          rows
    :param width_ratios: ratios of the width of the subfigures (must be a list with the same size as the number of
                          columns
    :param fig_size: size of the figure
    :param ax_range: the limits of the x and y axes
    :param y_label: the y label
    :param x_label: the x label
    :param ax_label_size: the fontsize of the axis labels
    :param shape: the shape of the array (rows, columns)
    :param plot_num: how many subplots we want to create (e.g. in case we have a 2x3 figure but only want 5 panels)
    :param all_x_labels: if True, all subplots in the bottom row of teh figure will have x labels, otherwise just the
                         middle one
    :param all_y_labels: if True, all subplots will have the y axis labeled, otherwise just the middle row
    :param legend_axis: if true, we add an additional column in which we can add the legend (in case it is too big to
                        fit within a subplot)
    :return:
    """
    # Loading the axis limits
    xmax, xmin, ymax, ymin = ax_range
    if log_xscale:
        if xmax < 0 or xmin < 0:
            xlog_type = 'symlog'
        else:
            xlog_type = 'log'
    if log_yscale:
        if ymax < 0 or ymin < 0:
            ylog_type = 'symlog'
        else:
            ylog_type = 'log'
    # Checking if we have the necessary settings for twin_x
    if add_twinx:
        assert twinx_ax_range is not None, 'Please specify axis range for twin axis'
        xmax_twin, xmin_twin, ymax_twin, ymin_twin = twinx_ax_range
        if log_twinxscale:
            if xmax_twin < 0 or xmin_twin < 0:
                twinxlog_type = 'symlog'
            else:
                twinxlog_type = 'log'
    # Creating the figure
    fig = plt.figure(figsize=fig_size)
    if legend_axis:
        grid = GridSpec(nrows=shape[0], ncols=shape[1] + 1, figure=fig, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    else:
        grid = GridSpec(nrows=shape[0], ncols=shape[1], figure=fig, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    ax, fig_index = [], 1
    for row in range(shape[0]):
        for column in range(shape[1]):
            ax_sub = fig.add_subplot(grid[row, column])
            # Setting the axis limits
            ax_sub.set_ylim((ymin, ymax))
            ax_sub.set_xlim((xmin, xmax))
            # Setting the tick parameters and setting the axis scale
            ax_sub.tick_params(axis='both', labelsize=ax_ticklabel_size)
            if log_xscale:
                ax_sub.set_xscale(xlog_type)
            if log_yscale:
                ax_sub.set_yscale(ylog_type)
            if x_time_axis:
                years = mdates.YearLocator()
                months = mdates.MonthLocator()
                yearsFmt = mdates.DateFormatter('%Y')
                ax_sub.xaxis.set_major_locator(years)
                ax_sub.xaxis.set_minor_locator(months)
                ax_sub.xaxis.set_major_formatter(yearsFmt)
            # Creating the twinx axis
            if add_twinx:
                if not twin_x_all_columns and column == shape[1] - 1:
                    twin_ax_sub = ax_sub.twinx()
                    twin_ax_sub.set_ylim((ymin_twin, ymax_twin))
                    twin_ax_sub.set_xlim((xmin_twin, xmax_twin))
                else:
                    twin_ax_sub = ax_sub.twinx()
                    twin_ax_sub.set_ylim((ymin_twin, ymax_twin))
                    twin_ax_sub.set_xlim((xmin_twin, xmax_twin))
                if log_twinxscale:
                    twin_ax_sub.set_yscale(twinxlog_type)
            # Labeling the x and y axes
            # Only add y labels if we are in the first column
            if column == 0:
                if all_y_labels or row == shape[0] // 2:
                    ax_sub.set_ylabel(y_label, fontsize=ax_label_size)
            else:
                ax_sub.tick_params(labelleft=False)
            if add_twinx:
                if column == shape[1] - 1 and all_y_labels or row == shape[0] // 2:
                    twin_ax_sub.set_ylabel(twinx_y_label, fontsize=ax_label_size)
            # Only add x labels if we are in the bottom row, and only to the middle one unless all_x_labels == True
            if row == (shape[0] - 1):
                if not all_x_labels and column % 2 is 1:
                    ax_sub.set_xlabel(x_label, fontsize=ax_label_size)
                elif all_x_labels:
                    ax_sub.set_xlabel(x_label, fontsize=ax_label_size)
                else:
                    ax_sub.tick_params(labelbottom=False)
            else:
                ax_sub.tick_params(labelbottom=False)
            # Adding the axis to the list and continuiing to the next plot
            ax.append(ax_sub)
            fig_index += 1
            if fig_index > plot_num:
                break
    if legend_axis:
        # Adding a legend axis
        ax_legend = fig.add_subplot(grid[:, -1])
        ax_legend.set_axis_off()
        ax.append(ax_legend)
    return tuple(ax)

