import settings
import utils
from netCDF4 import Dataset
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean
import pandas as pd


def General_input_scenario(scenario, figure_direc, figsize=(10, 8), fontsize=14):
    # Setting the folder within which we have the output
    output_direc = figure_direc + 'General/'
    utils.check_direc_exist(output_direc)

    # Getting the input location data
    file_dict = scenario.file_dict
    input_dict = file_dict['STARTFILES_filename']


    print_statement = 'We have {} particles in this release situation'.format(len(np.load(input_dict['lat'])))
    utils.print_statement(print_statement, to_print=True)

    # Now we have all the unique input locations, and the number of particles that are released there
    df = pd.DataFrame({'lat': np.load(input_dict['lat']), 'lon': np.load(input_dict['lon'])})
    df = df.groupby(['lat', 'lon']).size().reset_index().rename(columns={0: 'count'})

    # Getting the spatial domain for the figure
    spatial_domain = np.nanmin(file_dict['LON']), np.nanmax(file_dict['LON']), \
                     np.nanmin(file_dict['LAT']), np.nanmax(file_dict['LAT'])

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

    # Plotting the data
    plt.scatter(df['lon'], df['lat'], s=df['count'], zorder=1000, edgecolor='r', facecolor='none')

    file_name = output_direc + 'InputScenario.png'
    plt.savefig(file_name, bbox_inches='tight')
