import settings
import utils
from netCDF4 import Dataset
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean
import pandas as pd
import advection_scenarios.create_input_files as create_input_files


class General_input_scenario:
    def __init__(self, scenario, figure_direc):
        # Scenario specific variables
        self.scenario = scenario
        self.figure_direc = figure_direc
        self.input_file_prefix = create_input_files.create_input_files(prefix=settings.ADVECTION_DATA, repeat_dt=None,
                                                                       grid=np.ones(1), lon=np.ones(1), lat=np.ones(1))
        # Data variables
        self.output_direc = figure_direc + 'General/'
        utils.check_direc_exist(self.output_direc)
        self.file_dict = self.scenario.file_dict
        # Figure variables
        self.figure_size = (10, 8)
        self.figure_shape = (1, 1)
        self.ax_label_size = 14
        self.ax_ticklabel_size = 12
        self.spatial_domain = np.nanmin(self.file_dict['LON']), np.nanmax(self.file_dict['LON']), \
                              np.nanmin(self.file_dict['LAT']), np.nanmax(self.file_dict['LAT'])

    def plot(self):
        # Loading the data
        input_lon, input_lat = np.array([]), np.array([])
        for run in range(settings.RUN_RANGE):
            lon_file_name = self.input_file_prefix + '{}_run={}.npy'.format('lon', run)
            lat_file_name = self.input_file_prefix + '{}_run={}.npy'.format('lat', run)
            input_lon = np.append(input_lon, np.load(lon_file_name))
            input_lat = np.append(input_lat, np.load(lat_file_name))

        print_statement = 'We have {} particles in this release situation'.format(len(input_lon))
        utils.print_statement(print_statement, to_print=True)

        # Now we have all the unique input locations, and the number of particles that are released there
        df = pd.DataFrame({'lat': input_lat, 'lon': input_lon})
        df = df.groupby(['lat', 'lon']).size().reset_index().rename(columns={0: 'count'})

        # Creating the base figure
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(nrows=self.figure_shape[0], ncols=self.figure_shape[1] + 1, width_ratios=[1, 0.1])

        ax = []
        for rows in range(self.figure_shape[0]):
            for columns in range(self.figure_shape[1]):
                ax.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                      domain=self.spatial_domain, lat_grid_step=5, lon_grid_step=10,
                                                      resolution='10m'))

        # Plotting the data
        if settings.INPUT == 'LebretonDivision':
            sizes = df['count']/df['count'].min() * 100
        elif settings.INPUT == 'Lebreton':
            sizes = df['count'] / df['count'].min() * 10
        ax[0].scatter(df['lon'], df['lat'], s=sizes, zorder=1000, edgecolor='r', facecolor='none')

        file_name = self.output_direc + 'InputScenario_{}_{}.png'.format(settings.ADVECTION_DATA, settings.INPUT)
        plt.savefig(file_name, bbox_inches='tight')
