import settings as settings
import utils
from advection_scenarios import advection_files
import xarray as xr
import numpy as np
from scipy import io
import progressbar
import pandas as pd


class parcels_to_max_distance:
        def __init__(self, file_dict):
            self.parallel_step = settings.PARALLEL_STEP
            self.file_dict = file_dict
            self.scenario_name = settings.SCENARIO_NAME
            assert self.scenario_name in ['SizeTransport'], "The max distance function is not set up for {}".format(self.scenario_name)
            self.domain = (-180, 180, -90, 90)  # lon_min, lon_max, lat_min, lat_max
            self.temp_direc, self.output_direc = self.get_directories()
            # Creating the output dict
            self.output_dict = self.create_output_dict()

        def run(self):
            if self.parallel_step == 1:
                utils.print_statement('Nothing happens for parcels_to_max_distance when settings.PARALLEL_STEP == 1',
                                      to_print=True)
            elif self.parallel_step == 2:
                # Calculate the max distance over the first year of the simulation
                parcels_dataset = xr.load_dataset(self.file_dict[settings.RUN][0])
                self.output_dict['max_distance'] = parcels_dataset.distance2coast.max(dim='obs', skipna=True)

                # If the simulation is longer than 1 year, looping through all subsequent years and determining the overall
                # max distance from shore for each particle
                if settings.SIM_LENGTH > 1:
                    for restart in range(1, settings.SIM_LENGTH):
                        restart_dataset = xr.load_dataset(self.file_dict[settings.RUN][restart])
                        stack_max = np.vstack((self.output_dict['max_distance'], restart_dataset.distance2coast.max(dims='obs', skipna=True)))
                        self.output_dict['max_distance'] = stack_max

                # Save the starting lon and lat positions of each particle
                self.output_dict['release_lon'] = parcels_dataset.lon[:, 0]
                self.output_dict['release_lat'] = parcels_dataset.lat[:, 0]

                # Create a dataframe with the starting lon/lat and the max distance from shore
                df = pd.DataFrame.from_dict({'release_lon': self.output_dict['release_lon'],
                                             'release_lat': self.output_dict['release_lat'],
                                             'max_distance': self.output_dict['max_distance']})

                # Group the particles by release location and calculate the median max distance for each release site
                median = df.groupby(['release_lon', 'release_lat']).median()

                # Create a dictionary contains all the release sites and median distance
                release_dict, index = {}, 0
                for release_lon, item in median:
                    lon_sub = median.get_group(release_lon)
                    for release_lat, item in lon_sub:
                        print(release_lon, release_lat)
                        print(lon_sub[release_lat]['max_distance'])
                        print('\n')




        @staticmethod
        def create_output_dict():
            # Set array onto which we save the median max distance from shore from each release point, along with the LON
            # and LAT arrays
            advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                                advection_scenario=settings.ADVECTION_DATA,
                                                                repeat_dt=None)
            adv_file_dict = advection_scenario.file_names
            LON, LAT, GRID = adv_file_dict['LON'], adv_file_dict['LAT'], adv_file_dict['GRID']
            median_max_distance = np.zeros(GRID.shape, dtype=float)

            # Create the output dict
            output_dict = {"max_distance": None, "release_lon": None, "release_lat": None,
                           "median_max_distance": median_max_distance, "LON": LON, "LAT": LAT}
            return output_dict

        def get_directories(self):
            temp_direc = settings.SCRATCH_DIR
            output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/{}/'.format(self.scenario_name)
            utils.check_direc_exist(temp_direc)
            utils.check_direc_exist(output_direc)
            return temp_direc, output_direc

        def file_name(self):
            prefix = 'maximum_distance'
            file_name = utils.analysis_save_file_name(input_file=self.file_dict[settings.RUN][settings.RESTART],
                                                      prefix=prefix)
            return self.output_direc + file_name



# def parcels_to_max_distance(file_dict: dict, lon_min: float = -180, lon_max: float = 180, lat_min: float = -90,
#                             lat_max: float = 90):
#     domain = [lon_min, lon_max, lat_min, lat_max]
#     output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/'
#
#     # loop through the runs
#     for run in progressbar.progressbar(range(settings.RUN_RANGE)):
#         # Loop through the restart files
#         for restart in range(settings.SIM_LENGTH):
#             # Load the lon, lat, time, beach and weight data
#             parcels_file = file_dict[run][restart]
#             dataset = Dataset(parcels_file)
#             lon, lat = dataset.variables['lon'][:, :-1], dataset.variables['lat'][:, :-1]
#             distance = dataset.variables['distance'][:, :-1]
#             # Get the maximum distance over the course of each trajectory
#             max_distance_file = np.nanmax(distance, axis=1, keepdims=True)
#             if restart == 0:
#                 # If this is the beginning of the run, we want to see if the initial positions of the particles fall
#                 # within the area we have defined as the domain
#                 within_domain = utils.particles_in_domain(domain=domain, lon=lon[:, :1], lat=lat[:, :1])
#             # Making a cumulative array for the run
#             if restart == 0:
#                 max_distance_file = max_distance_file[within_domain]
#                 max_distance_run = max_distance_file.reshape((len(max_distance_file), 1))
#             else:
#                 max_distance_file = max_distance_file[within_domain]
#                 max_distance_run = np.hstack((max_distance_run, max_distance_file.reshape((len(max_distance_file), 1))))
#         if run == 0:
#             max_distance_simulation = np.nanmax(max_distance_run, axis=1, keepdims=True)
#         else:
#             max_distance_simulation = np.vstack(
#                 (max_distance_simulation, np.nanmax(max_distance_run, axis=1, keepdims=True)))
#     # Get the output dictionary
#     output_dict = {'max_dist': max_distance_simulation}
#     # Saving the computed concentration, where we have a few default names for the prefix
#     if lon_min == -180 and lon_max == 180 and lat_min == -90 and lat_max == 90:
#         prefix = 'maximum_distance-global'
#     else:
#         prefix = 'maximum_distance-lon_min={}-lon_max={}-lat_min={}-lat_max={}'.format(lon_min, lon_max, lat_min,
#                                                                                        lat_max)
#
#     output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
#     io.savemat(output_name, output_dict)
#     utils.print_statement("The maximum distance has been saved")
