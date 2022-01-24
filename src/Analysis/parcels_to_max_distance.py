import settings as settings
import utils
import xarray as xr
import numpy as np
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

                # If the simulation is longer than 1 year, looping through all subsequent years and determining the
                # overall max distance from shore for each particle
                if settings.SIM_LENGTH > 1:
                    for restart in range(1, settings.SIM_LENGTH):
                        restart_dataset = xr.load_dataset(self.file_dict[settings.RUN][restart])
                        print(restart_dataset)
                        print(parcels_dataset)
                        stack_max = np.vstack((self.output_dict['max_distance'],
                                               parcels_dataset.distance2coast.max(dims='obs', skipna=False)))
                        self.output_dict['max_distance'] = np.nanmax(stack_max, axis=1)
                        print(self.output_dict['max_distance'].shape)

                # Save the starting lon and lat positions of each particle
                self.output_dict['particle_release_lon'] = parcels_dataset.lon[:, 0]
                self.output_dict['particle_release_lat'] = parcels_dataset.lat[:, 0]

                # Create a dataframe with the starting lon/lat and the max distance from shore
                df = pd.DataFrame.from_dict({'particle_release_lon': self.output_dict['particle_release_lon'],
                                             'particle_release_lat': self.output_dict['particle_release_lat'],
                                             'max_distance': self.output_dict['max_distance']})

                # Group the particles by release location and calculate the median max distance for each release site
                median = df.groupby(['particle_release_lon', 'particle_release_lat']).median()

                # For each release site, save the lon/lat and median maximum distance from shore
                self.output_dict['site_lon'] = median.index.get_level_values('particle_release_lon').values
                self.output_dict['site_lat'] = median.index.get_level_values('particle_release_lat').values
                self.output_dict['median_max_distance'] = median['max_distance'].values

                # Saving the output
                utils.save_obj(filename=self.file_name(), item=self.output_dict)

        @staticmethod
        def create_output_dict():
            # Create the output dict
            output_dict = {"max_distance": None, "particle_release_lon": None, "particle_release_lat": None,
                           "median_max_distance": None, "site_lon": None, "site_lat": None}
            return output_dict

        def get_directories(self):
            temp_direc = settings.SCRATCH_DIR
            output_direc = utils.get_output_directory(server=settings.SERVER) + 'max_distance' \
                                                                                '/{}/'.format(self.scenario_name)
            utils.check_direc_exist(temp_direc)
            utils.check_direc_exist(output_direc)
            return temp_direc, output_direc

        def file_name(self):
            prefix = 'maximum_distance'
            file_name = utils.analysis_save_file_name(input_file=self.file_dict[settings.RUN][settings.RESTART],
                                                      prefix=prefix)
            return self.output_direc + file_name

