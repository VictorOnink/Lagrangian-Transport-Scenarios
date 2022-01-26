import settings as settings
import utils
import xarray as xr
import numpy as np
import pandas as pd


class parcels_to_bayesian:
    def __init__(self, file_dict):
        self.parallel_step = settings.PARALLEL_STEP
        self.file_dict = file_dict
        self.scenario_name = settings.SCENARIO_NAME
        assert self.scenario_name in ['SizeTransport'], "The max distance function is not set up for {}".format(
            self.scenario_name)
        self.domain = (-180, 180, -90, 90)  # lon_min, lon_max, lat_min, lat_max
        self.temp_direc, self.output_direc = self.get_directories()
        # Creating the output dict
        # self.output_dict = self.create_output_dict()

    def get_directories(self):
        temp_direc = settings.SCRATCH_DIR
        output_direc = utils.get_output_directory(server=settings.SERVER) + 'bayes/{}/'.format(self.scenario_name)
        utils.check_direc_exist(temp_direc)
        utils.check_direc_exist(output_direc)
        return temp_direc, output_direc

    def source_cluster(self):
        # Load the starting locations of all of the particles
        parcels_dataset = xr.load_dataset(self.file_dict[settings.RUN][0])

        # Group by the starting location and determine the number of particles released at each site
        df = pd.DataFrame.from_dict({'particle_release_lon': parcels_dataset.lon[:, 0],
                                     'particle_release_lat': parcels_dataset.lat[:, 0],
                                     'weight': np.ones(parcels_dataset.lat[:, 0].shape)})
        release_sites = df.groupby(['particle_release_lon', 'particle_release_lat']).count()
        print(release_sites)

    def run(self):
        pass

