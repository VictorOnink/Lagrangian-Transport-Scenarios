import settings as settings
import utils
import xarray as xr
import numpy as np
import pandas as pd
from advection_scenarios import advection_files


class parcels_to_bayesian:
    def __init__(self, file_dict):
        self.parallel_step = settings.PARALLEL_STEP
        self.file_dict = file_dict
        self.scenario_name = settings.SCENARIO_NAME
        assert self.scenario_name in ['SizeTransport'], "The max distance function is not set up for {}".format(
            self.scenario_name)
        self.cluster_size = 2  # size of the clustering in degrees lat/lon
        self.temp_direc, self.output_direc = self.get_directories()
        # Get the cluster locations and ID for each particle which cluster they are from
        self.release_cluster, self.cluster_dict = self.source_cluster()

    def run(self):
        pass

    def get_directories(self):
        temp_direc = settings.SCRATCH_DIR
        output_direc = utils.get_output_directory(server=settings.SERVER) + 'bayes/{}/'.format(self.scenario_name)
        utils.check_direc_exist(temp_direc)
        utils.check_direc_exist(output_direc)
        return temp_direc, output_direc

    def source_cluster(self):
        # Load the starting locations of all of the particles
        parcels_dataset = xr.load_dataset(self.file_dict[settings.RUN][0])
        p_lat = parcels_dataset.lat[:, 0]
        p_lon = parcels_dataset.lon[:, 0]

        # Group by the starting location and determine the number of particles released at each site
        df = pd.DataFrame.from_dict({'particle_release_lon': p_lon, 'particle_release_lat': p_lat,
                                     'weight': np.ones(p_lat.shape)})
        release_sites = df.groupby(['particle_release_lon', 'particle_release_lat']).count()

        # Get the cluster lon, lat and weight
        release_lon = release_sites.index.get_level_values('particle_release_lon').values
        release_lat = release_sites.index.get_level_values('particle_release_lat').values
        release_weight = release_sites['weight'].values

        # Create larger clusters, where we have 2x2 degree bins. For this, we first determine all the sites on this 2x2
        # grid that have input sites
        cluster_lon, cluster_lat = self.cluster_lon_lat()

        cluster_grid, _, _ = utils.histogram(lon_data=release_lon, lat_data=release_lat, bins_Lon=cluster_lon,
                                             bins_Lat=cluster_lat, weight_data=release_weight, area_correc=False)

        # Now, we create an array that for each particle (not each release site, but each particle) will indicate the
        # release cluster from which it originates
        release_cluster = np.zeros(parcels_dataset.lon[:, 0].shape)
        cluster_index = 0
        cluster_dict = {}
        for index_lon in range(0, cluster_lon.size - 1):
            for index_lat in range(0, cluster_lat.size - 1):
                lat_min, lat_max = cluster_lat[index_lat], cluster_lat[index_lat + 1]
                lon_min, lon_max = cluster_lon[index_lon], cluster_lon[index_lon + 1]
                selection = np.logical_and(np.logical_and(lat_min <= p_lat, p_lat < lat_max), np.logical_and(lon_min <= p_lon, p_lon < lon_max))
                if np.nanmax(selection) == 1:
                    release_cluster[selection] = cluster_index
                    # Within the cluster dict we save the mid lon/lat for each cluster, as well as the number of
                    # particles
                    cluster_dict[cluster_index] = ((lat_min + lat_max) / 2, (lon_min + lon_max) / 2, np.nansum(selection))
                    cluster_index += 1
        print(cluster_dict)
        return release_cluster, cluster_dict

    def cluster_lon_lat(self):
        # Load the dimensions
        adv_file_dict = advection_files.AdvectionFiles().file_names
        LON, LAT = adv_file_dict['LON'], adv_file_dict['LAT']

        # Get the lon/lat arrays for the clustering
        cluster_lon = np.arange(np.floor(LON.min()) - self.cluster_size, np.floor(LON.max()) + self.cluster_size,
                                step=self.cluster_size)
        cluster_lat = np.arange(np.floor(LAT.min()) - self.cluster_size, np.floor(LAT.max()) + self.cluster_size,
                                step=self.cluster_size)

        return cluster_lon, cluster_lat

