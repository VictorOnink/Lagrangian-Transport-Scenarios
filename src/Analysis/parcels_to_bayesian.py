import settings as settings
import utils
import xarray as xr
import numpy as np
import pandas as pd
from advection_scenarios import advection_files
import matplotlib.pyplot as plt
import visualization.visualization_utils as vUtils
import cmocean.cm as cmo
import matplotlib.colors as colors


class parcels_to_bayesian:
    def __init__(self, file_dict):
        self.parallel_step = settings.PARALLEL_STEP
        self.file_dict = file_dict
        self.scenario_name = settings.SCENARIO_NAME
        assert self.scenario_name in ['SizeTransport'], "The max distance function is not set up for {}".format(
            self.scenario_name)
        # size of the clustering in degrees lat/lon
        self.cluster_size = 1
        # Setting directories
        self.temp_direc, self.output_direc = self.get_directories()
        # Getting the LON/LAT/GRID arrays
        self.LON, self.LAT, self.GRID = self.get_lon_lat()
        # Get the cluster locations and ID for each particle which cluster they are from
        self.release_cluster, self.cluster_dict = self.source_cluster()
        self.cluster_number = np.nanmax(self.release_cluster)
        self.particle_number = self.release_cluster.size
        # Creating the output dictionary
        self.output_dict = self.create_output_dict()
        self.file_name = self.get_filename()

    def run(self):
        if self.parallel_step == 1:
            # First, we load the particle trajectories
            parcels_dataset = xr.load_dataset(self.file_dict[settings.RUN][settings.RESTART])
            particle_lon, particle_lat = parcels_dataset.lon, parcels_dataset.lat

            # Now, we loop through the various release clusters and calculate the annual averaged concentration for each
            # cluster
            for cluster_id in range(self.cluster_number.astype(int)):
                # Getting the trajectories for this cluster
                select_cluster = np.where(self.release_cluster == cluster_id)[0]
                lon_selection = particle_lon[select_cluster, :]
                lat_selection = particle_lat[select_cluster, :]
                # Flattening the arrays
                lon_selection = lon_selection.values.flatten()
                lat_selection = lat_selection.values.flatten()
                # Calculate the concentrations over all these positions
                concentration, binned_LAT, binned_LON = utils.histogram(lon_data=lon_selection, lat_data=lat_selection,
                                                                        bins_Lon=self.LON, bins_Lat=self.LAT)
                # Now, get all locations where the concentrations are greater than zero
                indices_LAT, indices_LON = np.where(concentration > 0)
                # Finally, adding the nonzero concentration points to the output dict
                for ind_la, ind_lo in zip(indices_LAT, indices_LON):
                    if (ind_la, ind_lo) in self.output_dict.keys():
                        self.output_dict[(ind_la, ind_lo)][cluster_id] = concentration[ind_la, ind_lo]
                        self.output_dict[(ind_la, ind_lo)]['total'] += concentration[ind_la, ind_lo]
                    else:
                        self.output_dict[(ind_la, ind_lo)] = {}
                        self.output_dict[(ind_la, ind_lo)][cluster_id] = concentration[ind_la, ind_lo]
                        self.output_dict[(ind_la, ind_lo)]['total'] = concentration[ind_la, ind_lo]
            # Saving the output
            utils.save_obj(self.file_name, self.output_dict)

        elif self.parallel_step == 2:
            utils.print_statement('Nothing happens for parcels_to_bayesian when settings.PARALLEL_STEP == 2',
                                  to_print=True)
        else:
            ValueError('settings.PARALLEL_STEP can not have a value of {}'.format(self.parallel_step))

    def get_filename(self):
        prefix = 'bayesian'
        file_name = utils.analysis_save_file_name(input_file=self.file_dict[settings.RUN][settings.RESTART],
                                                  prefix=prefix)
        return self.output_direc + file_name

    def get_directories(self):
        temp_direc = settings.SCRATCH_DIREC
        output_direc = settings.DATA_OUTPUT_DIREC + 'bayes/{}/'.format(self.scenario_name)
        utils.check_direc_exist(temp_direc)
        utils.check_direc_exist(output_direc)
        return temp_direc, output_direc

    def create_output_dict(self):
        output_dict = {'LON': 0.5 * self.LON[1:] + 0.5 * self.LON[:-1], 'LAT': 0.5 * self.LAT[1:] + 0.5 * self.LAT[:-1],
                       'cluster_dict': self.cluster_dict, 'norm': None, 'p_number': self.particle_number}
        return output_dict

    def source_cluster(self):
        # Load the starting locations of all of the particles
        parcels_dataset = xr.load_dataset(self.file_dict[settings.RUN][0])
        p_lat = parcels_dataset.lat[:, 0]
        p_lon = parcels_dataset.lon[:, 0]

        # Create larger clusters, where we have 2x2 degree bins. For this, we first determine all the sites on this 2x2
        # grid that have input sites
        cluster_lon, cluster_lat = self.cluster_lon_lat()

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
        return release_cluster, cluster_dict

    @staticmethod
    def get_lon_lat():
        adv_file_dict = advection_files.AdvectionFiles().file_names
        LON, LAT, GRID = adv_file_dict['LON'], adv_file_dict['LAT'], adv_file_dict['GRID']
        return LON, LAT, GRID

    def cluster_lon_lat(self):

        # Get the lon/lat arrays for the clustering
        cluster_lon = np.arange(np.floor(self.LON.min()) - self.cluster_size, np.floor(self.LON.max()) + self.cluster_size,
                                step=self.cluster_size)
        cluster_lat = np.arange(np.floor(self.LAT.min()) - self.cluster_size, np.floor(self.LAT.max()) + self.cluster_size,
                                step=self.cluster_size)

        return cluster_lon, cluster_lat

    def plot_clusters(self):
        # First, we need to set a number of parameters
        figure_direc = settings.FIGURE_OUTPUT_DIREC + '{}/'.format(settings.SCENARIO_NAME)
        figure_size = (10, 8)
        figure_shape = (1, 1)
        ax_label_size = 14
        adv_file_dict = advection_files.AdvectionFiles().file_names
        LON, LAT = adv_file_dict['LON'], adv_file_dict['LAT']
        spatial_domain = np.nanmin(LON), np.nanmax(LON), np.nanmin(LAT), np.nanmax(LAT)

        # Next, we create a map of our region
        fig = plt.figure(figsize=figure_size)
        gs = fig.add_gridspec(nrows=figure_shape[0], ncols=figure_shape[1] + 1, width_ratios=[1, 0.1])

        ax = []
        for rows in range(figure_shape[0]):
            for columns in range(figure_shape[1]):
                ax.append(vUtils.cartopy_standard_map(fig=fig, gridspec=gs, row=rows, column=columns,
                                                      domain=spatial_domain, lat_grid_step=5, lon_grid_step=10,
                                                      resolution='10m'))

        # Creating a cmap for the scatter plot, and creating a legend
        cmap_name = cmo.thermal
        norm = colors.LogNorm(vmin=1e-5, vmax=1e0)
        cbar_label, extend = r"P($S_i$)", None
        cmap = plt.cm.ScalarMappable(cmap=cmap_name, norm=norm)
        cax = fig.add_subplot(gs[:, -1])
        cbar = plt.colorbar(cmap, cax=cax, orientation='vertical', extend=extend)
        cbar.set_label(cbar_label, fontsize=ax_label_size)

        # Plotting the release midpoints
        lon_site, lat_site, number_site = [], [], []
        for site in self.cluster_dict.keys():
            la, lo, n = self.cluster_dict[site]
            lon_site.append(lo)
            lat_site.append(la)
            number_site.append(n / self.release_cluster.size)
        ax[0].scatter(lon_site, lat_site, c=number_site, cmap=cmap_name, norm=norm, zorder=1000, s=20)

        # Saving the figure
        utils.print_statement('This uses {} clusters'.format(self.cluster_number), to_print=True)
        plt.savefig(figure_direc + 'General/Input_clusters.png', bbox_inches='tight')




