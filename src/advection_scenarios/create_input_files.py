import settings
import utils
from advection_scenarios.create_distance_to_shore import create_distance_to_shore_land
import numpy as np
import glob
from netCDF4 import Dataset
import pandas as pd
import geopy.distance
from datetime import timedelta
import math
import fiona
from shapely import vectorized
from shapely.geometry import shape
import xarray
import progressbar
from parcels import Field


class create_input_files:
    def __init__(self, advection_prefix: str, grid: np.array, lon: np.array, lat: np.array, repeat_dt: timedelta):
        # Check if the directory for the input files exists
        utils.check_direc_exist(settings.INPUT_DIREC)
        # Input parameters
        self.input_prefix = None
        self.advection_prefix = advection_prefix
        self.repeat_dt = repeat_dt
        self.lon_inputs, self.lat_inputs, self.plastic_inputs = None, None, None
        self.inputs_grid = None
        self.distance = None
        self.coastal = None
        self.inputs_coastal_grid = None
        self.releases = None
        self.particle_number, self.particle_weight, self.particle_remain_num, self.particle_remain = None, None, None, None
        self.particle_lat, self.particle_lon, self.particle_weight = None, None, None
        # Grid parameters
        self.GRID = grid
        self.LON = lon
        self.LAT = lat

    def create_files(self):
        # Get the prefix for the input files, then check if the files exist
        self.input_prefix = self.get_input_prefix()
        if len(glob.glob(self.input_prefix + "*")) > 0:
            utils.print_statement("The input files {} are already present".format(self.input_prefix))
            return self.input_prefix
        else:
            utils.print_statement("We need to create the input files {}".format(self.input_prefix))

            # Calculate the number of particle releases per year
            self.releases = self.number_of_releases()

            # Get the unprocessed lon/lat coordinates of all the plastic sources
            if settings.INPUT in ["Jambeck"]:
                self.lon_inputs, self.lat_inputs, self.plastic_inputs, self.input_max, self.input_min, self.distance = self.get_lon_lat_weights_Jambeck()
            elif settings.INPUT in ["Lebreton", "LebretonDivision", "LebretonKaandorpInit"]:
                self.lon_inputs, self.lat_inputs, self.plastic_inputs, self.input_max, self.input_min = self.get_lon_lat_weights_Lebreton()
            elif settings.INPUT in ["Point_Release"]:
                self.lon_inputs, self.lat_inputs, self.plastic_inputs, self.input_max, self.input_min = self.get_lon_lat_weights_PointRelease()
            elif settings.INPUT in ["Uniform"]:
                self.get_lon_lat_weights_Uniform()
                return self.input_prefix

            # Only keep particles within the domain
            self.lon_inputs, self.lat_inputs, self.plastic_inputs = self.within_domain()
            utils.print_statement("We have {} input locations".format(len(self.lon_inputs)), to_print=True)

            # Get the inputs onto the grid of the advection data
            self.inputs_grid = self.histogram(lon_data=self.lon_inputs, lat_data=self.lat_inputs,
                                              weight_data=self.plastic_inputs)

            # Remove any sources further than 50 kilometers from land
            if settings.INPUT in ["Jambeck"]:
                print(self.inputs_grid.shape)
                print(self.distance.shape)
                self.inputs_grid[self.distance > 50] = 0

            # Getting all ocean cells that are adjacent to coastal ocean cells (which are in turn adjacent to land). The
            # particles will be placed one cell away from land. For brevity, cells one removed from land are referred to
            # as coastal within this code
            self.coastal = self.get_coastal_cells()

            # Moving all plastic sources onto the coastal cells, and the computing the input per release_dt
            self.inputs_coastal_grid = self.input_to_nearest_coastal()
            self.inputs_coastal_grid /= self.releases()

            # Calculating the number of particles that are released per cell, along with the weights of these cells
            self.particle_number, self.particle_weight, self.particle_remain_num, self.particle_remain = self.number_weights_releases()

            # Getting arrays of all the individual particles and their weights
            self.particle_lat, self.particle_lon, self.particle_weight = self.particle_grid_to_list()

            # Determine how much of the total plastic is accounted for with the given self.input_max and self.input_min
            total_input = np.sum(self.inputs_coastal_grid)
            missing_percent = np.divide(total_input - np.sum(self.particle_weight), total_input) * 100
            utils.print_statement("The particles account for {}% of the total inputs".format(100 - missing_percent),
                                  to_print=True)

            # Dividing the particles into runs
            run_number = self.split_to_runs(particle_lat=self.particle_lat, particle_lon=self.particle_lon,
                                            particle_weight=self.particle_weight, input_prefix=self.input_prefix)
            utils.print_statement("The {} input files have been created.".format(run_number), to_print=True)

    def get_input_prefix(self):
        """
        Set the prefix of the input files
        :return:
        """
        if settings.INPUT in ["Jambeck", "Lebreton", "LebretonDivision", "LebretonKaandorpInit"]:
            return settings.INPUT_DIREC + settings.INPUT + "_{}_{}_".format(self.advection_prefix, settings.STARTYEAR)
        elif settings.INPUT in ["Point_Release"]:
            str_format = (self.advection_prefix, settings.STARTYEAR, settings.INPUT_LAT, settings.INPUT_LON)
            return settings.INPUT_DIREC + settings.INPUT + "_{}_{}_{}_{}_".format(*str_format)
        elif settings.INPUT == "Uniform":
            str_format = (self.advection_prefix, settings.STARTYEAR)
            return settings.INPUT_DIREC + settings.INPUT + "_{}_{}_".format(*str_format)
        else:
            ValueError("Perhaps take another look at what input you are using? {} is not valid".format(settings.INPUT))

    def number_of_releases(self):
        if self.repeat_dt is None:
            releases = 1
        else:
            releases = math.floor(timedelta(days=365) / self.repeat_dt) + 1
        return releases

    def get_lon_lat_weights_Jambeck(self) -> tuple:
        """
        Getting the lon/lat/weights data for the Jambeck input scenario
        :return:
        """
        # Get the population data
        dataset = Dataset(settings.INPUT_DIREC + "gpw_v4_population_count_adjusted_rev11_2pt5_min.nc")
        var_name = "UN WPP-Adjusted Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes"
        population = np.array(dataset.variables[var_name][2, :, :])
        lon_population = dataset.variables["longitude"][:]
        lat_population = dataset.variables["latitude"][:]
        Lon_population, Lat_population = np.meshgrid(lon_population, lat_population)
        # Get the mismanaged plastic fraction
        mismanaged = self.get_mismanaged_fraction_Jambeck(dataset=dataset)
        # Get the distance from land to shore
        distance_file = settings.INPUT_DIREC + self.advection_prefix + "_distance_to_coast_land.nc"
        utils.print_statement(distance_file)
        distance = self.get_distance_to_shore(filename=distance_file)
        # The yearly mismanaged plastic
        transport_to_ocean = 0.15  # percentage of mismanaged plastic that reaches the ocean
        kg_to_tons = 1000
        mismanaged_total = np.multiply(mismanaged, population) * 365 * transport_to_ocean / kg_to_tons
        # Get everything in column arrays to load
        lon_inputs = Lon_population[mismanaged_total > 0].flatten()
        lat_inputs = Lat_population[mismanaged_total > 0].flatten()
        plastic_inputs = mismanaged_total[mismanaged_total > 0].flatten()
        # Get the max/min mass for each parcels particle
        input_max, input_min = settings.INPUT_MAX, settings.INPUT_MIN
        return lon_inputs, lat_inputs, plastic_inputs, input_max, input_min, distance

    def get_distance_to_shore(self, filename: str):
        if utils.check_file_exist(filename):
            utils.print_statement("The distance to shore file already exists")
        else:
            create_distance_to_shore_land(output_name=filename, grid=self.GRID, lon=self.LON, lat=self.LAT)
        return self.slicing_correction(Dataset(filename).variables["distance"][:])

    @staticmethod
    def get_lon_lat_weights_Lebreton() -> tuple:
        """
        Getting the lon/lat/weights data for input scenarios based on the Lebreton et al. (2017) river input estimates
        :return:
        """
        lebData = pd.read_csv(settings.INPUT_DIREC + "PlasticRiverInputs.csv")
        lon_inputs = np.array(lebData["X"])
        lat_inputs = np.array(lebData["Y"])
        plastic_inputs = np.array(lebData["i_low"])
        # Get the max/min mass for each parcels particle
        input_max, input_min = settings.INPUT_MAX, settings.INPUT_MIN
        return lon_inputs, lat_inputs, plastic_inputs, input_max, input_min

    @staticmethod
    def get_lon_lat_weights_PointRelease(self) -> tuple:
        """
        Getting the lon/lat/weights data for point release scenarios
        :return:
        """
        lon_inputs = np.ones((settings.PARTICLE_NUMBER, 1)) * settings.INPUT_LON
        lat_inputs = np.ones((settings.PARTICLE_NUMBER, 1)) * settings.INPUT_LAT
        plastic_inputs = np.ones((settings.PARTICLE_NUMBER, 1))
        # Get the max/min mass for each parcels particle
        input_max, input_min = 1.0, 0.0
        return lon_inputs, lat_inputs, plastic_inputs, input_max, input_min

    def get_lon_lat_weights_Uniform(self):
        """
        Getting the lon, lat positions for a uniform particle release, assuming all particle weights are even
        :return:
        """
        lon_min, lon_max, lat_min, lat_max = np.min(self.LON), np.max(self.LON), np.min(self.LAT), np.max(self.LAT)
        release_grid = np.mgrid[lon_min:lon_max:settings.RELEASE_GRID, lat_min:lat_max:settings.RELEASE_GRID]
        n = release_grid[0].size
        lon_inputs, lat_inputs = np.reshape(release_grid[0], n), np.reshape(release_grid[1], n)
        # Remove cells on land
        Land = self.land_mark_checker()
        [lon_inputs, lat_inputs] = [np.array([lo for lo, la in zip(lon_inputs, lat_inputs) if Land[0, 0, la, lo] == 0.0]),
                                    np.array([la for lo, la in zip(lon_inputs, lat_inputs) if Land[0, 0, la, lo] == 0.0])]
        _ = self.split_to_runs(particle_lat=lat_inputs, particle_lon=lon_inputs, particle_weight=None,
                               output_prefix=self.input_prefix)

    def land_mark_checker(self):
        mask = np.ma.getmask(self.GRID)
        land = Field("Land", mask, lon=self.lon, lat=self.lat, transpose=False, mesh="spherical")
        return land

    def within_domain(self):
        lon_max, lon_min = np.max(self.LON), np.min(self.LON)
        lat_max, lat_min = np.max(self.LAT), np.min(self.LAT)
        # Check which cells are within the domain and non-zero inputs
        domain = (self.lon_inputs <= lon_max) & (self.lon_inputs >= lon_min) & (self.lat_inputs >= lat_min) & \
                 (self.lat_inputs <= lat_max) & (self.plastic_inputs > 0)
        str_format = (len(self.lon_inputs), settings.INPUT, np.sum(domain * 1))
        utils.print_statement(
            "Of the original {} input sites in the {} scenario, {} are within the domain".format(*str_format))
        return self.lon_inputs[domain], self.lat_inputs[domain], self.plastic_inputs[domain]

    def histogram(self, lon_data, lat_data, weight_data):
        masses = np.zeros(self.GRID.shape)
        counts = np.zeros(self.GRID.shape)
        for i in range(np.array(lon_data).shape[0]):
            if weight_data[i] > 0:
                lat_selec = np.argmin(np.abs(lat_data[i] - self.LAT))
                lon_selec = np.argmin(np.abs(lon_data[i] - self.LON))
                masses[lat_selec, lon_selec] += weight_data[i]
                counts[lat_selec, lon_selec] += 1
        return masses  # weight / km^2

    def get_cells_adjacent_to_land(self) -> np.array:
        """
        This returns an array where all ocean cells that are directly touching (horizontally, vertically or adjacent)
        a land cell are marked with 1, all other cells are 0
        :return:
        """
        mask = self.GRID.mask
        coastal = np.zeros(mask.shape, dtype=bool)
        for lat_index in range(mask.shape[0]):
            for lon_index in range(mask.shape[1]):
                if not mask[lat_index, lon_index]:
                    for lat_step in [-1, 0, 1]:
                        for lon_step in [-1, 0, 1]:
                            if mask[(lat_index + lat_step) % mask.shape[0], (lon_index + lon_step) % mask.shape[1]]:
                                coastal[lat_index, lon_index] = True
        return coastal

    def get_coastal_cells(self) -> np.array:
        """
        This returns an array where all ocean cells that are directly touching (horizontally, vertically or adjacent)
        a land-adjacent ocean cell cell are marked with 1, all other cells are 0. The land adjacent ocean cells are
        computed in self.get_cells_adjacent_to_land()
        :return:
        """
        mask = self.GRID.mask
        coastal = self.get_cells_adjacent_to_land()
        # We want the next layer of ocean cells after the coastal cells
        coastal_extended = np.zeros(mask.shape, dtype=bool)
        for lat_index in range(mask.shape[0]):
            for lon_index in range(mask.shape[1]):
                if not mask[lat_index, lon_index] and coastal[lat_index, lon_index]:
                    for lat_step in [-1, 0, 1]:
                        for lon_step in [-1, 0, 1]:
                            if coastal[(lat_index + lat_step) % mask.shape[0], (lon_index + lon_step) % mask.shape[1]]:
                                coastal_extended[lat_index, lon_index] = True
        return coastal_extended

    def input_to_nearest_coastal(self):
        # Find the positions with non-zero inputs
        non_zero_inputs = np.where(self.inputs_grid > 0)
        inputs_coastal_grid = np.zeros(self.coastal.shape)
        # Dimensions of the self.coastal array
        N_lat, N_lon = self.coastal.shape
        # Looping through the non-zero points to find the nearest one with
        for point in range(len(non_zero_inputs[0])):
            lat_point, lon_point = non_zero_inputs[0][point], non_zero_inputs[1][point]
            # If the input is already in the coastal cell
            if self.coastal[lat_point, lon_point]:
                inputs_coastal_grid[lat_point, lon_point] += self.inputs_grid[lat_point, lon_point]
            # Else find the nearest coastal cell
            else:
                step = 1
                coastal_lon, coastal_lat, coastal_distance = [], [], []
                while len(coastal_lon) < 1:
                    for lat_step in [-step, step]:
                        for lon_step in range(-step, step + 1):
                            if self.coastal[(lat_point + lat_step) % N_lat, (lon_point + lon_step) % N_lon]:
                                coastal_lat.append((lat_point + lat_step) % N_lat)
                                coastal_lon.append((lon_point + lon_step) % N_lon)
                    step += 1

                # Now find the coastal cell that is geographically the closest to the input
                nearest_dist = geopy.distance.distance((self.lat[lat_point], self.lon[lon_point]),
                                                       (self.lat[coastal_lat[0]], self.lon[coastal_lon[0]]))
                nearest_lat, nearest_lon = coastal_lat[0], coastal_lon[0]
                if len(coastal_lon) > 0:
                    for candidate in range(len(coastal_lon)):
                        distance = geopy.distance.distance((self.lat[lat_point], self.lon[lon_point]),
                                                           (self.lat[coastal_lat[candidate]], self.lon[coastal_lon[candidate]]))
                        if distance < nearest_dist:
                            nearest_dist = distance
                            nearest_lat, nearest_lon = coastal_lat[candidate], coastal_lon[candidate]
                inputs_coastal_grid[nearest_lat, nearest_lon] += self.inputs_grid[lat_point, lon_point]
        return inputs_coastal_grid

    def number_weights_releases(self):
        particle_number = np.zeros(self.inputs_coastal_grid.shape)
        particle_weight = np.zeros(self.inputs_coastal_grid.shape)
        particle_remain = np.zeros(self.inputs_coastal_grid.shape)
        particle_remain_num = np.zeros(self.inputs_coastal_grid.shape)
        # First all the inputs lower than the cutoff and greater than minimum
        selection = (self.inputs_coastal_grid < self.input_max) & (self.inputs_coastal_grid >= self.input_min)
        particle_number[selection] += 1
        particle_weight[selection] += self.inputs_coastal_grid[selection]
        # Next we consider the cells where we have inputs greater than the cutoff
        selection = self.inputs_coastal_grid >= self.input_max
        particle_number[selection] += np.floor(self.inputs_coastal_grid[selection] / self.input_max)
        particle_weight[selection] += self.input_max
        # Finally, we have the remainder for when inputs are not multiples of the cutoff, and so we have a remainder
        selection = (self.inputs_coastal_grid - np.multiply(particle_number, particle_weight)) > self.input_min
        particle_remain_num[selection] += 1
        particle_remain[selection] += (self.inputs_coastal_grid[selection] - np.multiply(particle_number, particle_weight)[selection])
        return particle_number.astype("int"), particle_weight, particle_remain_num.astype("int"), particle_remain

    def particle_grid_to_list(self):
        particle_lat, particle_lon, particle_mass = [], [], []
        # Getting the dimensions of self.particle_number
        N_lat, N_lon = self.particle_number.shape
        for lat_index in range(N_lat):
            for lon_index in range(N_lon):
                if self.particle_number[lat_index, lon_index] > 0:
                    for reps in range(self.particle_number[lat_index, lon_index]):
                        particle_lat.append(self.lat[lat_index])
                        particle_lon.append(self.lon[lon_index])
                        particle_mass.append(self.particle_weight[lat_index, lon_index])
                if self.particle_remain_num[lat_index, lon_index] > 0:
                    particle_lat.append(self.lat[lat_index])
                    particle_lon.append(self.lon[lon_index])
                    particle_mass.append(self.particle_remain[lat_index, lon_index])
        particle_lat, particle_lon, particle_mass = np.array(particle_lat), np.array(particle_lon), np.array(particle_mass)
        non_zero = particle_mass > 0
        return particle_lat[non_zero], particle_lon[non_zero], particle_mass[non_zero]

    @staticmethod
    def get_mismanaged_fraction_Jambeck(dataset: Dataset):
        mismanaged_file = settings.INPUT_DIREC + "Jambeck_mismanaged_grid.nc"
        if utils.check_file_exist(mismanaged_file):
            utils.print_statement("The mismanaged grid already exists")
            return Dataset(mismanaged_file).variables["mismanaged_plastic"][:]
        else:
            utils.print_statement("We need to generate the mismanaged grid")
            # Load the grid of the population data
            lon_pop, lat_pop = dataset.variables["longitude"][:], dataset.variables["latitude"][:]
            Lon, Lat = np.meshgrid(lon_pop, lat_pop)
            # Load the Jambeck estimates of mismanaged plastic per capita
            jambeck_excel = pd.read_excel(settings.INPUT_DIREC + "Jambeck2010data.xlsx")
            jambeck_country = list(jambeck_excel["Country"])
            jambeck_data = jambeck_excel["Mismanaged plastic waste [kg/person/day]7"]
            # Initialize grid for mismanaged plastic estimates
            mismanaged_grid = np.zeros(Lon.shape)
            # Getting the country shapefiles
            countries = fiona.open(settings.INPUT_DIREC +
                                   "country_shapefile/gpw_v4_national_identifier_grid_rev11_30_sec.shp")
            # Looping through all the countries, and marking which of the grid cells fall within which
            for country_index in progressbar.progressbar(range(len(countries))):
                if countries[country_index]["properties"]["NAME0"] in jambeck_country:
                    country_geometry = shape(countries[country_index]["geometry"])
                    country_mask = vectorized.contains(country_geometry, Lon, Lat)
                    if country_index in [145, 146]:
                        # 122 is the Netherlands Antilles value
                        mismanaged_grid[country_mask] = jambeck_data[122]
                    else:
                        mismanaged_grid[country_mask] = jambeck_data[
                            jambeck_country.index(countries[country_index]["properties"]["NAME0"])]
            # Saving the entire distance field
            utils.print_statement("Starting to save the mismanaged grid")
            coords = [("lat", lat_pop), ("lon", lon_pop)]
            misman = xarray.DataArray(mismanaged_grid, coords=coords)
            dcoo = {"lat": lat_pop, "lon": lon_pop}
            dset = xarray.Dataset({"mismanaged_plastic": misman}, coords=dcoo)
            dset.to_netcdf(mismanaged_file)
            utils.print_statement("The mismanaged grid has now been created and saved for future use")
            return mismanaged_grid

    @staticmethod
    def slicing_correction(array: np.array):
        """
        This assures that array is returned as 2D array
        :param array:
        :return:
        """
        if len(array.shape) == 2:
            return array
        elif len(array.shape) == 3:
            return array[0, :, :]
        elif len(array.shape) == 4:
            return array[0, 0, :, :]
        else:
            ValueError("What sort of data are you using that it has {} dimensions?".format(len(array.shape)))

    @staticmethod
    def split_to_runs(particle_lat: np.array, particle_lon: np.array, particle_weight: None,
                      input_prefix: str) -> int:
        """
        This splits all the particles within the particle_lat/particle_lon/particle_weight arrays into runs of
        settings.INPUT_DIV particles each.
        :param particle_lat: array containing all the Lon positions
        :param particle_lon: array containing all the Lat positions
        :param particle_weight: array containing the particle weights, but can be None if all particles are assumed to
                                have the same weighting factor
        :param input_prefix: prefix of the input files
        :return:
        """
        # Determine the number of runs
        run_number = len(particle_lat) // settings.INPUT_DIV + 1

        # Get all the variables that need to be split into indivual run files
        if particle_weight is not None:
            var_dict = {"lon": particle_lon, "lat": particle_lat, "weight": particle_weight}
        else:
            var_dict = {"lon": particle_lon, "lat": particle_lat}

        # Split the runs
        for run in range(run_number):
            for variables in var_dict:
                var_run = var_dict[variables][run * settings.INPUT_DIV: (run + 1) * settings.INPUT_DIV]
                np.save(input_prefix + "{}_run={}.npy".format(variables, run), var_run)

        return int(run_number)
