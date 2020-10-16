import settings
import utils

import glob
from netCDF4 import Dataset
import os

import advection_scenarios.create_boundary_current as create_boundary_current

class AdvectionFiles:
    """
    We create a class, of which the variables are the names of the files we need to load. Therefore, the fieldset
    factory can just call upon this class to see which files it needs to read

    We also check all the relevant files exist, and if it doesn't, we create it,
    """

    def __init__(self, server, stokes, advection_scenario):
        self.server = server
        self.stokes = stokes
        self.data_dir = utils._get_data_directory(server=server)
        self.input_dir = utils._get_input_directory(server=server)
        self.advection_scenario = advection_scenario
        os.system('echo "The advection scenario is "' + self.advection_scenario)

    def file_names(self):
        if self.advection_scenario == 'HYCOM_GLOBAL':
            prefix = 'HYCOM_GLOBAL'
            # The core UV velocity fields
            UV_filenames = glob.glob(self.data_dir + "HYCOM/HYCOM_Surface*2000-01-01.nc") + \
                           glob.glob(self.data_dir + "HYCOM/HYCOM_Surface*{}*.nc".format(
                               settings.START_YEAR + settings.RESTART))
            # Remove duplicates
            UV_filenames = list(dict.fromkeys(UV_filenames))
            UV_filenames.sort()
            UV_variables = {'U': 'water_u', 'V': 'water_v'}
            UV_dimensions = {'U': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'},
                             'V': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'},
                             }

            # The stokes drift fields
            STOKES_filenames = glob.glob(self.data_dir + "WaveWatchIIIstokes/ww3.{}*_uss.nc".format(
                settings.START_YEAR + settings.RESTART))
            STOKES_filenames.sort()
            STOKES_variables = {'Ust': 'uuss', 'Vst': 'vuss'}
            STOKES_dimensions = {'Ust': {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'},
                                 'Vst': {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'},
                                 }

            # The sea surface elevation fields
            ELEV_filenames = glob.glob(
                self.data_dir + "HYCOM/HYCOM_SeaEleve_3h_{}*.nc".format(settings.START_YEAR + settings.RESTART))
            ELEV_filenames.sort()
            ELEV_variables = {'eta': 'surf_el'}
            ELEV_dimensions = {'time': 'time', 'lat': 'lat', 'lon': 'lon'}

            # The surface winds
            WIND_filenames = glob.glob(self.data_dir + "Wind/ERA5-wind10m*.nc")
            WIND_filenames.sort()
            WIND_variables = {'u10': 'u10', 'v10': 'v10'}
            WIND_dimensions = {'time': 'time', 'lat': 'latitude', 'lon': 'longitude'}

            # The shore type
            COAST_TYPE_filename = self.input_dir + 'HYCOM_GLOBAL_Luijendijk_sandy.npy'

            # Basic grid data
            dataset = Dataset(UV_filenames[0])
            GRID = dataset.variables['water_u'][0, 0, :, :]
            LON = dataset.variables['lon'][:]
            LAT = dataset.variables['lat'][:]

        if self.advection_scenario == 'HYCOM_CARIBBEAN':
            prefix = 'HYCOM_CARIBBEAN'
            # The core UV velocity fields
            UV_filenames = glob.glob(self.data_dir + "HYCOM_CARIBBEAN/hycom_gomu_501_2010010100_t000.nc.nc4") + \
                           glob.glob(self.data_dir + "HYCOM_CARIBBEAN/hycom_gomu_501_{}*.nc.nc4".format(
                               settings.START_YEAR + settings.RESTART))
            # Remove duplicates
            UV_filenames = list(dict.fromkeys(UV_filenames))
            UV_filenames.sort()
            UV_variables = {'U': 'water_u', 'V': 'water_v'}
            UV_dimensions = {'U': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'},
                             'V': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'},
                             }

            # The stokes drift fields
            STOKES_filenames = glob.glob(self.data_dir + "WaveWatchIIIstokes/ww3.{}*_uss.nc".format(
                settings.START_YEAR + settings.RESTART))
            STOKES_filenames.sort()
            STOKES_variables = {'Ust': 'uuss', 'Vst': 'vuss'}
            STOKES_dimensions = {'Ust': {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'},
                                 'Vst': {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'},
                                 }

            # The sea surface elevation fields
            ELEV_filenames = glob.glob(self.data_dir + "HYCOM_CARIBBEAN/hycom_gomu_501_{}*.nc.nc4".format(
                settings.START_YEAR + settings.RESTART))
            ELEV_filenames.sort()
            ELEV_variables = {'eta': 'surf_el'}
            ELEV_dimensions = {'time': 'time', 'lat': 'lat', 'lon': 'lon'}

            # The surface winds
            WIND_filenames = glob.glob(self.data_dir + "Wind/ERA5-wind10m*.nc")
            WIND_filenames.sort()
            WIND_variables = {'u10': 'u10', 'v10': 'v10'}
            WIND_dimensions = {'time': 'time', 'lat': 'latitude', 'lon': 'longitude'}

            # The shore type
            COAST_TYPE_filename = self.input_dir + 'HYCOM_GLOBAL_Luijendijk_sandy.npy'

            # Basic grid data
            dataset = Dataset(UV_filenames[0])
            GRID = dataset.variables['water_u'][0, 0, :, :]
            LON = dataset.variables['lon'][:]
            LAT = dataset.variables['lat'][:]

        # The border current
        BORDER_filename = self.input_dir + prefix + '_boundary_velocities.nc'
        if utils._check_file_exist(BORDER_filename):
            os.system('echo "border current file exists"')
        else:
            os.system('echo "The border current file does not yet exist"')
            create_boundary_current.create_border_current(output_name=BORDER_filename, filenames=UV_filenames,
                                                          variables=UV_variables, dimensions=UV_dimensions,
                                                          grid=GRID)
            if utils._check_file_exist(BORDER_filename):
                os.system('echo "The border current for "'+str(prefix)+' has been created')
            else:
                os.system('echo "The border current file still does not exist"')

        # Distance to shore (from ocean cells to the nearest land cell)
        DISTANCE_filename = self.input_dir + prefix + '_distance2coast.nc'
        if utils._check_file_exist(DISTANCE_filename):
            os.system('echo "distance to coast file exists"')
        else:
            os.system('echo "The distance to coast file does not yet exist"')

        # Land ID
        LANDID_filename = self.input_dir + prefix + '_land_cell_identifier.npy'
        if utils._check_file_exist(LANDID_filename):
            os.system('echo "LANDID file exists"')
        else:
            os.system('echo "The LANDID file does not yet exist"')

        # Putting everything in a dictionary so we can output it
        variable_names = ['UV_filenames', 'UV_variables', 'UV_dimensions', 'STOKES_filenames', 'STOKES_variables',
                          'STOKES_dimensions', 'ELEV_filenames', 'ELEV_variables', 'ELEV_dimensions',
                          'WIND_filenames', 'WIND_variables', 'WIND_dimensions', 'COAST_TYPE_filename', 'GRID', 'LON',
                          'LAT', 'BORDER_filename', 'DISTANCE_filename', 'LANDID_filename']
        variables = [UV_filenames, UV_variables, UV_dimensions, STOKES_filenames, STOKES_variables, STOKES_dimensions,
                     ELEV_filenames, ELEV_variables, ELEV_dimensions, WIND_filenames, WIND_variables, WIND_dimensions,
                     COAST_TYPE_filename, GRID, LON, LAT, BORDER_filename, DISTANCE_filename, LANDID_filename]

        file_dict = {}
        for var in range(len(variables)):
            file_dict[variable_names[var]] = variables[var]
        return file_dict
