import settings
import utils
import glob
from netCDF4 import Dataset
import os

import advection_scenarios.create_boundary_current as create_boundary_current
import advection_scenarios.create_distance_to_shore as create_distance_to_shore
import advection_scenarios.create_land_ID as create_land_ID
import advection_scenarios.create_grid_spacing as create_grid_spacing
import advection_scenarios.create_input_files as create_input_files
import advection_scenarios.create_tidal_Kz_files as create_tidal_Kz_files


class AdvectionFiles:
    """
    We create a class, of which the variables are the names of the files we need to load. Therefore, the fieldset
    factory can just call upon this class to see which files it needs to read

    We also check all the relevant files exist, and if it doesn't, we create it,
    """

    def __init__(self, advection_scenario=settings.ADVECTION_DATA, repeat_dt=None):
        self.data_dir = settings.DATA_DIREC
        self.input_dir = settings.DATA_INPUT_DIREC
        self.advection_scenario = advection_scenario
        self.repeat_dt = repeat_dt
        os.system('echo "The advection scenario is "' + self.advection_scenario)

    @property
    def file_names(self):
        file_dict = {}
        if self.advection_scenario == 'HYCOM_GLOBAL':
            prefix = 'HYCOM_GLOBAL'
            # The core UV velocity fields
            UV_filenames = glob.glob(self.data_dir + "HYCOM/HYCOM_Surface*2000-01-01.nc") + \
                           glob.glob(self.data_dir + "HYCOM/HYCOM_Surface*{}*.nc".format(
                               settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            UV_filenames = list(dict.fromkeys(UV_filenames))
            UV_filenames.sort()
            UV_variables = {'U': 'water_u', 'V': 'water_v'}
            UV_dimensions = {'U': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'},
                             'V': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'}}
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='UV_filenames', variable=UV_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='UV_variables', variable=UV_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='UV_dimensions',
                                         variable=UV_dimensions)

            # The stokes drift fields
            STOKES_filenames = glob.glob(self.data_dir + "WaveWatchIIIstokes/ww3.{}*_uss.nc".format(
                settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            STOKES_filenames.sort()
            STOKES_variables = {'Ust': 'uuss', 'Vst': 'vuss'}
            STOKES_dimensions = {'Ust': {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'},
                                 'Vst': {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'},
                                 }
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='STOKES_filenames',
                                         variable=STOKES_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='STOKES_variables',
                                         variable=STOKES_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='STOKES_dimensions',
                                         variable=STOKES_dimensions)

            # The sea surface elevation fields
            ELEV_filenames = glob.glob(
                self.data_dir + "HYCOM/HYCOM_SeaEleve_3h_{}*.nc".format(settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            ELEV_filenames.sort()
            ELEV_variables = {'eta': 'surf_el'}
            ELEV_dimensions = {'time': 'time', 'lat': 'lat', 'lon': 'lon'}
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='ELEV_filenames',
                                         variable=ELEV_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='ELEV_variables',
                                         variable=ELEV_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='ELEV_dimensions',
                                         variable=ELEV_dimensions)

            # The surface winds
            WIND_filenames = glob.glob(self.data_dir + "Wind/ERA5*{}.nc".format(settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            WIND_filenames.sort()
            WIND_variables = {'u10': 'u10', 'v10': 'v10'}
            WIND_dimensions = {'time': 'time', 'lat': 'latitude', 'lon': 'longitude'}
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='WIND_filenames',
                                         variable=WIND_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='WIND_variables',
                                         variable=WIND_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='WIND_dimensions',
                                         variable=WIND_dimensions)

            # The shore type
            COAST_TYPE_filename = self.input_dir + 'HYCOM_GLOBAL_Luijendijk_sandy.npy'
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='COAST_TYPE_filename',
                                         variable=COAST_TYPE_filename)

            # Basic grid data
            dataset = Dataset(UV_filenames[0])
            GRID = dataset.variables['water_u'][0, 0, :, :]
            LON = dataset.variables['lon'][:]
            LAT = dataset.variables['lat'][:]
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='GRID', variable=GRID)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='LON', variable=LON)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='LAT', variable=LAT)

        if self.advection_scenario == 'HYCOM_CARIBBEAN':
            prefix = 'HYCOM_CARIBBEAN'
            # The core UV velocity fields
            UV_filenames = glob.glob(self.data_dir + "HYCOM_CARIBBEAN/hycom_gomu_501_2010010100_t000.nc.nc4") + \
                           glob.glob(self.data_dir + "HYCOM_CARIBBEAN/hycom_gomu_501_{}*.nc.nc4".format(
                               settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            # Remove duplicates
            UV_filenames = list(dict.fromkeys(UV_filenames))
            UV_filenames.sort()
            UV_variables = {'U': 'water_u', 'V': 'water_v'}
            UV_dimensions = {'U': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'},
                             'V': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'},
                             }
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='UV_filenames', variable=UV_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='UV_variables', variable=UV_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='UV_dimensions',
                                         variable=UV_dimensions)

            # The stokes drift fields
            STOKES_filenames = glob.glob(self.data_dir + "WaveWatchIIIstokes/ww3.{}*_uss.nc".format(
                settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            STOKES_filenames.sort()
            STOKES_variables = {'Ust': 'uuss', 'Vst': 'vuss'}
            STOKES_dimensions = {'Ust': {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'},
                                 'Vst': {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'},
                                 }
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='STOKES_filenames',
                                         variable=STOKES_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='STOKES_variables',
                                         variable=STOKES_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='STOKES_dimensions',
                                         variable=STOKES_dimensions)

            # The sea surface elevation fields
            ELEV_filenames = glob.glob(self.data_dir + "HYCOM_CARIBBEAN/hycom_gomu_501_{}*.nc.nc4".format(
                settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            ELEV_filenames.sort()
            ELEV_variables = {'eta': 'surf_el'}
            ELEV_dimensions = {'time': 'time', 'lat': 'lat', 'lon': 'lon'}
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='ELEV_filenames',
                                         variable=ELEV_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='ELEV_variables',
                                         variable=ELEV_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='ELEV_dimensions',
                                         variable=ELEV_dimensions)

            # The surface winds
            WIND_filenames = glob.glob(self.data_dir + "Wind/ERA5*{}.nc".format(settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            WIND_filenames.sort()
            WIND_variables = {'u10': 'u10', 'v10': 'v10'}
            WIND_dimensions = {'time': 'time', 'lat': 'latitude', 'lon': 'longitude'}
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='WIND_filenames',
                                         variable=WIND_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='WIND_variables',
                                         variable=WIND_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='WIND_dimensions',
                                         variable=WIND_dimensions)

            # The shore type
            COAST_TYPE_filename = self.input_dir + 'HYCOM_GLOBAL_Luijendijk_sandy.npy'
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='COAST_TYPE_filename',
                                         variable=COAST_TYPE_filename)

            # Basic grid data
            dataset = Dataset(UV_filenames[0])
            GRID = dataset.variables['water_u'][0, 0, :, :]
            LON = dataset.variables['lon'][:]
            LAT = dataset.variables['lat'][:]
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='GRID', variable=GRID)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='LON', variable=LON)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='LAT', variable=LAT)

        if self.advection_scenario == 'CMEMS_MEDITERRANEAN':
            prefix = 'CMEMS_MEDITERRANEAN'
            # The core UV velocity fields
            UV_filenames = glob.glob(self.data_dir + "CMEMS_MED/20100101_d-CMCC--RFVL-MFSe3r1-MED-b20200901_re-sv01.00.nc") + \
                           glob.glob(self.data_dir + "CMEMS_MED/{}*--RFVL*.nc".format(
                               settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            # Remove duplicates
            UV_filenames = list(dict.fromkeys(UV_filenames))
            UV_filenames.sort()
            UV_variables = {'U': 'uo', 'V': 'vo'}
            UV_dimensions = {'U': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'},
                             'V': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'}}
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='UV_filenames', variable=UV_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='UV_variables', variable=UV_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='UV_dimensions',
                                         variable=UV_dimensions)

            # The stokes drift fields
            STOKES_filenames = glob.glob(self.data_dir + "WAVE_MED/{}*--WAVE*.nc".format(
                settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            STOKES_filenames.sort()
            STOKES_variables = {'Ust': 'VSDX', 'Vst': 'VSDY'}
            STOKES_dimensions = {'Ust': {'time': 'time', 'lat': 'lat', 'lon': 'lon'},
                                 'Vst': {'time': 'time', 'lat': 'lat', 'lon': 'lon'},
                                 }
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='STOKES_filenames',
                                         variable=STOKES_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='STOKES_variables',
                                         variable=STOKES_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='STOKES_dimensions',
                                         variable=STOKES_dimensions)

            # The peak wave period fields
            PERIOD_filenames = glob.glob(self.data_dir + "WAVE_MED/{}*--WAVE*.nc".format(
                settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            PERIOD_filenames.sort()
            PERIOD_variables = {'WP': 'VTPK'}
            PERIOD_dimensions = {'WP': {'time': 'time', 'lat': 'lat', 'lon': 'lon'},
                                 }
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='PERIOD_filenames',
                                         variable=PERIOD_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='PERIOD_variables',
                                         variable=PERIOD_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='PERIOD_dimensions',
                                         variable=PERIOD_dimensions)

            # Ocean Temperature
            TEMP_filenames = glob.glob(self.data_dir + "CMEMS_MED/{}*--TEMP*.nc".format(
                settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            TEMP_filenames.sort()
            TEMP_variables = {'cons_temperature': 'thetao'}
            TEMP_dimensions = {'cons_temperature': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'}}
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='TEMP_filenames',
                                         variable=TEMP_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='TEMP_variables',
                                         variable=TEMP_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='TEMP_dimensions',
                                         variable=TEMP_dimensions)

            # Ocean Salinity
            SALINITY_filenames = glob.glob(self.data_dir + "CMEMS_MED/{}*--PSAL*.nc".format(
                settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            SALINITY_filenames.sort()
            SALINITY_variables = {'abs_salinity': 'so'}
            SALINITY_dimensions = {'abs_salinity': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'}}
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='SALINITY_filenames',
                                         variable=SALINITY_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='SALINITY_variables',
                                         variable=SALINITY_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='SALINITY_dimensions',
                                         variable=SALINITY_dimensions)

            # Ocean bathymetry
            BATH_filenames = self.data_dir + "CMEMS_MED/MED-MFC_006_004_mask_bathy.nc"
            BATH_variables = {'DEPTH': 'deptho', 'LON': 'longitude', 'LAT': 'latitude'}
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='BATH_filenames',
                                         variable=BATH_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='BATH_variables',
                                         variable=BATH_variables)

            # The surface winds
            WIND_filenames = glob.glob(self.data_dir + "Wind/ERA5*{}.nc".format(settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            WIND_filenames.sort()
            WIND_variables = {'u10': 'u10', 'v10': 'v10'}
            WIND_dimensions = {'time': 'time', 'lat': 'latitude', 'lon': 'longitude'}
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='WIND_filenames',
                                         variable=WIND_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='WIND_variables',
                                         variable=WIND_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='WIND_dimensions',
                                         variable=WIND_dimensions)

            # Basic grid data
            dataset = Dataset(UV_filenames[0])
            GRID = dataset.variables[UV_variables['U']][0, 0, :, :]
            LON = dataset.variables['lon'][:]
            LAT = dataset.variables['lat'][:]
            DEPTH = dataset.variables['depth'][:]
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='GRID', variable=GRID)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='LON', variable=LON)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='LAT', variable=LAT)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='DEPTH', variable=DEPTH)

            # The Mixed Layer Depth (MLD)
            MLD_filenames = glob.glob(self.data_dir + "CMEMS_MED/{}*--AMXL*.nc".format(
                settings.STARTYEAR + settings.RESTART * settings.BACKWARD_MULT))
            MLD_filenames.sort()
            MLD_variables = {'MLD': 'mlotst'}
            MLD_dimensions = {'MLD': {'time': 'time', 'lat': 'lat', 'lon': 'lon'}}
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='MLD_filenames', variable=MLD_filenames)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='MLD_variables', variable=MLD_variables)
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='MLD_dimensions', variable=MLD_dimensions)

            # Vertical Kz due to internal tides (TKZ)
            TIDE_Kz_filenames = self.input_dir + prefix + "_Kz_TIDAL.nc"
            file_dict = add_to_file_dict(file_dict=file_dict, variable_name='TIDE_Kz_filenames', variable=TIDE_Kz_filenames)
            if not utils.check_file_exist(TIDE_Kz_filenames):
                utils.print_statement("The tidal Kz file does not yet exist")
                create_tidal_Kz_files.create_tidal_Kz_files(file_name=TIDE_Kz_filenames, LON=LON, LAT=LAT, DEPTH=DEPTH,
                                                            BATH_filenames=BATH_filenames,
                                                            BATH_variables=BATH_variables)
            else:
                utils.print_statement("The tidal Kz file already exists")

        # The border current
        BORDER_filename = self.input_dir + prefix + '_boundary_velocities.nc'
        file_dict = add_to_file_dict(file_dict=file_dict, variable_name='BORDER_filename',
                                     variable=BORDER_filename)
        if utils.check_file_exist(BORDER_filename):
            utils.print_statement("border current file exists")
        else:
            utils.print_statement("The border current file does not yet exist")
            create_boundary_current.create_border_current(output_name=BORDER_filename, filenames=UV_filenames,
                                                          variables=UV_variables, dimensions=UV_dimensions,
                                                          grid=GRID)
            if utils.check_file_exist(BORDER_filename):
                utils.print_statement("The border current for {} has been created".format(prefix))
            else:
                utils.print_statement("The border current file still does not exist")

        # Distance to shore (from ocean cells to the nearest land cell)
        DISTANCE_filename = self.input_dir + prefix + '_distance2coast.nc'
        file_dict = add_to_file_dict(file_dict=file_dict, variable_name='DISTANCE_filename',
                                     variable=DISTANCE_filename)
        if utils.check_file_exist(DISTANCE_filename):
            utils.print_statement("distance to coast file exists")
        else:
            utils.print_statement("The distance to coast file does not yet exist")
            create_distance_to_shore.create_distance_to_shore(output_name=DISTANCE_filename, grid=GRID, lon=LON,
                                                              lat=LAT)
            if utils.check_file_exist(BORDER_filename):
                utils.print_statement("The border current for {} has been created".format(prefix))
            else:
                utils.print_statement("The border current file still does not exist")

        # Land ID
        LANDID_filename = self.input_dir + prefix + '_land_cell_identifier.nc'
        file_dict = add_to_file_dict(file_dict=file_dict, variable_name='LANDID_filename',
                                     variable=LANDID_filename)
        if utils.check_file_exist(LANDID_filename):
            utils.print_statement("LANDID file exists")
        else:
            utils.print_statement("The LANDID file does not yet exist")
            create_land_ID.create_land_ID(output_name=LANDID_filename, grid=GRID, lon=LON, lat=LAT)
            if utils.check_file_exist(BORDER_filename):
                utils.print_statement("The LANDID file for {} has been created".format(prefix))
            else:
                utils.print_statement("The LANDID file still does not exist")

        # Grid size
        GRIDSPACING_filename = self.input_dir + prefix + '_grid_spacing.nc'
        file_dict = add_to_file_dict(file_dict=file_dict, variable_name='GRIDSPACING_filename',
                                     variable=GRIDSPACING_filename)
        if utils.check_file_exist(GRIDSPACING_filename):
            utils.print_statement("Grid spacing file exists")
        else:
            utils.print_statement("The grid spacing file does not yet exist")
            create_grid_spacing.create_grid_spacing(output_name=GRIDSPACING_filename, grid=GRID, lon=LON, lat=LAT)
            if utils.check_file_exist(GRIDSPACING_filename):
                utils.print_statement("The grid spacing file for {} has been created".format(prefix))
            else:
                utils.print_statement("The grid spacing file still does not exist")

        # Checking for the input files
        output_prefix = create_input_files.create_input_files(advection_prefix=prefix, grid=GRID, lon=LON, lat=LAT,
                                                              repeat_dt=self.repeat_dt).create_files()
        STARTFILES_filename = {}
        for variable in ['lon', 'lat', 'weight']:
            file_name = output_prefix + '{}_run={}.npy'.format(variable, settings.RUN)
            if utils.check_file_exist(file_name):
                STARTFILES_filename[variable] = file_name
        file_dict = add_to_file_dict(file_dict=file_dict, variable_name='STARTFILES_filename',
                                     variable=STARTFILES_filename)

        return file_dict


def add_to_file_dict(file_dict, variable_name: str, variable):
    file_dict[variable_name] = variable
    return file_dict

