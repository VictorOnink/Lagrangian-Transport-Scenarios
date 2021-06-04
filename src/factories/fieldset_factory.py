import os
import settings as settings
from netCDF4 import Dataset
import numpy as np
from parcels import GeographicPolar, Geographic, FieldSet, Field
import math
import utils


class FieldSetFactory():
    """"""

    @classmethod
    def create_fieldset(cls, file_dict: dict, stokes: int,
                        stokes_depth: bool = False,
                        border_current: bool = False,
                        diffusion: bool = False,
                        landID: bool = False,
                        distance: bool = False,
                        wind: bool = False,
                        sea_elev: bool = False,
                        salinity: bool = False,
                        temperature: bool = False,
                        bathymetry: bool = False,
                        vicinity: bool = False,
                        beach_timescale: bool = False,
                        resus_timescale: bool = False,
                        wind_min: bool = False,
                        MLD: bool = False,
                        KPP_mixing: bool = False,
                        TIDAL_mixing: bool = False,
                        seabed_resuspension: bool = False,
                        coastal_zone: bool = True,
                        grid_spacing: bool = True,
                        halo: bool = True
                        ):
        fieldset = _get_base_fieldset(file_dict=file_dict)
        if stokes is 0:
            fieldset = _add_stokes_drift(fieldset=fieldset, file_dict=file_dict)
        if stokes_depth:
            _add_stokes_depth_depen(fieldset=fieldset, file_dict=file_dict)
        if border_current:
            _add_border_current(fieldset=fieldset, file_dict=file_dict)
        if diffusion:
            _add_diffusion(fieldset=fieldset, file_dict=file_dict)
        if landID:
            _add_land_ID_field(fieldset=fieldset, file_dict=file_dict)
        if distance:
            _add_distance2shore_field(fieldset=fieldset, file_dict=file_dict)
        if wind:
            _add_wind_field(fieldset=fieldset, file_dict=file_dict)
        if sea_elev:
            _add_sea_elevation_field(fieldset=fieldset, file_dict=file_dict)
        if salinity:
            _add_salinity_field(fieldset=fieldset, file_dict=file_dict)
        if temperature:
            _add_temperature_field(fieldset=fieldset, file_dict=file_dict)
        if bathymetry:
            _add_bathymetry_field(fieldset=fieldset, file_dict=file_dict)
        if vicinity:
            _add_vicinity_constant(fieldset=fieldset)
        if beach_timescale:
            _add_beach_timescale_field(fieldset=fieldset)
        if resus_timescale:
            _add_resus_timescale_field(fieldset=fieldset, file_dict=file_dict)
        if wind_min:
            _add_min_resuspension_wind_constant(fieldset=fieldset)
        if MLD:
            _add_MLD_field(fieldset=fieldset, file_dict=file_dict)
        if KPP_mixing:
            _add_KPP_wind_mixing(fieldset=fieldset, file_dict=file_dict)
        if TIDAL_mixing:
            _add_TIDAL_mixing(fieldset=fieldset, file_dict=file_dict)
        if seabed_resuspension:
            _add_seabed_resuspension(fieldset=fieldset)
        if coastal_zone:
            _add_coastal_zone_boundary(fieldset=fieldset)
        if grid_spacing:
            _add_grid_spacing(fieldset=fieldset, file_dict=file_dict)
        if halo:
            _add_halo(fieldset)
        return fieldset


def _get_base_fieldset(file_dict: dict) -> FieldSet:
    """

    :param data_dir:
    :return:
    """
    # Defining the folders in which all the data is stored on the different servers
    os.system('echo "Creating the main fieldset"')
    # Loading in the surface currents, where we always load in the 2000-01-01 file to ensure that time is always given
    # relative to the same starting point, whereas we then only load in the files for the specific year that the
    # simulation runs in. This speeds up the fieldset creation somewhat.
    _check_presence(variable='UV_filenames', file_dict=file_dict)
    filenames = {'U': file_dict['UV_filenames'],
                 'V': file_dict['UV_filenames'],
                 }
    fieldset = FieldSet.from_netcdf(filenames, file_dict['UV_variables'], file_dict['UV_dimensions'],
                                    allow_time_extrapolation=True)
    return fieldset


def _add_stokes_drift(fieldset: FieldSet, file_dict: dict):
    """
    Do we include stokes drift yes or no.
    :param fieldset:
    :param input_dir:
    """
    os.system('echo "Adding Stokes drift"')
    _check_presence(variable='STOKES_filenames', file_dict=file_dict)
    filenames = {'Ust': file_dict['STOKES_filenames'], 'Vst': file_dict['STOKES_filenames']}
    fieldset_stoke = FieldSet.from_netcdf(filenames, file_dict['STOKES_variables'], file_dict['STOKES_dimensions'],
                                          allow_time_extrapolation=True)

    fieldset_stoke.Ust.units = GeographicPolar()
    fieldset_stoke.Vst.units = Geographic()
    # Adding the Stokes drift fields to the general fieldset
    fieldset.add_field(fieldset_stoke.Ust)
    fieldset.add_field(fieldset_stoke.Vst)
    return fieldset


def _add_stokes_depth_depen(fieldset: FieldSet, file_dict: dict):
    os.system('echo "Adding Stokes drift depth dependence data"')
    _check_presence(variable='PERIOD_filenames', file_dict=file_dict)
    filenames = {'WP': file_dict['PERIOD_filenames']}
    fieldset_period = FieldSet.from_netcdf(filenames, file_dict['PERIOD_variables'], file_dict['PERIOD_dimensions'],
                                           allow_time_extrapolation=True)
    fieldset.add_field(fieldset_period.WP)


def _add_border_current(fieldset: FieldSet, file_dict: dict):
    """

    :param fieldset:
    :param input_dir:
    """
    os.system('echo "Adding the border current"')
    _check_presence(variable='BORDER_filename', file_dict=file_dict)
    datasetBor = Dataset(file_dict['BORDER_filename'])
    borU = datasetBor.variables['border_u'][:]
    borV = datasetBor.variables['border_v'][:]
    # # Adding the actual field
    fieldset.add_field(Field('borU', borU, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))
    fieldset.add_field(Field('borV', borV, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))
    # making sure the units are interpreted as m s^-1
    fieldset.borU.units = GeographicPolar()
    fieldset.borV.units = Geographic()


def _add_diffusion(fieldset: FieldSet, file_dict: dict):
    """

    :param fieldset:
    :param input_dir:
    """
    os.system('echo "Adding diffusion"')
    kh = settings.K_HOR  # m^2 s^-1
    _check_presence(variable='GRID', file_dict=file_dict)
    mask = file_dict['GRID'].mask
    kh_f = kh * np.ones(mask.shape)
    kh_f[mask == True] = 0
    fieldset.add_field(Field('Kh_zonal', kh_f, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))
    fieldset.add_field(Field('Kh_meridional', kh_f, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))
    # In case we use the Euler-Maruyama scheme for diffusion, we need a parameter `fieldset.dres`, to set the resolution
    # for the central difference gradient approximation. The value I've taken from the parcels diffusion tutorial, and
    # should be sufficiently smaller than the spatial resolution of any of the data I work with
    fieldset.add_constant('dres', 0.00005)


def _add_land_ID_field(fieldset: FieldSet, file_dict: dict):
    """

    :param fieldset:
    :param input_dir:
    """
    os.system('echo "Adding land/water boolean field"')
    _check_presence(variable='LANDID_filename', file_dict=file_dict)
    dataset = Dataset(file_dict['LANDID_filename'])
    landID = dataset.variables['land_ID'][:]
    fieldset.add_field(Field('landID', landID, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))


def _add_distance2shore_field(fieldset: FieldSet, file_dict: dict):
    """

    :param fieldset:
    :param input_dir:
    :return:
    """
    os.system('echo "Adding distance to shore"')
    _check_presence(variable='DISTANCE_filename', file_dict=file_dict)
    datasetCoast = Dataset(file_dict['DISTANCE_filename'])
    distance = datasetCoast.variables['distance'][:]
    fieldset.add_field(Field('distance2shore', distance, lon=file_dict['LON'], lat=file_dict['LAT'],
                             mesh='spherical'))


def _add_wind_field(fieldset: FieldSet, file_dict: dict):
    """

    :param fieldset:
    :param input_dir:
    :return:
    """
    os.system('echo "Adding 10m winds"')
    _check_presence(variable='WIND_filenames', file_dict=file_dict)
    filenames = {'u10': file_dict['WIND_filenames'],
                 'v10': file_dict['WIND_filenames']}
    # Creating a fieldset for the wind data
    fieldset_wind = FieldSet.from_netcdf(filenames, file_dict['WIND_variables'], file_dict['WIND_dimensions'],
                                         allow_time_extrapolation=True)
    # Adding the wind fields to the general fieldset
    fieldset.add_field(fieldset_wind.u10)
    fieldset.add_field(fieldset_wind.v10)


def _add_MLD_field(fieldset: FieldSet, file_dict: dict):
    os.system('echo "Adding Mixed Layer Depth"')
    _check_presence(variable='MLD_filenames', file_dict=file_dict)
    filenames = {'MLD': file_dict['MLD_filenames']}
    fieldset_MLD = FieldSet.from_netcdf(filenames, file_dict['MLD_variables'], file_dict['MLD_dimensions'],
                                        allow_time_extrapolation=True)
    # Adding the MLD field to the general fieldset
    fieldset.add_field(fieldset_MLD.MLD)


def _add_sea_elevation_field(fieldset: FieldSet, file_dict: dict):
    os.system('echo "Adding sea surface elevation"')
    _check_presence(variable='ELEV_filenames', file_dict=file_dict)
    filenames = {'eta': file_dict['ELEV_filenames']}
    # Creating a fieldset for the wind data
    fieldset_sea = FieldSet.from_netcdf(filenames, file_dict['ELEV_variables'], file_dict['ELEV_dimensions'],
                                        allow_time_extrapolation=True)
    # Adding the wind fields to the general fieldset
    fieldset.add_field(fieldset_sea.eta)


def _add_salinity_field(fieldset: FieldSet, file_dict: dict):
    os.system('echo "Adding ocean salinity"')
    _check_presence(variable='SALINITY_filenames', file_dict=file_dict)
    filenames = {'abs_salinity': file_dict['SALINITY_filenames']}
    # Creating a fieldset for the salinity data
    fieldset_sal = FieldSet.from_netcdf(filenames, file_dict['SALINITY_variables'], file_dict['SALINITY_dimensions'],
                                        allow_time_extrapolation=True)
    # Adding the wind fields to the general fieldset
    fieldset.add_field(fieldset_sal.abs_salinity)


def _add_temperature_field(fieldset: FieldSet, file_dict: dict):
    os.system('echo "Adding ocean temperature"')
    _check_presence(variable='TEMP_filenames', file_dict=file_dict)
    filenames = {'cons_temperature': file_dict['TEMP_filenames']}
    # Creating a fieldset for the temperature data
    fieldset_temp = FieldSet.from_netcdf(filenames, file_dict['TEMP_variables'], file_dict['TEMP_dimensions'],
                                         allow_time_extrapolation=True)
    # Adding the temperature field to the general fieldset
    fieldset.add_field(fieldset_temp.cons_temperature)


def _add_bathymetry_field(fieldset: FieldSet, file_dict: dict):
    os.system('echo "Adding ocean bathymetry"')
    _check_presence(variable='BATH_filenames', file_dict=file_dict)
    dataset = Dataset(file_dict['BATH_filenames'])
    bathymetry = np.array(dataset.variables[file_dict['BATH_variables']['DEPTH']][:])
    bathymetry[dataset.variables[file_dict['BATH_variables']['DEPTH']][:].mask] = np.nanmin(np.abs(file_dict['DEPTH']))
    fieldset.add_field(Field('bathymetry', bathymetry, lon=dataset.variables[file_dict['BATH_variables']['LON']][:],
                             lat=dataset.variables[file_dict['BATH_variables']['LAT']][:], mesh='spherical'))


def _add_vicinity_constant(fieldset: FieldSet):
    """

    :param fieldset:
    :return:
    """
    # The vicinity timescale
    fieldset.add_constant('vic', settings.VICINITY)


def _add_beach_timescale_field(fieldset: FieldSet):
    # Here only the beaching probability is a global constant, the resuspension
    # probability will instead be represented using a field
    p_b = math.exp(-settings.TIME_STEP.total_seconds() / (settings.SHORE_TIME * 86400.))
    fieldset.add_constant('p_beach', p_b)


def _compute_shore_resus_Field(file_dict: dict):
    """

    :param input_dir:
    :return:
    """
    if settings.SCENARIO_NAME == 'ShoreDependentResuspension':
        _check_presence(variable='COAST_TYPE_filename', file_dict=file_dict)
        s = np.load(file_dict['COAST_TYPE_filename'])
        if settings.SHORE_DEP == 0:
            resus_field = settings.RESUS_TIME * (0.75 + 0.25 * s)
        if settings.SHORE_DEP == 1:
            resus_field = settings.RESUS_TIME * (0.25 + 0.75 * s)
        return resus_field


def _add_resus_timescale_field(fieldset: FieldSet, file_dict: dict):
    """

    :param fieldset:
    :param input_dir:
    :return:
    """
    if settings.SCENARIO_NAME == 'ShoreDependentResuspension':
        p_r = np.exp(-settings.TIME_STEP.total_seconds() / (_compute_shore_resus_Field(file_dict) * 86400.))
        fieldset.add_field(Field('p_resus', p_r, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))
    if settings.SCENARIO_NAME == 'SizeTransport':
        lambda_R = utils.get_resuspension_timescale()
        p_r = math.exp(-settings.TIME_STEP.total_seconds() / (lambda_R * 86400.))
        fieldset.add_constant('p_resus', p_r)
    else:
        p_r = math.exp(-settings.TIME_STEP.total_seconds() / (settings.RESUS_TIME * 86400.))
        fieldset.add_constant('p_resus', p_r)


def _add_min_resuspension_wind_constant(fieldset: FieldSet):
    # The minimum offshore wind speed for resuspension to be possible
    fieldset.add_constant('Wmin', settings.WMIN / 10)


def _add_coastal_zone_boundary(fieldset: FieldSet):
    # The distance from the nearest coastline that is defines the coastal zone within which beaching occurs
    fieldset.add_constant('Coastal_Boundary', settings.COAST_D)


def _add_KPP_wind_mixing(fieldset: FieldSet, file_dict: dict):
    """
    Adding a number of constants that are required for computing the wind-driven mixing according to the KPP scheme
    """
    os.system('echo "Adding KPP wind mixing parameters"')
    fieldset.add_constant('G', settings.G)  # Acceleration due to gravity
    fieldset.add_constant('VK', settings.VK)  # von Karman coefficient
    fieldset.add_constant('RHO_A', settings.RHO_A)  # Density of air
    fieldset.add_constant('PHI', settings.PHI)  # Stability function
    fieldset.add_constant('BETA', settings.BETA)  # Wave age based on 10m wind speed
    fieldset.add_constant('BETA_STAR', settings.BETA_STAR)  # Wave age based on frictional air velocity
    fieldset.add_constant('SURF_Z',
                          np.nanmin(np.abs(file_dict['DEPTH'])))  # Correction in case the surface depth is not z=0


def _add_TIDAL_mixing(fieldset: FieldSet, file_dict: dict):
    os.system('echo "Adding TIDAL Kz mixing"')
    _check_presence(variable='TIDE_Kz_filenames', file_dict=file_dict)
    dataset = Dataset(file_dict['TIDE_Kz_filenames'])
    TIDAL_Kz = dataset.variables['TIDAL_Kz'][:]
    TIDAL_dKz = dataset.variables['TIDAL_dKz'][:]
    fieldset.add_field(Field('TIDAL_Kz', TIDAL_Kz, lon=file_dict['LON'], lat=file_dict['LAT'], depth=file_dict['DEPTH'],
                             mesh='spherical'))
    fieldset.add_field(Field('TIDAL_dKz', TIDAL_dKz, lon=file_dict['LON'], lat=file_dict['LAT'],
                             depth=file_dict['DEPTH'], mesh='spherical'))
    fieldset.add_constant('K_Z_BULK', settings.K_Z_BULK)


def _add_grid_spacing(fieldset: FieldSet, file_dict: dict):
    # Adding in the lon and lat grid spacing for use in the initial scattering of particles on the first time step
    os.system('echo "Adding lon and lat grid spacing"')
    _check_presence(variable='GRIDSPACING_filename', file_dict=file_dict)
    datasetCoast = Dataset(file_dict['GRIDSPACING_filename'])
    dlon = datasetCoast.variables['lon_spacing'][:]
    dlat = datasetCoast.variables['lat_spacing'][:]
    fieldset.add_field(Field('dlon', dlon, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))
    fieldset.add_field(Field('dlat', dlat, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))


def _add_halo(fieldset: FieldSet):
    os.system('echo "Finally, the periodic halo"')
    fieldset.add_periodic_halo(zonal=True)


def _check_presence(variable: str, file_dict: dict):
    if variable not in file_dict.keys():
        str_format = (variable, variable)
        os.system('echo "The variable {} is not within file_dict. Check that {} is in the advection scenario"'.format(
            *str_format))


def _add_seabed_resuspension(fieldset: FieldSet):
    # Adding in the constant critical shear stress for particle resuspension from the sea bed
    os.system('echo "Adding seabed resuspension constant"')
    fieldset.add_constant('SEABED_CRIT', settings.SEABED_CRIT)