import os
import settings as settings
from netCDF4 import Dataset
import numpy as np
from parcels import GeographicPolar, Geographic, FieldSet, Field
import math
from advection_scenarios import advection_files

class FieldSetFactory():
    """"""

    @classmethod
    def create_fieldset(cls, server: int, stokes: int,
                        border_current: bool = False,
                        diffusion: bool = False,
                        landID: bool = False,
                        distance: bool = False,
                        wind: bool = False,
                        sea_elev: bool = False,
                        vicinity: bool = False,
                        beach_timescale: bool = False,
                        resus_timescale: bool = False,
                        wind_min: bool = False,
                        coastal_zone: bool = True,
                        grid_spacing: bool = True,
                        halo: bool = True

                        ):
        """

        :param server:
        :param stokes:
        :param border_current:
        :param diffusion:
        :param landID:
        :param distance:
        :param wind:
        :param sea_elev:
        :param vicinity:
        :param beach_timescale:
        :param resus_timescale:
        :param wind_min:
        :param halo:
        :return:
        """
        advection_scenario = advection_files.AdvectionFiles(server=server, stokes=stokes,
                                                            advection_scenario=settings.ADVECTION_DATA)
        file_dict = advection_scenario.file_names()

        fieldset = _get_base_fieldset(file_dict=file_dict)
        if stokes == 0:
            fieldset = _add_stokes_drift(fieldset=fieldset, file_dict=file_dict)
        if border_current:
            _add_border_current(fieldset=fieldset, file_dict=file_dict)
        if diffusion:
            _add_diffusion(fieldset=fieldset, file_dict=file_dict)
        if landID:
            _add_landID_field(fieldset=fieldset, file_dict=file_dict)
        if distance:
            _add_distance2shore_field(fieldset=fieldset, file_dict=file_dict)
        if wind:
            _add_wind_field(fieldset=fieldset, file_dict=file_dict)
        if sea_elev:
            _add_seaElevation_field(fieldset=fieldset, file_dict=file_dict)
        if vicinity:
            _add_vicinity_constant(fieldset=fieldset)
        if beach_timescale:
            _add_beachTimescale_field(fieldset=fieldset)
        if resus_timescale:
            _add_resusTimescale_field(fieldset=fieldset, file_dict=file_dict)
        if wind_min:
            _add_minResuspensionWind_constant(fieldset=fieldset)
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
    filenames = {'Ust': file_dict['STOKES_filenames'], 'Vst': file_dict['STOKES_filenames']}
    fieldset_stoke = FieldSet.from_netcdf(filenames, file_dict['STOKES_variables'], file_dict['STOKES_dimensions'],
                                          allow_time_extrapolation=True)

    fieldset_stoke.Ust.units = GeographicPolar()
    fieldset_stoke.Vst.units = Geographic()
    # Adding the Stokes drift fields to the general fieldset
    fieldset = FieldSet(U=fieldset.U+fieldset_stoke.Ust,
                        V=fieldset.V+fieldset_stoke.Vst)
    return fieldset


def _add_border_current(fieldset: FieldSet, file_dict: dict):
    """

    :param fieldset:
    :param input_dir:
    """
    os.system('echo "Adding the border current"')
    datasetBor = Dataset(file_dict['BORDER_filename'])
    borU = datasetBor.variables['border_u'][:]
    borV = datasetBor.variables['border_v'][:]
    # Normalizing the border current so that the total current is always 1m/s
    # borMag = np.sqrt(np.square(borU) + np.square(borV))
    # borMag[borMag == 0] = 1
    # borU = np.divide(borU, borMag)
    # borV = np.divide(borV, borMag)
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
    kh = settings.KH_HOR  # m^2 s^-1
    mask = file_dict['GRID'].mask
    kh_f = kh * np.ones(mask.shape)
    kh_f[mask == True] = 0
    fieldset.add_field(Field('Kh_zonal', kh_f, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))
    fieldset.add_field(Field('Kh_meridional', kh_f, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))


def _add_landID_field(fieldset: FieldSet, file_dict: dict):
    """

    :param fieldset:
    :param input_dir:
    """
    os.system('echo "Adding land/water boolean field"')
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
    datasetCoast = Dataset(file_dict['DISTANCE_filename'])
    distance = datasetCoast.variables['distance'][:]
    fieldset.add_field(Field('distance2shore', distance, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))


def _add_wind_field(fieldset: FieldSet, file_dict: dict):
    """

    :param fieldset:
    :param input_dir:
    :return:
    """
    os.system('echo "Adding 10m winds"')
    filenames = {'u10': file_dict['WIND_filenames'],
                 'v10': file_dict['WIND_filenames']}
    # Creating a fieldset for the wind data
    fieldset_wind = FieldSet.from_netcdf(filenames, file_dict['WIND_variables'], file_dict['WIND_dimensions'],
                                         allow_time_extrapolation=True)
    fieldset_wind.u10.units = GeographicPolar()
    fieldset_wind.v10.units = Geographic()
    # Adding the wind fields to the general fieldset
    fieldset.add_field(fieldset_wind.u10)
    fieldset.add_field(fieldset_wind.v10)


def _add_seaElevation_field(fieldset: FieldSet, file_dict: dict):
    os.system('echo "Adding sea surface elevation"')
    filenames = {'eta': file_dict['ELEV_filenames']}
    # Creating a fieldset for the wind data
    fieldset_sea = FieldSet.from_netcdf(filenames, file_dict['ELEV_variables'], file_dict['ELEV_dimensions'],
                                        allow_time_extrapolation=True)
    # Adding the wind fields to the general fieldset
    fieldset.add_field(fieldset_sea.eta)


def _add_vicinity_constant(fieldset: FieldSet):
    """
    
    :param fieldset: 
    :return: 
    """
    # The vicinity timescale
    fieldset.add_constant('vic', settings.VICINITY)


def _add_beachTimescale_field(fieldset: FieldSet):
    # Here only the beaching probability is a global constant, the resuspension
    # probability will instead be represented using a field
    p_b = math.exp(-settings.TIME_STEP.total_seconds() / (settings.SHORE_TIME * 86400.))
    fieldset.add_constant('p_beach', p_b)


def _compute_ShoreResus_Field(file_dict: dict):
    """

    :param input_dir:
    :return:
    """
    if settings.SCENARIO_NAME == 'ShoreDependentResuspension':
        s = np.load(file_dict['COAST_TYPE_filename'])
        if settings.SHORE_DEP == 0:
            resus_field = settings.RESUS_TIME * (0.75 + 0.25 * s)
        if settings.SHORE_DEP == 1:
            resus_field = settings.RESUS_TIME * (0.25 + 0.75 * s)
        return resus_field


def _add_resusTimescale_field(fieldset: FieldSet, file_dict: dict):
    """

    :param fieldset:
    :param input_dir:
    :return:
    """
    if settings.SCENARIO_NAME == 'Stochastic':
        p_r = math.exp(-settings.TIME_STEP.total_seconds() / (settings.RESUS_TIME * 86400.))
        fieldset.add_constant('p_resus', p_r)

    elif settings.SCENARIO_NAME == 'ShoreDependentResuspension':
        p_r = np.exp(-settings.TIME_STEP.total_seconds() / (_compute_ShoreResus_Field(file_dict) * 86400.))
        fieldset.add_field(Field('p_resus', p_r, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))


def _add_minResuspensionWind_constant(fieldset: FieldSet):
    # The minimum offshore wind speed for resuspension to be possible
    fieldset.add_constant('Wmin', settings.WMIN / 10)


def _add_coastal_zone_boundary(fieldset: FieldSet):
    # The distance from the nearest coastline that is defines the coastal zone within which beaching occurs
    fieldset.add_constant('Coastal_Boundary', settings.COAST_D)

def _add_grid_spacing(fieldset: FieldSet, file_dict: dict):
    # Adding in the lon and lat grid spacing for use in the initial scattering of particles on the first time step
    os.system('echo "Adding lon and lat grid spacing"')
    datasetCoast = Dataset(file_dict['GRIDSPACING_filename'])
    dlon = datasetCoast.variables['lon_spacing'][:]
    dlat = datasetCoast.variables['lat_spacing'][:]
    fieldset.add_field(Field('dlon', dlon, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))
    fieldset.add_field(Field('dlat', dlat, lon=file_dict['LON'], lat=file_dict['LAT'], mesh='spherical'))


def _add_halo(fieldset: FieldSet):
    os.system('echo "Finally, the periodic halo"')
    fieldset.add_periodic_halo(zonal=True)
