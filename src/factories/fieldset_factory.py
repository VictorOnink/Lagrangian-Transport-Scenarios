import os
import src.settings as settings
from netCDF4 import Dataset
import numpy as np
from parcels import GeographicPolar, Geographic,FieldSet,Field
import math
from datetime import timedelta


class FieldSetFactory():
    """"""
    @classmethod
    def create_fieldset(cls, server: int, stokes: int,
                        stokes_drift: bool = False,
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
                        halo: bool = True
                        ):
        """

        :param server:
        :param stokes:
        :param stokes_drift:
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
        data_dir = _get_data_directory(server=server)
        input_dir = _get_input_directory(server=server)

        fieldset = _get_base_fieldset(data_dir=data_dir)
        if stokes_drift:
            _add_stokes_drift(fieldset=fieldset, input_dir=input_dir)
        if border_current:
            _add_border_current(fieldset=fieldset, input_dir=input_dir)
        if diffusion:
            _add_diffusion(fieldset=fieldset, input_dir=input_dir)
        if landID:
            _add_landID_field(fieldset=fieldset, input_dir=input_dir)
        if distance:
            _add_distance2shore_field(fieldset=fieldset, input_dir=input_dir)
        if wind:
            _add_wind_field(fieldset=fieldset,input_dir=input_dir)
        if sea_elev:
            _add_seaElevation_field(fieldset=fieldset,input_dir=input_dir)
        if vicinity:
            _add_vicinity_constant(fieldset=fieldset)
        if beach_timescale:
            _add_beachTimescale_field(fieldset=fieldset)
        if resus_timescale:
            _add_resusTimescale_field(fieldset=fieldset)
        if wind_min:
            _add_minResuspensionWind_constant(fieldset=fieldset)
        if halo:
            _add_halo(fieldset)
        return fieldset


def _get_data_directory(server: int) -> str:
    """

    :param server:
    :return:
    """
    return settings.DATA_DIR_SERVERS[server]


def _get_input_directory(server: int) -> str:
    """

    :param server:
    :return:
    """
    return settings.DATA_INPUT_DIR_SERVERS[server]


def _get_base_fieldset(data_dir: str) -> FieldSet:
    """

    :param data_dir:
    :return:
    """
    # Defining the folders in which all the data is stored on the different servers
    os.system('echo "Creating the main fieldset"')
    # Loading in the surface currents
    filenames = {'U': data_dir + "HYCOM/HYCOM_Surface*2010*.nc",
                 'V': data_dir + "HYCOM/HYCOM_Surface*2010*.nc",
                 }
    variables = {'U': 'water_u',
                 'V': 'water_v',
                 }
    dimensions = {'U': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'},
                  'V': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'},
                  }

    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=True)
    return fieldset


def _add_stokes_drift(fieldset: FieldSet, input_dir:str):
    """
    Do we include stokes drift yes or no.
    :param fieldset:
    :param input_dir:
    """
    os.system('echo "Adding Stokes drift"')
    filenames = {
                 'Ust': data_dir + "WaveWatchIIIstokes/ww3.2010*_uss.nc",
                 'Vst': data_dir + "WaveWatchIIIstokes/ww3.2010*_uss.nc",
                 }
    variables = {
                 'Ust': 'uuss',
                 'Vst': 'vuss',
                 }
    dimensions = {
                  'Ust': {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'},
                  'Vst': {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'},
                  }
    fieldset_stoke = FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=True)
    fieldset_stoke.Ust.units = GeographicPolar()
    fieldset_stoke.Vst.units = Geographic()
    # Adding the Stokes drift fields to the general fieldset
    fieldset.add_field(fieldset_stoke.Ust)
    fieldset.add_field(fieldset_stoke.Vst)


def _add_border_current(fieldset: FieldSet, input_dir:str):
    """

    :param fieldset:
    :param input_dir:
    """
    os.system('echo "Adding the border current"')
    datasetBor = Dataset(input_dir + 'boundary_velocities_HYCOM.nc')
    borU = datasetBor.variables['MaskUvel'][0, :, :]
    borV = datasetBor.variables['MaskVvel'][0, :, :]
    # Normalizing the border current so that the total current is always 1m/s
    borMag = np.sqrt(np.square(borU) + np.square(borV))
    borMag[borMag == 0] = 1
    borU = np.divide(borU, borMag)
    borV = np.divide(borV, borMag)
    # Adding the actual field
    lonBor, latBor = datasetBor.variables['lon'][:], datasetBor.variables['lat'][:]
    fieldset.add_field(Field('borU', borU, lon=lonBor, lat=latBor, mesh='spherical'))
    fieldset.add_field(Field('borV', borV, lon=lonBor, lat=latBor, mesh='spherical'))
    # making sure the units are interpreted as m s^-1
    fieldset.borU.units = GeographicPolar()
    fieldset.borV.units = Geographic()



def _add_diffusion(fieldset: FieldSet, input_dir: str):
    """

    :param fieldset:
    :param input_dir:
    """
    os.system('echo "Adding diffusion"')
    kh = 10  # m^2 s^-1, following Lacerda et al. (2019) and Liubertseva et al. (2018)
    dataset = Dataset(input_dir + 'HYCOM/HYCOM_Surface_3h_2000-01-01.nc')
    uo = dataset.variables['water_u'][0, 0, :, :]
    lat_kh = dataset.variables['lat'][:]
    lon_kh = dataset.variables['lon'][:]
    kh_f = kh * np.ones(uo.shape)
    kh_f[uo.mask == True] = 0
    fieldset.add_field(Field('Kh_zonal', kh_f, lon=lon_kh, lat=lat_kh, mesh='spherical'))
    fieldset.add_field(Field('Kh_meridional', kh_f, lon=lon_kh, lat=lat_kh, mesh='spherical'))


def _add_landID_field(fieldset: FieldSet, input_dir: str):
    """

    :param fieldset:
    :param input_dir:
    """
    os.system('echo "Adding land/water boolean field"')
    landID=np.load(input_dir+'land_cell_identifier.npy')
    datasetBor = Dataset(input_dir + 'boundary_velocities_HYCOM.nc')
    lonBor, latBor = datasetBor.variables['lon'][:], datasetBor.variables['lat'][:]
    fieldset.add_field(Field('landID', landID,lon=lonBor,lat=latBor,mesh='spherical'))

def _add_distance2shore_field(fieldset: FieldSet, input_dir: str):
    """

    :param fieldset:
    :param input_dir:
    :return:
    """
    os.system('echo "Adding distance to shore"')
    datasetCoast=Dataset(input_dir+'distance2coast.nc')
    distance=datasetCoast.variables['distance'][0,:,:]
    lonD,latD=datasetCoast.variables['lon'][:],datasetCoast.variables['lat'][:]
    fieldset.add_field(Field('distance2shore', distance,lon=lonD,lat=latD,mesh='spherical'))



def _add_wind_field(fieldset: FieldSet, input_dir: str):
    """

    :param fieldset:
    :param input_dir:
    :return:
    """
    os.system('echo "Adding 10m winds"')
    windfiles = glob.glob(input_dir + "Wind/ERA5-wind10m*.nc")
    windfiles.sort()
    filenames = {'u10': windfiles,
                 'v10': windfiles}
    variables = {'u10': 'u10', 'v10': 'v10'}
    dimensions = {'time': 'time', 'lat': 'latitude', 'lon': 'longitude'}
    # Creating a fieldset for the wind data
    fieldset_wind = FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=True)
    fieldset_wind.u10.units = GeographicPolar()
    fieldset_wind.v10.units = Geographic()
    # Adding the wind fields to the general fieldset
    fieldset.add_field(fieldset_wind.u10)
    fieldset.add_field(fieldset_wind.v10)

def _add_seaElevation_field(fieldset: FieldSet, input_dir: str):
    os.system('echo "Adding sea surface elevation"')
    elevfiles = glob.glob(input_dir + "HYCOM/HYCOM_SeaEleve_3h_20*.nc")
    elevfiles.sort()
    filenames = {'eta': elevfiles}
    variables = {'eta': 'surf_el'}
    dimensions = {'time': 'time', 'lat': 'lat', 'lon': 'lon'}
    # Creating a fieldset for the wind data
    fieldset_sea = FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=True)
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

def _compute_ShoreResus_Field(input_dir: str):
    """

    :param SHORE_DEP:
    :param SCENARIO:
    :return:
    """
    if settings.SCENARIO == 'ShoreDependentResuspension':
        s = np.load(input_dir + 'coastline_sand_vs_not_sand.npy')
        if settings.SHORE_DEP == 0:
            resusCor = settings.RESUS_TIME*(0.75+0.25*s)
        if settings.SHORE_DEP == 1:
            resusCor = settings.RESUS_TIME*(0.25+0.75*s)
        return resusCor

def _add_resusTimescale_field(fieldset: FieldSet, input_dir: str):
    """

    :param fieldset:
    :param input_dir:
    :return:
    """
    if settings.SCENARIO_NAME == 'Stochastic':
        p_r = math.exp(-settings.TIME_STEP.total_seconds() / (settings.RESUS_TIME * 86400.))
        fieldset.add_constant('p_resus', p_r)

    elif settings.SCENARIO_NAME == 'ShoreDependentResuspension':
        p_r = np.exp(-settings.TIME_STEP.total_seconds() / (_compute_ShoreResus_Field(input_dir) * 86400.))
        dataset = Dataset(input_dir + 'HYCOM/HYCOM_Surface_3h_2000-01-01.nc')
        lat_s = dataset.variables['lat'][:]
        lon_s = dataset.variables['lon'][:]
        fieldset.add_field(Field('p_resus', p_r, lon=lon_s, lat=lat_s, mesh='spherical'))

def _add_minResuspensionWind_constant(fieldset: FieldSet):
    # The minimum offshore wind speed for resuspension to be possible
    fieldset.add_constant('Wmin', settings.WMIN / 10)

def _add_halo(fieldset: FieldSet):
    os.system('echo "Finally, the periodic halo"')
    fieldset.add_periodic_halo(zonal=True)













