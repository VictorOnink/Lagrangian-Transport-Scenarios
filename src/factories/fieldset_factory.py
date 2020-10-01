import os
import src.settings as settings
from netCDF4 import Dataset
import numpy as np
from parcels import GeographicPolar, Geographic,FieldSet,Field


class FieldSetFactory():
    """"""
    @classmethod
    def create_fieldset(cls, server: int, stokes: int,
                        stokes_drift: bool = False,
                        border_current: bool = False,
                        diffusion: bool = False,
                        bool_field: bool = False):
        """

        :param server:
        :param stokes:
        :param stokes_drift:
        :param border_current:
        :param diffusion:
        :param bool_field:
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
        if bool_field:
            _add_bool_field(fieldset=fieldset, input_dir=input_dir)
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
    # Loading in the surface currents and Stokes drift
    filenames = {'U': data_dir + "HYCOM/HYCOM_Surface*2010*.nc",
                 'V': data_dir + "HYCOM/HYCOM_Surface*2010*.nc",
                 'Ust': data_dir + "WaveWatchIIIstokes/ww3.2010*_uss.nc",
                 'Vst': data_dir + "WaveWatchIIIstokes/ww3.2010*_uss.nc",
                 }
    variables = {'U': 'water_u',
                 'V': 'water_v',
                 'Ust': 'uuss',
                 'Vst': 'vuss',
                 }
    dimensions = {'U': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'},
                  'V': {'time': 'time', 'depth': 'depth', 'lat': 'lat', 'lon': 'lon'},
                  'Ust': {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'},
                  'Vst': {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'},
                  }

    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=True)
    return fieldset


def _add_stokes_drift(fieldset: FieldSet, input_dir:str):
    """

    :param fieldset:
    :param input_dir:
    """
    os.system('echo "Adding Stokes drift"')
    fieldset.Ust.units = GeographicPolar()
    fieldset.Vst.units = Geographic()


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
    # making sure the units are interpreted as 1 m/s
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


def _add_bool_field(fieldset: FieldSet, input_dir: str):
    """

    :param fieldset:
    :param input_dir:
    """
    os.system('echo "Adding land/water boolean field"')
    landID=np.load(input_dir+'land_cell_identifier.npy')
    datasetBor = Dataset(input_dir + 'boundary_velocities_HYCOM.nc')
    lonBor, latBor = datasetBor.variables['lon'][:], datasetBor.variables['lat'][:]
    fieldset.add_field(Field('landID', landID,lon=lonBor,lat=latBor,mesh='spherical'))













