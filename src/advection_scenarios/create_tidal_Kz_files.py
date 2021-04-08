import numpy as np
from numpy import array
import scipy.interpolate as interpolate
import xarray
from netCDF4 import Dataset
import settings
import utils


def create_tidal_Kz_files(file_name: str, LON: array, LAT: array, DEPTH: array, BATH_filenames: str,
                          BATH_variables: dict):
    """
    Determining Kz from the estimates of production of turbulent kinetic energy, as outlined in
    de Lavergne et al. (2020) and based on code shared by Clement Vic. The input file for TIDAL_filename can be
    downloaded at https://www.seanoe.org/data/00619/73082/
    """
    # Loading the global data from de Lavergne et al. (2020)
    TIDAL_filename = utils.get_input_directory(server=settings.SERVER) + 'global_tidal_energy_dissipation.nc'
    TIDAL_data = {}
    for key in Dataset(TIDAL_filename).variables.keys():
        TIDAL_data[key] = Dataset(TIDAL_filename).variables[key][:]

    # Load the bathymetry data
    bathymetry = Dataset(BATH_filenames).variables[BATH_variables['DEPTH']][:]
    bathymetry[bathymetry.mask] = 0

    # Computing Kz on the TIDAL_data grid according to Kv = gamma * epsilon / N^2
    gamma = 0.2  # Mixing efficiency
    TIDAL_Kz = np.divide(gamma * TIDAL_data['epsilon_tid'], TIDAL_data['buoyancy_frequency_squared'])

    # The TIDAL_data gridding isn't regular in the z-direction. We will first interpolate the TIDAL_Kz fields onto the
    # DEPTH levels
    TIDAL_Kz_inter = interpolate_to_DEPTH(TIDAL_Kz, TIDAL_data, DEPTH)

    # # Now we interpolate the Kz values onto the model grid defined by LON, LAT and DEPTH
    # GRID_Kz = interpolate_to_GRID(TIDAL_Kz_inter, DEPTH, LON, LAT, TIDAL_data)

    # Due to very low N^2 values (corresponding to very weak stratification), there are some regions where Kz is
    # unfeasibly high (Kz > 100 m^2/s). Therefore, I cap GRID_Kz at 0.1 m^2/s. This only affects 0.08% of all the cells
    # in the TIDAL_data, and prevents numerical issues later on.
    TIDAL_Kz[TIDAL_Kz > 1e-1] = 1e-1

    # Computing the TIDAL_Kz gradient
    print(TIDAL_Kz.shape)
    print(DEPTH.shape)
    TIDAL_dKz = np.gradient(TIDAL_Kz, DEPTH, axis=0)
    print('max {}, min {}, mean {}'.format(np.nanmax(TIDAL_dKz), np.nanmin(TIDAL_dKz), np.nanmean(TIDAL_dKz)))

    # Saving the field to a .nc file
    coords = [('time', np.array([0])), ('depth', DEPTH), ('lat', LAT), ('lon', LON)]
    dcoo = {'time': np.array([0]), 'depth': DEPTH, 'lat': LAT, 'lon': LON}
    dset = xarray.Dataset({'TIDAL_Kz': xarray.DataArray(TIDAL_Kz[np.newaxis, ...], coords=coords)}, coords=dcoo)
    dset.to_netcdf(file_name)


def interpolate_to_DEPTH(TIDAL_Kz: array, TIDAL_data: dict, DEPTH: array):
    TIDAL_depth = TIDAL_data['depth_midpoint']
    TIDAL_depth[TIDAL_depth.mask] = np.nanmax(TIDAL_depth)
    TIDAL_Kz[TIDAL_Kz.mask] = 0
    TIDAL_Kz_inter = np.zeros((DEPTH.shape[0], TIDAL_Kz.shape[1], TIDAL_Kz.shape[2]))
    for lat in range(TIDAL_Kz.shape[1]):
        for lon in range(TIDAL_Kz.shape[2]):
            inter_1D = interpolate.interp1d(TIDAL_depth[:, lat, lon], TIDAL_Kz[:, lat, lon], bounds_error=False)
            TIDAL_Kz_inter[:, lat, lon] = inter_1D(np.array(DEPTH))
    # TIDAL_data only starts at z=5m, so we will assume homogenous mixing there for now to prevent issues coming up
    # later. This affects the first two rows of the array
    TIDAL_min_depth = np.nanmin(TIDAL_depth, axis=(0, 1, 2))
    for z, z_level in enumerate(DEPTH):
        if z_level < TIDAL_min_depth:
            TIDAL_Kz_inter[z, :, :] = TIDAL_Kz[0, :, :]

    return TIDAL_Kz_inter


def interpolate_to_GRID(TIDAL_Kz_inter: array, DEPTH: array, LON: array, LAT: array, TIDAL_data: dict):
    GRID_Kz = np.zeros((DEPTH.shape[0], LAT.shape[0], LON.shape[0]))
    # Set all the masked TIDAL_Kz_inter values to 0
    TIDAL_Kz_inter[np.isnan(TIDAL_Kz_inter)] = 0
    for z_level in range(DEPTH.shape[0]):
        T_LAT, T_LON = np.meshgrid(TIDAL_data['lat'], TIDAL_data['lon'], sparse=True)
        # For some reason the interpolation function requires the transpose of the 2D array, but I'm not 100% sure why
        inter_2D = interpolate.interp2d(T_LAT, T_LON, TIDAL_Kz_inter[z_level, :, :].T)
        GRID_Kz[z_level, :, :] = inter_2D(LON, LAT)
    return GRID_Kz
