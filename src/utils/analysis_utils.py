import numpy as np
import os

import settings


def surface_area_grid(lat, lon):
    """
    Calculate the area of each grid cell
    Area is in square kilometers
    Input
    -----------
    lat: vector of latitude in degrees
    lon: vector of longitude in degrees
    Output
    -----------
    area: grid-cell area in square-kilometers with dimensions, [lat,lon]
    Notes
    -----------
    Based on the function in
    https://github.com/chadagreene/CDT/blob/master/cdt/cdtarea.m
    This particular function was copied from
    https://towardsdatascience.com/the-correct-way-to-average-the-globe-92ceecd172b7
    """
    xlon, ylat = np.meshgrid(lon, lat)
    R = earth_radius(ylat)

    dlat = np.deg2rad(np.gradient(ylat, axis=0))
    dlon = np.deg2rad(np.gradient(xlon, axis=1))

    dy = dlat * R
    dx = dlon * R * np.cos(np.deg2rad(ylat))

    area = dy * dx

    return area


def earth_radius(lat):
    '''
    calculate radius of Earth assuming oblate spheroid
    defined by WGS84
    Input
    ---------
    lat: vector or latitudes in degrees
    Output
    ----------
    r: vector of radius in kilometers
    Notes
    -----------
    WGS84: https://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350.2-a/Chapter%203.pdf
    This particular function was copied from
    https://towardsdatascience.com/the-correct-way-to-average-the-globe-92ceecd172b7
    '''
    # define oblate spheroid from WGS84
    a = 6378137
    b = 6356752.3142
    e2 = 1 - (b ** 2 / a ** 2)
    # convert from geodecic to geocentric
    # see equation 3-110 in WGS84
    lat = np.deg2rad(lat)
    lat_gc = np.arctan((1 - e2) * np.tan(lat))
    # radius equation
    # see equation 3-107 in WGS84
    r = (
            (a * (1 - e2) ** 0.5)
            / (1 - (e2 * np.cos(lat_gc) ** 2)) ** 0.5
    )
    return r / 1000.


def distance_between_points(lon1, lat1, lon2, lat2, units='km'):
    """
    This function computes the distance between two sets of points using the Haversine approximation. This does assume
    the Earth is a perfect sphere, but this is an approximation we'll accept for a relatively cheap computation

    This function is designed to work both with float lon/lat values and numpy arrays of lon/lat values, but then we do
    need to have that the sizes of lon1, lat1, lon2 and lat2 are the same

    Based on: https://stackoverflow.com/questions/19412462/getting-distance-between-two-points-based-on-latitude-longitude/43211266#43211266
    :param lon1: array or float value containing first set of lon values
    :param lat1: array of float value containing first set of lat values
    :param lon2: array or float value containing second set of lon values
    :param lat2: array or float value containing second set of lat values
    :param units: distance in km or m
    :return:
    """
    # Check variable types
    variable_assertion = 'The variable types do not match'
    assert type(lon1) == type(lat1) and type(lon1) == type(lon2) and type(lon1) == type(lat2), variable_assertion

    # Check if the numpy arrays are the correct size if the coordinates are indeed numpy arrays
    if type(lon1) == np.ndarray:
        size_assertion = 'The coordinate arrays do not have the same sizes'
        assert lon1.size == lon2.size and lon1.size == lat1.size and lon1.size == lat2.size, size_assertion

    # Check units are either km or m
    assert units in ['km', 'm'], "Please choose either 'km' or 'm' for your units"

    # Computing the difference in coordinates in radians
    lon1_rad, lat1_rad = np.deg2rad(lon1), np.deg2rad(lat1)
    lon2_rad, lat2_rad = np.deg2rad(lon2), np.deg2rad(lat2)
    dlon, dlat = np.abs(lon1_rad - lon2_rad) % (2 * np.pi), np.abs(lat1_rad - lat2_rad) % np.pi

    # Converting that to a physical distance
    radius = earth_radius(lat=0)
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2) ** 2
    angle = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distance = radius * angle

    # Returning the distance
    if units == 'km':
        return distance
    else:
        return distance * 1000.0


def histogram(lon_data, lat_data, bins_Lon, bins_Lat, weight_data=None, area_correc=True, operation='sum'):
    """
    :parammoperation: 'sum' = get array with sum of weights per cell, divided by km^2
                      'mean' = get array with mean of weights per cell, divided by km^2
                      'count' = get array with counts of occurrences within each cell, divided by km^2
    :param area_correc: if we want to correct the concentration for the surface area of the cell
    :param lon_data: Nx1 or (N,) array
    :param lat_data: Nx1 or (N,) array
    :param bins_Lon: LONx1 array
    :param bins_Lat: LATx1 array
    :param weight_data: Nx1 or (N,) array
    :return:
    """
    # Calculating the bin midpoints
    bin_mid_lon = 0.5 * bins_Lon[1:] + 0.5 * bins_Lon[:-1]
    bin_mid_lat = 0.5 * bins_Lat[1:] + 0.5 * bins_Lat[:-1]

    # If we only want the counts, then we set the weights to one
    if operation == 'count':
        weight_data = np.ones(lon_data.shape, dtype=float)

    # Calculating the bin concentration
    bin_concentrations, _, _ = np.histogram2d(lat_data, lon_data, bins=[bins_Lat, bins_Lon], weights=weight_data)

    # If we want the mean, then we divide the bin_concentrations by the number of instances per bin
    if operation == 'mean':
        bin_counts, _, _ = np.histogram2d(lon_data, lat_data, bins=[bins_Lon, bins_Lat],
                                          weights=np.ones(lon_data.shape, dtype=float))
        bin_concentrations = np.divide(bin_concentrations, bin_counts)

    if area_correc:
        bin_concentrations = np.divide(bin_concentrations, surface_area_grid(bin_mid_lat, bin_mid_lon))

    return bin_concentrations, bin_mid_lat, bin_mid_lon


def analysis_save_file_name(input_file: str, prefix: str, suffix=None, split=None):
    _, file_name = os.path.split(input_file)
    if split is None:
        split = file_name_string_split()
    file_name = prefix + '_' + file_name.split(split)[0]
    if suffix is not None:
        file_name += suffix
    return file_name


def file_name_string_split():
    if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
        return '_y='
    else:
        return '_r='


def particles_in_domain(domain, lon, lat):
    lon_min, lon_max, lat_min, lat_max = domain
    # Select only particles that are within the domain that we are interested in
    select = (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)
    # Now return the boolean array indicating which particles are within the domain of interest
    return select


def dict_key_vertical_concentration(restart, month):
    return 'year_{}_month_{}'.format(restart, month)


def analysis_simulation_year_key(simulation_years):
    return 'year_{}'.format(simulation_years)


def init_size_key(size):
    return 'size_{:.1E}'.format(size)


def size_range(size_class_number=None, init_size=settings.INIT_SIZE, single_size_class=None, units='m'):
    if single_size_class is not None:
        assert type(single_size_class) == int, 'The size class number must be an integer'
        output = init_size * 2 ** -single_size_class
    else:
        output = np.ones(size_class_number, dtype=float) * init_size
        for size_class in range(size_class_number):
            output[size_class] *= 2 ** -size_class
    if units == 'mm':
        output *= 1e3
    return output
