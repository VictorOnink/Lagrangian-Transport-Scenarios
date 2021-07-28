import numpy as np
import os


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
    dlon, dlat = lon1_rad - lon2_rad, lat1_rad - lat2_rad
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


def histogram(lon_data, lat_data, bins_Lon, bins_Lat, weight_data=0,
              area_correc=True, operation='sum'):
    """
    operation -> 'sum' = get array with sum of weights per cell, divided by km^2
                 'mean' = get array with mean of weights per cell, divided by km^2
                 'count' = get array with counts of occurences within each cell, divided by km^2
    :param lon_data: Nx1 or (N,) array
    :param lat_data: Nx1 or (N,) array
    :param bins_Lon: LONx1 array
    :param bins_Lat: LATx1 array
    :param weight_data: 0 if we are only interested in counts, else Nx1 or (N,) array
    :return:
    """
    lon_data, lat_data = lon_data.reshape(np.size(lon_data)), lat_data.reshape(np.size(lat_data))
    weight_data = weight_data.reshape(np.size(weight_data))
    masses = np.zeros((len(bins_Lat), len(bins_Lon)))
    counts = np.zeros((len(bins_Lat), len(bins_Lon)))
    if operation == 'count':
        for i in range(np.array(lon_data).shape[0]):
            if weight_data[i] > 0:
                lat_selec, lon_selec = np.argmin(np.abs(lat_data[i] - bins_Lat)), np.argmin(
                    np.abs(lon_data[i] - bins_Lon))
                counts[lat_selec, lon_selec] += 1
        if area_correc:
            counts = np.divide(counts, surface_area_grid(bins_Lat, bins_Lon))
        return counts  # counts / km^2
    else:
        for i in range(np.array(lon_data).shape[0]):
            if weight_data[i] > 0:
                lat_selec, lon_selec = np.argmin(np.abs(lat_data[i] - bins_Lat)), np.argmin(
                    np.abs(lon_data[i] - bins_Lon))
                masses[lat_selec, lon_selec] += weight_data[i]
                counts[lat_selec, lon_selec] += 1
        if operation == 'mean':
            masses[counts > 0] = np.divide(masses[counts > 0], counts[counts > 0])
        if area_correc:
            masses = np.divide(masses, surface_area_grid(bins_Lat, bins_Lon))
        return masses  # weight / km^2


def analysis_save_file_name(input_file: str, prefix: str, suffix=None):
    _, file_name = os.path.split(input_file)
    file_name = prefix + '_' + file_name.split('_r=')[0]
    if suffix is not None:
        file_name += suffix
    return file_name


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