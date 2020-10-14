import numpy as np
import os


def AreaCalc(size_Lat, size_Lon):  # Calculate surface area of grid cells
    deg2rd = np.pi / 180.
    r = 6378.1
    lon_bins = np.linspace(-180, 180, size_Lon + 1)
    lat_bins = np.linspace(-80, 80, size_Lat + 1)
    Area = np.array([[deg2rd * (lon_bins[i + 1] - lon_bins[i]) * (np.sin(deg2rd * lat_bins[j + 1])
                                                                  - np.sin(deg2rd * lat_bins[j])) for i in
                      range(len(lon_bins) - 1)]
                     for j in range(len(lat_bins) - 1)])
    Area = r * r * Area
    return Area  # km^2


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
    if operation != 'count':
        for i in range(np.array(lon_data).shape[0]):
            if weight_data[i] > 0:
                lat_selec, lon_selec = np.argmin(np.abs(lat_data[i] - bins_Lat)), np.argmin(
                    np.abs(lon_data[i] - bins_Lon))
                masses[lat_selec, lon_selec] += weight_data[i]
                counts[lat_selec, lon_selec] += 1
        if operation == 'mean':
            masses[counts > 0] = np.divide(masses[counts > 0], counts[counts > 0])
        if area_correc == True:
            masses = np.divide(masses, AreaCalc(size_Lat=len(bins_Lat), size_Lon=len(bins_Lon)))
        return masses  # weight / km^2
    elif operation == 'count':
        for i in range(np.array(lon_data).shape[0]):
            if weight_data[i] > 0:
                lat_selec, lon_selec = np.argmin(np.abs(lat_data[i] - bins_Lat)), np.argmin(
                    np.abs(lon_data[i] - bins_Lon))
                counts[lat_selec, lon_selec] += 1
        if area_correc == True:
            counts = np.divide(counts, AreaCalc(size_Lat=len(bins_Lat), size_Lon=len(bins_Lon)))
        return counts  # counts / km^2


def _analysis_save_file_name(input_file: str, prefix: str, out_type: str = '.mat'):
    _, file_name = os.path.split(input_file)
    return prefix + '_' + file_name.split('_r=')[0] + out_type


def _particles_in_domain(domain, lon, lat, weight=0, beach=0, time=0,
                         return_sub_array=False):
    lon_min, lon_max, lat_min, lat_max = domain
    # Select only particles that are within the domain that we are interested in
    select = (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)
    # And now we return the time, weight and beach arrays of the particles that are within
    # the domain of interest
    if return_sub_array==True:
        return weight[select], beach[select], time[select]
    else:
        return select
