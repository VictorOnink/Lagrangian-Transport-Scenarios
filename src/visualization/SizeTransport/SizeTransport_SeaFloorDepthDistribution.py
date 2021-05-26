import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import numpy as np
import string
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from advection_scenarios import advection_files


def SizeTransport_SeaFloorDepthDistribution(scenario, figure_direc, size_list, rho_list, figsize=(10, 10), fontsize=12):
    # Getting the size of the domain that we want to plot for
    advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                        advection_scenario='CMEMS_MEDITERRANEAN',
                                                        repeat_dt=None)
    adv_file_dict = advection_scenario.file_names

    # Setting the folder within which we have the output, and where we have the saved timeslices
    output_direc = figure_direc + 'mix/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/{}/'.format('SizeTransport')
    utils.check_direc_exist(output_direc)

    # Loading in the data
    date = datetime(2013, 12, 31, 12, 0, 0).strftime("%Y-%m-%d-%H-%M-%S")
    prefix = 'timeslices_{}'.format(date)

    timeseries_dict = {}
    for index, size in enumerate(size_list):
        data_dict = vUtils.SizeTransport_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                                   size=size, rho=rho_list[index])
        timeseries_dict[size] = {}
        for variable in ['lon', 'lat', 'beach', 'z']:
            timeseries_dict[size][variable] = data_dict[variable]

    # Beach histogram:
    depth_bins = np.arange(0, np.nanmax(adv_file_dict['DEPTH']), 1.0)
    depth_histogram_dict = {}
    for size in size_list:
        # Loading in the particle depths, and the beach states
        beach_state, all_depths = timeseries_dict[size]['beach'], timeseries_dict[size]['z']
        # Selecting just the particle depths for those at the seabed
        seabed_depths = all_depths[beach_state == 3]
        # Getting the histogram for the number of particles in which depth bins
        seabed_depths, _ = np.histogram(a=seabed_depths, bins=depth_bins)
        # Normalizing the histogram by the total number of seabed particles
        seabed_depths /= np.nansum(seabed_depths)
        # Saving the histogram into the dictionary
        depth_histogram_dict[size] = seabed_depths