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
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeslices/{}/'.format('SizeTransport')
    utils.check_direc_exist(output_direc)

    # Loading in the data
    date = datetime(2012, 12, 31, 12, 0, 0).strftime("%Y-%m-%d-%H-%M-%S")
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
        depth_histogram, _ = np.histogram(seabed_depths, bins=depth_bins)
        # Normalizing the histogram by the total number of seabed particles
        depth_histogram = depth_histogram.astype(float)
        depth_histogram /= np.nansum(depth_histogram)
        depth_histogram *= 100.
        # Saving the histogram into the dictionary
        depth_histogram_dict[size] = depth_histogram

    # Creating the figure structure
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=1, ncols=1)
    ax = fig.add_subplot(gs[0, 0])

    ax.set_ylabel(r'Fraction of Total (%)', fontsize=fontsize)
    ax.set_ylim([0, 100])
    ax.set_xlabel('Depth (m)', fontsize=fontsize)
    ax.set_xlim([0, 200])

    # Plotting the data
    for index_size, size in enumerate(size_list):
        ax.plot(depth_bins[:-1], depth_histogram_dict[size], linestyle='-',
                color=vUtils.discrete_color_from_cmap(index_size, subdivisions=len(size_list)),
                label=size_label(size))
    # And adding in a legend
    ax.legend(fontsize=fontsize, loc='upper right')

    file_name = output_direc + 'Seabed_depth_histogram.png'
    plt.savefig(file_name, bbox_inches='tight')


def size_label(size):
    return r'r = {:.1E} mm'.format(size * 1e3)