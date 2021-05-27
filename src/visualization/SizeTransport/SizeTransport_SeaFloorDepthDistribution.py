import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import numpy as np
import string
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from advection_scenarios import advection_files


def SizeTransport_SeaFloorDepthDistribution(scenario, figure_direc, size_list, rho_list, figsize=(10, 10), fontsize=16,
                                            histogram=False, cumulative=False):
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
    if histogram:
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
    if cumulative:
        depth_cumulative_dict = {}
        for size in size_list:
            # Loading in the particle depths, and the beach states
            beach_state, all_depths = timeseries_dict[size]['beach'], timeseries_dict[size]['z']
            # Selecting just the particle depths for those at the seabed
            seabed_depths = all_depths[beach_state == 3]
            # Initializing the cumulative counts array
            depth_cumulative_dict[size] = np.zeros(depth_bins.shape, dtype=float)
            # Now looping through the depths, and seeing how many particles are stuck at each depth
            for index, depth in enumerate(depth_bins):
                depth_cumulative_dict[size][index] += np.nansum(seabed_depths < depth)
            # Normalizing by the total number of particles on the sea bed
            depth_cumulative_dict[size] /= len(seabed_depths)
            depth_cumulative_dict[size] *= 100.


    # Creating the figure structure
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=1, ncols=1)
    ax = fig.add_subplot(gs[0, 0])

    ax.set_ylabel(r'Fraction of Total (%)', fontsize=fontsize)
    ax.labelsize(fontsize)
    ax.set_ylim([0, 105])
    ax.set_xlabel('Depth (m)', fontsize=fontsize)
    if histogram:
        ax.set_xlim([0, 50])
    elif cumulative:
        ax.set_xscale('log')
        ax.set_xlim([1e0, 1e3])


    # Plotting the data
    if histogram:
        for index_size, size in enumerate(size_list):
            ax.plot(depth_bins[:-1], depth_histogram_dict[size], linestyle='-',
                    color=vUtils.discrete_color_from_cmap(index_size, subdivisions=len(size_list)),
                    label=size_label(size))
    elif cumulative:
        for index_size, size in enumerate(size_list):
            ax.plot(depth_bins, depth_cumulative_dict[size], linestyle='-',
                    color=vUtils.discrete_color_from_cmap(index_size, subdivisions=len(size_list)),
                    label=size_label(size))
    # And adding in a legend
    if histogram:
        ax.legend(fontsize=fontsize, loc='upper right')
    elif cumulative:
        ax.legend(fontsize=fontsize, loc='lower right')

    if histogram:
        file_name = output_direc + 'Seabed_depth_histogram.png'
    elif cumulative:
        file_name = output_direc + 'Seabed_depth_cumulative.png'
    plt.savefig(file_name, bbox_inches='tight')


def size_label(size):
    return r'r = {:.2f} mm'.format(size * 1e3)