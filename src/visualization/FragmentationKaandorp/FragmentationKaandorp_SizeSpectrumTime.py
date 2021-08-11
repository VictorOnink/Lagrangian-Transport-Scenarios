import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def FragmentationKaandorp_SizeSpectrumTime(figure_direc, scenario, shore_time, lambda_frag_list, figsize=(12, 10),
                                           fontsize=12):
    # Setting the folder within which we have the output, and where we have the saved data
    output_direc = figure_direc + 'size_distribution/'
    utils.check_direc_exist(output_direc)
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/FragmentationKaandorp/'

    # Loading in the data
    prefix = 'size_distribution'
    size_dict = {}
    for lambda_frag in lambda_frag_list:
        data_dict = vUtils.FragmentationKaandorp_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                                           shore_time=shore_time, lambda_frag=lambda_frag, rho=920)
        size_dict[lambda_frag] = {}
        for month in data_dict.keys():
            if month != 'size_bins':
                size_dict[lambda_frag][month] = data_dict[month]
    size_bins = data_dict['size_bins']
    for keys in size_dict.keys():
        utils.print_statement(keys, to_print=True)

    # Creating the figure
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=2, ncols=4, width_ratios=[1, 1, 1, 0.2])
    gs.update(wspace=0.2, hspace=0.2)

    ax_list = []
    for row in range(gs.nrows):
        for column in range(gs.ncols - 1):
            ax = fig.add_subplot(gs[row, 0])
            if column == 0:
                ax.set_ylabel(r'Particles', fontsize=fontsize)
            if row == (gs.nrows - 1):
                ax.set_xlabel(r'Size (m)', fontsize=fontsize)
            ax.set_xlim([size_bins.min(), size_bins.max()])
            ax.set_ylim([0, 1e4])
            ax.tick_params(which='major', length=7)
            ax.tick_params(which='minor', length=3)
            ax.set_yscale('log')
            ax.set_xscale('log')
            if row != (gs.nrows - 1):
                ax.set_xticklabels([])
            ax_list.append(ax)

    file_name = output_direc + 'SizeSpectrumTime-ST={}.png'.format(shore_time)
    plt.savefig(file_name, bbox_inches='tight')
