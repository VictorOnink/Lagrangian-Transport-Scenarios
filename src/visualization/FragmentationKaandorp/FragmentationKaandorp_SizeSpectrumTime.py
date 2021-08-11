import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def FragmentationKaandorp_SizeSpectrumTime(figure_direc, scenario, shore_time, lambda_frag_list):
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