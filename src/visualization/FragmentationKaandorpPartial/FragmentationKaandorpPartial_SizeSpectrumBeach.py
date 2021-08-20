import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import string
import numpy as np


def FragmentationKaandorpPartial_SizeSpectrumBeach(figure_direc, scenario, shore_time, lambda_frag_list, density,
                                                  fig_size=(18, 12), x_label='Size (m)', y_label='Number of Particles',
                                                  ax_ticklabel_size=12, ax_label_size=14, legend_size=14):
    # Setting the folder within which we have the output, and where we have the saved data
    output_direc = figure_direc + 'size_distribution/'
    utils.check_direc_exist(output_direc)
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/FragmentationKaandorpPartial/'

    # Loading in the data
    prefix = 'size_distribution'
    size_dict = {}
    for lambda_frag in lambda_frag_list:
        data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=scenario, prefix=prefix,
                                                                  data_direc=data_direc, shore_time=shore_time,
                                                                  lambda_frag=lambda_frag, rho=density)
        size_dict[lambda_frag] = {}
        for month in data_dict.keys():
            if month != 'size_bins':
                size_dict[lambda_frag][month] = data_dict[month]
    size_bins = data_dict['size_bins'][:-1]
