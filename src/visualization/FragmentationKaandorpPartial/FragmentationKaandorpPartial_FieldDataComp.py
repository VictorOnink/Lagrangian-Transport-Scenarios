import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import string
import numpy as np


def FragmentationKaandorpPartial_FieldDataComp(figure_direc, scenario, shore_time, lambda_frag, density,
                                               fig_size=(10, 14), x_label='Size (m)', y_label='Number of Particles',
                                               ax_ticklabel_size=12, ax_label_size=14, legend_size=14):
    # Setting the folder within which we have the output, and where we have the saved data
    output_direc = figure_direc + 'size_distribution/'
    utils.check_direc_exist(output_direc)
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/FragmentationKaandorpPartial/'

    # Getting the sizes of the size classes
    sizes_class = utils.size_range(size_class_number=settings.SIZE_CLASS_NUMBER)

    # Loading the data
    prefix = 'size_distribution'
    data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=scenario, prefix=prefix,
                                                              data_direc=data_direc, shore_time=shore_time,
                                                              lambda_frag=lambda_frag, rho=density, )
    time_index = data_dict['final_index']
    beach_state_list = ['total', 'beach', 'seabed', 'adrift', 'adrift_2m']