import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import numpy as np
import string
from datetime import datetime, timedelta
import matplotlib.dates as mdates


def FragmentationKaandorpPartial_timeseries(scenario, figure_direc, shore_time, lambda_frag, rho,
                                           figsize=(12, 10), ax_label_size=18, tick_label_size=16):
    # Setting the folder within which we have the output, and where we have the saved timeslices
    output_direc = figure_direc + 'timeseries/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/FragmentationKaandorpPartial/'
    utils.check_direc_exist(output_direc)

    # Loading in the data
    prefix = 'timeseries'
    timeseries_dict = {}
    # beach_state_list = ['beach', 'afloat', 'seabed', 'removed', 'total']
    beach_state_list = ['beach', 'afloat', 'seabed', 'total']
    for size_class in range(settings.SIZE_CLASS_NUMBER):
        timeseries_dict[size_class] = {}
        data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=scenario, prefix=prefix, data_direc=data_direc,
                                                                  shore_time=shore_time, lambda_frag=lambda_frag, rho=rho)
        for beach_state in beach_state_list:
            timeseries_dict[size_class][beach_state] = data_dict[beach_state][size_class]
        timeseries_dict[size_class]['time_raw'] = data_dict['time'][size_class]

    utils.print_statement('This worked', to_print=True)
