import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string


def FragmentationKaandorpPartial_Concentration(scenario, figure_direc, rho, shore_time, beach_state, simulation_year,
                                               lambda_frag,
                                               figsize=(20, 10), ax_ticklabel_size=12, ax_label_size=14):
    # Setting the folder within which we have the output, and where we have the saved timeslices
    output_direc = figure_direc + 'concentrations/'
    data_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format('SizeTransport')

    # Loading in the data
    prefix = 'horizontal_concentration'
    concentration_dict = {}
    if simulation_year == 'average':
        key_concentration = "overall_concentration"
    else:
        key_concentration = utils.analysis_simulation_year_key(simulation_year)
    data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=scenario, prefix=prefix,
                                                              data_direc=data_direc, shore_time=shore_time,
                                                              lambda_frag=lambda_frag, rho=rho, postprocess=True)
    concentration_dict = data_dict[key_concentration][beach_state]
    Lon, Lat = np.meshgrid(data_dict['lon'], data_dict['lat'])

    # Normalizing the concentration by the lowest non-zero concentration over all the sizes
    normalization_factor = 1e10
    for size in concentration_dict.keys():
        concentration = concentration_dict[size]
        min_non_zero = np.nanmin(concentration[concentration > 0])
        if min_non_zero < normalization_factor:
            normalization_factor = min_non_zero
    for size in concentration_dict.keys():
        concentration_dict[size] /= normalization_factor
