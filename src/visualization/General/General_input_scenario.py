import settings
import utils
from netCDF4 import Dataset
import numpy as np
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean
import pandas as pd

def General_input_scenario(scenario, figure_direc, figsize=(10, 8), fontsize=14):
    # Setting the folder within which we have the output
    output_direc = figure_direc + 'General/'
    utils.check_direc_exist(output_direc)

    # Getting the bathymetry data
    file_dict = scenario.file_dict
    input_dict = file_dict['STARTFILES_filename']
    df = pd.DataFrame({'lat': np.load(input_dict['lat']), 'lon': np.load(input_dict['lon'])})
    df = df.groupby(['lat', 'lon']).size()
    print(df)
    for index in df.index:
        print('{}, {}, {}'.format(index, df['lat'][index], df['lon'][index]))
