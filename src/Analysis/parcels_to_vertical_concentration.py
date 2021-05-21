import settings as settings
import utils
from advection_scenarios import advection_files
from netCDF4 import Dataset
import numpy as np
import progressbar
import os


def parcels_to_vertical_concentration(file_dict: dict):
    # Get the depth range for the region in question, and using that to set the bins
    advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                        advection_scenario=settings.ADVECTION_DATA,
                                                        repeat_dt=None)
    adv_file_dict = advection_scenario.file_names
    depth_min, depth_max = np.nanmin(adv_file_dict['DEPTH']), np.nanmax(adv_file_dict['DEPTH'])
    step = 0.5
    depth_bins = np.arange(depth_min, depth_max + step, step)

    # Getting the directory saving the output files
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'concentrations/{}/'.format(settings.SCENARIO_NAME)
    utils.check_direc_exist(output_direc)

    # Create the output dictionary, and a dictionary to keep track of particle counts for the normalization
    output_dict = {'depth': depth_bins}
    counts_dict = {}
    for simulation_year in range(settings.SIM_LENGTH):
        key_year = utils.analysis_simulation_year_key(simulation_year)
        output_dict[key_year] = {}
        counts_dict[key_year] = {}
        for month in range(12):
            output_dict[key_year][month] = np.zeros(len(depth_bins) - 1, dtype=np.float32)
            counts_dict[key_year][month] = 0.0

    # Looping through all the simulation years and runs
    for run in progressbar.progressbar(range(settings.RUN_RANGE)):
        # Loop through the restart files
        for restart in range(settings.SIM_LENGTH):
            key_year = utils.analysis_simulation_year_key(restart)
            # Load the depth data, with monthly intervals
            parcels_file = file_dict[run][restart]
            dataset = Dataset(parcels_file)
            depth = dataset.variables['z'][:]
            beach = dataset.variables['beach'][:]
            time_interval = depth.shape[1] // 365 * 31
            depth = depth[:, ::time_interval]
            beach = beach[:, ::time_interval]
            # Start looping through the different time intervals
            for month in range(depth.shape[1]):
                # Checking to make sure we don't have any nan values, or particles with beach==2 (indicating a particle
                # that was previously removed)
                non_nan_values = (~np.isnan(depth[:, month])) & (beach[:, month] != 2)
                depth_month = depth[:, month][non_nan_values]
                # Now, calculating the vertical histogram
                histogram_counts, _ = np.histogram(depth_month, bins=depth_bins)
                output_dict[key_year][month] += histogram_counts
                # Updating the number of particles used within a certain profile so that we can normalize later
                counts_dict[key_year][month] += depth_month.size

    # Now, we loop through all the profiles and normalize them by the number of particles within the profile
    for key_year in counts_dict.keys():
        for key_month in counts_dict[key_year].keys():
            output_dict[key_year][key_month] /= counts_dict[key_year][key_month]

    # Saving the output
    # Saving the computed concentration
    prefix = 'vertical_concentration'
    output_name = output_direc + utils._analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
    utils.save_obj(output_name, output_dict)
    os.system('echo "The vertical concentration has been saved"')




