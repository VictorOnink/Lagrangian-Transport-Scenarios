from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
from progressbar import ProgressBar
import settings
import utils


def parcels_to_timeslicing(file_dict: dict):
    """
    For easy animation making, we make slices of the particle lon, lat, depth and beach status so that we can make
    animations even when loading all the files in their entirety would overwhelm my computer RAM
    :param file_dict: a dictionary containing all file names for a specific simulation
    :return:
    """
    # Getting the directory saving the output files
    output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeslices/{}/'.format(
        settings.SCENARIO_NAME)
    utils.check_direc_exist(output_direc)
    # Setting the datatime object to which all times are set
    reference_time = datetime(settings.START_YEAR, 1, 1, 12)
    # loop through the runs
    for run in range(settings.RUN_RANGE):
        # Loop through the restart files
        for restart in range(settings.SIM_LENGTH):
            # Getting the parcels output file
            parcels_file = file_dict[run][restart]
            dataset = Dataset(parcels_file)
            # Loading the particle lon, lat, depth, beach and time
            full_dict = {'lon': None, 'lat': None, 'z': None, 'beach': None, 'time': None}
            for key in full_dict.keys():
                full_dict[key] = dataset.variables[key][:, :-1]
            # Going through all the time steps
            for timeslice in full_dict['time'][0, :]:
                slice_dict = {}
                selection = timeslice == full_dict['time']
                # Looping for the variables for which I want the time slices
                for key in ['lon', 'lat', 'z', 'beach']:
                    slice_dict[key] = full_dict[key][selection]
                # Setting the output name
                date = (reference_time + timedelta(seconds=timeslice)).strftime("%Y-%m-%d-%H-%M-%S")
                prefix = 'timeslices_{}'.format(date)
                output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
                # If this is a new file, then save this
                if not utils.check_file_exist(output_name + '.pkl'):
                    utils.save_obj(filename=output_name, item=slice_dict)
                # If the file already exists, then this is a previous run that we want to then append to
                else:
                    previous_run_dict = utils.load_obj(filename=output_name)
                    for key in previous_run_dict.keys():
                        previous_run_dict[key] = np.append(previous_run_dict[key], slice_dict[key])
                    utils.save_obj(filename=output_name, item=previous_run_dict)

