from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
from progressbar import ProgressBar
import settings
import utils

if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
    def parcels_to_timeslicing(file_dict: dict):
        """
        For easy animation making, we make slices of the particle lon, lat, depth and beach status so that we can make
        animations even when loading all the files in their entirety would overwhelm my computer RAM

        For the FragmentationKaandorpPartial scenario, we need to consider that we need to load files for different run,
        restart, month and starting year values
        :param file_dict: a dictionary containing all file names for a specific simulation
        :return:
        """
        # Getting the directory saving the output files
        output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeslices/{}/'.format(
            settings.SCENARIO_NAME)
        utils.check_direc_exist(output_direc)
        # Setting the datatime object to which all times are set
        reference_time = datetime(settings.STARTYEAR, 1, 1, 12)
        # loop through the runs
        pbar = ProgressBar()

        for ind_year, year in pbar(enumerate(range(settings.STARTYEAR, settings.STARTYEAR + settings.SIM_LENGTH))):
            for month in range(1, 13):
                for run in range(0, settings.RUN_RANGE):
                    # Loop through the restart files
                    for restart in range(0, settings.SIM_LENGTH - ind_year):
                        # Getting the parcels output file
                        parcels_file = file_dict['parcels'][year][month][run][restart]
                        print(parcels_file)
                        dataset = Dataset(parcels_file)
                        # Loading the particle lon, lat, depth, beach and time
                        full_dict = {'lon': None, 'lat': None, 'z': None, 'beach': None, 'time': None,
                                     'size_class': None}
                        for key in full_dict.keys():
                            full_dict[key] = dataset.variables[key][:, :-1]
                        # Going through all the time steps
                        for timeslice in full_dict['time'][0, :]:
                            slice_dict = {}
                            selection = timeslice == full_dict['time']
                            # Looping for the variables for which I want the time slices
                            for key in ['lon', 'lat', 'z', 'beach', 'size_class']:
                                slice_dict[key] = full_dict[key][selection]
                            # Setting the output name
                            date = (reference_time + timedelta(seconds=timeslice)).strftime("%Y-%m-%d-%H-%M-%S")
                            prefix = 'timeslices_{}'.format(date)
                            output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict['parcels'][settings.STARTYEAR][settings.STARTMONTH][0][0],
                                                                                       prefix=prefix)
                            # If the output file already exists, append the new slices to the previously saved one
                            if utils.check_file_exist(output_name):
                                print('it exists for year {} month {} run {} date {}'.format(year, month, run, date))
                                previous_run_dict = utils.load_obj(filename=output_name)
                                for key in previous_run_dict.keys():
                                    previous_run_dict[key] = np.append(previous_run_dict[key], slice_dict[key])
                                utils.save_obj(filename=output_name, item=previous_run_dict)
                            # Otherwise create a new file to save
                            else:
                                utils.save_obj(filename=output_name, item=slice_dict)

else:
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
        reference_time = datetime(settings.STARTYEAR, 1, 1, 12)
        # loop through the runs
        pbar = ProgressBar()
        for run in pbar(range(settings.RUN_RANGE)):
            # Loop through the restart files
            for restart in range(settings.SIM_LENGTH):
                # Getting the parcels output file
                parcels_file = file_dict[run][restart]
                dataset = Dataset(parcels_file)
                # Loading the particle lon, lat, depth, beach and time
                full_dict = {'lon': None, 'lat': None, 'z': None, 'beach': None, 'time': None}
                for key in full_dict.keys():
                    full_dict[key] = dataset.variables[key][:, :-1]
                if 'size_class' in dataset.variables.keys():
                    full_dict['size_class'] = dataset.variables['size_class'][:, :-1]
                # Going through all the time steps
                for timeslice in full_dict['time'][0, :]:
                    slice_dict = {}
                    selection = timeslice == full_dict['time']
                    # Looping for the variables for which I want the time slices
                    for key in ['lon', 'lat', 'z', 'beach']:
                        slice_dict[key] = full_dict[key][selection]
                    # For the fragmentation runs, getting the particle size
                    if 'size_class' in dataset.variables.keys():
                        slice_dict['size_class'] = full_dict['size_class'][selection]
                    # Setting the output name
                    date = (reference_time + timedelta(seconds=timeslice)).strftime("%Y-%m-%d-%H-%M-%S")
                    prefix = 'timeslices_{}'.format(date)
                    output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0],
                                                                               prefix=prefix)
                    # If this is run=0, then save this as a new file
                    if run == 0:
                        utils.save_obj(filename=output_name, item=slice_dict)
                    # Otherwise, load and append the new slices to the new
                    else:
                        previous_run_dict = utils.load_obj(filename=output_name)
                        for key in previous_run_dict.keys():
                            previous_run_dict[key] = np.append(previous_run_dict[key], slice_dict[key])
                        utils.save_obj(filename=output_name, item=previous_run_dict)

