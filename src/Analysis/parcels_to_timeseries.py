import settings as settings
import utils
from netCDF4 import Dataset
import numpy as np
from copy import deepcopy

if settings.SCENARIO_NAME in ['FragmentationKaandorpPartial']:
    def parcels_to_timeseries(file_dict: dict, lon_min: float = -180, lon_max: float = 180, lat_min: float = -90,
                              lat_max: float = 90):
        """
        This has the same basic function as the parcels_to_timeseries, but it looks at the particlenumber variable and the
        various size bins as implemented in the FragmentationKaandorpPartial scenario
        :param file_dict:
        :param lon_min:
        :param lon_max:
        :param lat_min:
        :param lat_max:
        :return:
        """
        domain = lon_min, lon_max, lat_min, lat_max
        output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/{}/'.format(
            settings.SCENARIO_NAME)
        utils.check_direc_exist(output_direc)
        # Get the time axis
        time_list = np.array([], dtype=np.int)
        restart_list = np.array([], dtype=np.int)
        for restart in range(settings.SIM_LENGTH):
            parcels_file = file_dict[0][restart]
            parcels_dataset = Dataset(parcels_file)
            time = parcels_dataset.variables['time'][:, :-1]
            for time_case in np.unique(time):
                if np.nansum(time_case == time) > 10 and time_case != np.nan:
                    time_list = np.append(time_list, time_case)
                    restart_list = np.append(restart_list, restart)

        # Initializing the arrays of the timeseries
        size_dict = {}
        for size_class in range(settings.SIZE_CLASS_NUMBER):
            size_dict[size_class] = np.zeros(time_list.shape, dtype=float)
        beach_state_dict = {'beach': deepcopy(size_dict),
                            'adrift': deepcopy(size_dict),
                            'seabed': deepcopy(size_dict),
                            'removed': deepcopy(size_dict),
                            'total': deepcopy(size_dict),
                            'time': time_list}
        beach_label_dict = {'beach': 1, 'adrift': 0, 'seabed': 3, 'removed': 2}

        utils.print_statement("Start running through the restart and run files")
        # loop through the runs
        for run in range(settings.RUN_RANGE):
            # Loop through the restart files
            for restart in range(settings.SIM_LENGTH):
                # Load the lon, lat, time, beach and weight data
                parcels_file = file_dict[run][restart]
                parcels_dataset = Dataset(parcels_file)
                full_data_dict = {}
                for variable in ['lon', 'lat', 'beach', 'time', 'particle_number', 'size_class']:
                    full_data_dict[variable] = parcels_dataset.variables[variable][:, :-1]

                # Just get the particles within the domain, which we do by setting all values not within the domain to
                # nan. These will therefore not be taken into account in the calculations of total counts/weights
                within_domain = utils.particles_in_domain(domain=domain, lon=full_data_dict['lon'],
                                                          lat=full_data_dict['lat'])
                for variable in full_data_dict.keys():
                    full_data_dict[variable][within_domain == False] = np.nan

                # Now, looping through all the time steps and adding up the total weight/counts of particles within each
                # of the beach state domains at each time step
                for time_index, time_value in enumerate(time_list):
                    if restart_list[time_index] == restart:
                        time_selection = full_data_dict['time'] == time_value
                        if np.nansum(full_data_dict['time'] == time_value) > 0:
                            time_dict = {}
                            for variable in ['beach', 'particle_number', 'size_class']:
                                time_dict[variable] = full_data_dict[variable][time_selection]
                            for beach_state in beach_label_dict.keys():
                                beach_selection = time_dict['beach'] == beach_label_dict[beach_state]
                                for size_class in range(settings.SIZE_CLASS_NUMBER):
                                    size_selection = time_dict['size_class'] == size_class
                                    beach_state_dict[beach_state][size_class][time_index] += np.nansum(
                                        time_dict['particle_number'][(beach_selection) & (size_selection)])
                            for size_class in range(settings.SIZE_CLASS_NUMBER):
                                size_selection = time_dict['size_class'] == size_class
                                beach_state_dict['total'][size_class][time_index] += np.nansum(
                                    time_dict['particle_number'][size_selection])

        # Saving the output
        prefix = 'timeseries'
        output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
        utils.save_obj(output_name, beach_state_dict)
        utils.print_statement("The timeseries has been saved")

else:
    def parcels_to_timeseries(file_dict: dict, lon_min: float = -180, lon_max: float = 180, lat_min: float = -90,
                              lat_max: float = 90):
        domain = lon_min, lon_max, lat_min, lat_max
        output_direc = utils.get_output_directory(server=settings.SERVER) + 'timeseries/{}/'.format(
            settings.SCENARIO_NAME)
        utils.check_direc_exist(output_direc)
        # Get the time axis
        time_list = np.array([], dtype=np.int)
        restart_list = np.array([], dtype=np.int)
        for restart in range(settings.SIM_LENGTH):
            parcels_file = file_dict[0][restart]
            parcels_dataset = Dataset(parcels_file)
            time_list = np.append(time_list, parcels_dataset.variables['time'][0, :-1])
            restart_list = np.append(restart_list, np.ones(parcels_dataset.variables['time'][0, :-1].shape) * restart)

        # Initializing the arrays of the timeseries
        beach_state_dict = {'beach': np.zeros(time_list.shape, dtype=float),
                            'adrift': np.zeros(time_list.shape, dtype=float),
                            'seabed': np.zeros(time_list.shape, dtype=float),
                            'removed': np.zeros(time_list.shape, dtype=float),
                            'total': np.zeros(time_list.shape, dtype=float),
                            'time': time_list}
        beach_label_dict = {'beach': 1, 'adrift': 0, 'seabed': 3, 'removed': 2}

        utils.print_statement("Start running through the restart and run files")
        # loop through the runs
        for run in range(settings.RUN_RANGE):
            # Loop through the restart files
            for restart in range(settings.SIM_LENGTH):
                # Load the lon, lat, time, beach and weight data
                parcels_file = file_dict[run][restart]
                parcels_dataset = Dataset(parcels_file)
                full_data_dict = {}
                for variable in ['lon', 'lat', 'beach', 'time']:
                    full_data_dict[variable] = parcels_dataset.variables[variable][:, :-1]
                if 'weights' in parcels_dataset.variables.keys():
                    full_data_dict['weights'] = parcels_dataset.variables['weights'][:, :-1] * settings.BUOYANT
                else:
                    full_data_dict['weights'] = np.ones(full_data_dict['lon'].shape, dtype=float)

                # Just get the particles within the domain, which we do by setting all values not within the domain to
                # nan. These will therefore not be taken into account in the calculations of total counts/weights
                within_domain = utils.particles_in_domain(domain=domain, lon=full_data_dict['lon'],
                                                          lat=full_data_dict['lat'])
                for variable in full_data_dict.keys():
                    full_data_dict[variable][within_domain == False] = np.nan

                # Now, looping through all the time steps and adding up the total weight/counts of particles within each
                # of the beach state domains at each time step
                for index, time_value in enumerate(time_list):
                    if restart_list[index] == restart:
                        time_selection = full_data_dict['time'] == time_value
                        if np.nansum(full_data_dict['time'] == time_value) > 0:
                            time_dict = {}
                            for variable in ['beach', 'weights']:
                                time_dict[variable] = full_data_dict[variable][time_selection]
                            for beach_state in beach_label_dict.keys():
                                beach_state_dict[beach_state][index] += np.nansum(
                                    time_dict['weights'][time_dict['beach'] == beach_label_dict[beach_state]])
                            beach_state_dict['total'][index] += np.nansum(full_data_dict['time'] == time_value)
        # Saving the output
        prefix = 'timeseries'
        output_name = output_direc + utils.analysis_save_file_name(input_file=file_dict[0][0], prefix=prefix)
        utils.save_obj(output_name, beach_state_dict)
        utils.print_statement("The timeseries has been saved")
