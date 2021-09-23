import settings as settings
import utils
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar
from copy import deepcopy


def parcels_to_particle_number(base_file, output_file, restart_file):
    # Load the relevant data from the base_file, which we take as lambda_frag=388. The split event is independent of
    # lambda_frag, and so it doesn't matter what is taken as the base to calculate the
    base_dataset = Dataset(base_file)
    var_list = ['parent', 'to_split', 'time', 'age', 'size_class']
    base_dict = dict.fromkeys(var_list)
    for variable in var_list:
        base_dict[variable] = base_dataset.variables[variable][:]
    particle_number = base_dict[variable].shape[0]
    time_step_number = base_dict[variable].shape[1]
    final_t = base_dict[variable].shape[1] - 1

    # Computing the particle mass if we account for the gradual removal of mass
    time_step_age = base_dict['age'] // settings.TIME_STEP.total_seconds()
    mass_removed = np.power(1 - settings.P_SINK, time_step_age)

    # Creating an output dictionary containing an array for the particle mass, and for the particle mass when accounting
    # for the gradual mass loss set by P_SINK.
    # Note: This accounting for mass loss needs to be updated whenever we have a split event
    base = np.ones(base_dict[variable].shape, dtype=np.float32)
    output_dict = {'particle_mass': deepcopy(base),
                   'particle_mass_sink': np.multiply(deepcopy(base), mass_removed)}

    if settings.RESTART > 0:
        assert utils.check_file_exist(restart_file), "The restart file {} doesn't exist".format(restart_file)
        restart_dict = utils.load_obj(restart_file)
        for variable in output_dict.keys():
            variable_restart = restart_dict[variable][:, -1]
            for p_ind, number in enumerate(variable_restart):
                output_dict[variable][p_ind, :] = number

    # Computing the fragmentation variable f
    f = 60 / settings.LAMBDA_FRAG

    # Looping through all the particles, identifying the split points
    pbar = ProgressBar()
    for p_id in pbar(range(particle_number)):
        # In the case of including the gradual sink, we multiply the initial values in the array by the gradual
        # Getting the index of all cases where a particle is splitting, where we first check if there is a splitting
        # event
        if base_dict['to_split'][p_id, :].max() == 1:
            split_cases = np.where(base_dict['to_split'][p_id, :] == 1)[0]
            # Looping through all the split cases
            for t_ind in split_cases:
                # calculating the particle number after each splitting event
                mass_fraction = utils.mass_per_size_class(k=0, f=f)
                if t_ind == final_t:
                    output_dict['particle_mass'][p_id, -1] = output_dict['particle_mass'][p_id, t_ind] * mass_fraction
                    output_dict['particle_mass_sink'][p_id, -1] = output_dict['particle_mass_sink'][p_id, t_ind] * mass_fraction
                    # Getting the ID of all the new created particles, that would not have been created at the next time
                    # step since there isn't a next step in the simulation.
                    c_id, _ = np.where((base_dict['time'] == base_dict['time'][p_id, t_ind]) & (base_dict['parent'] == p_id))
                else:
                    output_dict['particle_mass'][p_id, (t_ind + 1):] = output_dict['particle_mass'][p_id, t_ind] * mass_fraction
                    # For the particle_mass_sink, we first account for the loss of mass due to the fragmentation event,
                    # and then we correct for the continued mass loss over time due to P_SINK
                    remaining_time_steps = np.arange(0, time_step_number - (t_ind + 1))
                    sink_correction = np.power(1 - settings.P_SINK, remaining_time_steps)
                    output_dict['particle_mass_sink'][p_id, (t_ind + 1):] = output_dict['particle_mass_sink'][p_id, t_ind] * mass_fraction
                    # output_dict['particle_mass_sink'][p_id, (t_ind + 1):] *= sink_correction
                    # Getting the ID of all the new created particles, where we need to add the +1 to the time index since
                    # the new particles are only technically present in the next time step
                    c_id, _ = np.where((base_dict['time'] == base_dict['time'][p_id, t_ind + 1]) & (base_dict['parent'] == p_id))
                if c_id.size > 0:
                    # Looping through the new particles
                    for index_id in range(1, c_id.size):
                        new_particle_mass = utils.particle_number_per_size_class(k=index_id, f=f)
                        output_dict['particle_mass'][c_id[index_id], :] = output_dict['particle_mass'][p_id, t_ind] * new_particle_mass
                        output_dict['particle_mass_sink'][c_id[index_id], :] = output_dict['particle_mass_sink'][p_id, t_ind] * new_particle_mass

    # Calculating the particle number from the particle masses
    output_dict['particle_number'] = np.multiply(output_dict['particle_mass'], np.power(2, settings.DN * base_dict['size_class']))
    output_dict['particle_number_sink'] = np.multiply(output_dict['particle_mass_sink'], np.power(2, settings.DN * base_dict['size_class']))
    # Saving the output
    utils.save_obj(output_file, output_dict)

    for i in range(particle_number):
        str_format = i, output_dict['particle_mass'][i, -1], output_dict['particle_mass_sink'][i, -1], output_dict['particle_number'][i, -1], output_dict['particle_number_sink'][i, -1]
        print_statement = 'id {}, mass {} mass_sink {}, count {} count_sink {}'.format(*str_format)
        utils.print_statement(print_statement, to_print=True)
