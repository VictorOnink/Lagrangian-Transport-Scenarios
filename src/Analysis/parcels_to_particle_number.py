import settings as settings
import utils
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar


def parcels_to_particle_number(base_file, output_file):
    # Load the relevant data from the base_file, which we take as lambda_frag=388. The split event is independent of
    # lambda_frag, and so it doesn't matter what is taken as the base to calculate the
    base_dataset = Dataset(base_file)
    var_list = ['parent', 'to_split', 'time']
    base_dict = dict.fromkeys(var_list)
    for variable in var_list:
        base_dict[variable] = base_dataset.variables[variable][:, :-1]
    particle_number = base_dict[variable].shape[0]

    # Creating an output dictionary containing an array for the particle number
    output_dict = {'particle_number': np.ones(base_dict[variable].shape, dtype=np.float32)}

    # Computing the fragmentation variable f
    f = 60 / settings.LAMBDA_FRAG

    # Looping through all the particles, identifying the split points
    pbar = ProgressBar()
    for p_id in pbar(range(particle_number)):
        # Getting the index of all cases where a particle is splitting, where we first check if there is a splitting
        # event
        if base_dict['to_split'][p_id, :].max() == 1:
            split_cases = np.where(base_dict['to_split'][p_id, :] == 1)[0]
            split_times = base_dict['time'][p_id, :][split_cases]
            # Looping through all the split cases
            for t_ind in split_cases:
                # calculating the particle number after each splitting event
                number_fraction = utils.particle_number_per_size_class(k=0, f=f)
                output_dict['particle_number'][p_id, (t_ind + 1):] = output_dict['particle_number'][p_id, t_ind] * number_fraction
                # Getting the ID of all the new created particles, where we need to add the +1 to the time index since
                # the new particles are only technically present in the next time step
                c_id, _ = np.where((base_dict['time'] == base_dict['time'][p_id, t_ind + 1]) & (base_dict['parent'] == p_id))
                # Looping through the new particles
                for index_id in range(1, c_id.size):
                    new_particle_number = utils.particle_number_per_size_class(k=index_id, f=f)
                    output_dict['particle_number'][c_id[index_id], :] = output_dict['particle_number'][p_id, t_ind] * new_particle_number

    # Saving the output
    utils.save_obj(output_file, output_dict)

