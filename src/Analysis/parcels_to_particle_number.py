import settings as settings
import utils
from netCDF4 import Dataset
import numpy as np
from progressbar import ProgressBar
from copy import deepcopy


class parcels_to_particle_number:
    def __init__(self, base_file, output_file, restart_file):
        # The names of the parcels base_file, the output file of the analysis and restart particle number file for when
        # restart > 0
        self.base_file, self.output_file, self.restart_file = base_file, output_file, restart_file
        self.restart = settings.RESTART
        # Model parameters
        self.dt = settings.TIME_STEP.total_seconds()
        self.f = 60 / settings.LAMBDA_FRAG
        self.DN = settings.DN
        # Analysis parameters
        self.var_list = ['parent', 'to_split', 'time', 'age', 'size_class']
        self.mass_list = ['particle_mass', 'particle_mass_sink']
        self.number_list = ['particle_number', 'particle_number_sink']
        self.mass_number_list = utils.flatten_list_of_lists([self.mass_list, self.number_list])
        self.output_dict = {}
        self.particle_number, self.time_step_number, self.final_t = None, None, None

    def restart_initialization(self, output_dict):
        if self.restart > 0:
            assert utils.check_file_exist(self.restart_file), "Restart file {} doesn't exist".format(self.restart_file)
            restart_dict = utils.load_obj(self.restart_file)
            for variable in self.mass_list:
                variable_restart = restart_dict[variable][:, -1]
                for p_ind, number in enumerate(variable_restart):
                    output_dict[variable][p_ind, :] = number
        return output_dict

    def run(self):
        # Loading the variables from the parcels base file (set to always have lambda_frag=388)
        base_dataset = Dataset(self.base_file)
        base_dict = dict.fromkeys(self.var_list)
        for variable in base_dict.keys():
            base_dict[variable] = base_dataset.variables[variable][:]

        # The number of particles, the number of time steps and the final time index
        self.particle_number, self.time_step_number = base_dict[variable].shape
        self.final_t = self.time_step_number - 1

        # Computing the mass removal fraction
        age_in_timesteps = base_dict['age'] // self.dt
        mass_remainder = np.power(1 - settings.P_SINK, age_in_timesteps)

        # Initializing the output dict with arrays equal in size to the variables in self.base_file
        # Note: The mass loss accounting needs to be updated whenever we have a split event
        base_array = np.ones(base_dict[variable].shape, dtype=np.float32)
        for variable in self.mass_number_list:
            self.output_dict[variable] = deepcopy(base_array)

        # Initializing the particle numbers for when restart > 0
        self.output_dict = self.restart_initialization(output_dict=self.output_dict)

        # Adding the age correction to the particle mass array
        self.output_dict['particle_mass_sink'] = np.multiply(self.output_dict['particle_mass_sink'], mass_remainder)

        # Looping through the particles, identifying the split points
        pbar = ProgressBar()
        for p_id in pbar(range(self.particle_number)):
            # Checking if the particle has a split event, and then determining the split events
            if base_dict['to_split'][p_id, :].max() == 1:
                split_cases = np.where(base_dict['to_split'][p_id, :] == 1)[0]
                # We make an array that tracks the particles created in previous split events
                previous_split = [p_id]
                # Looping through the split events, calculating the particle mass after each split event
                for t_ind in split_cases:
                    mass_fraction = utils.mass_per_size_class(k=0, f=self.f)
                    if t_ind == self.final_t:
                        for variable in self.mass_list:
                            self.output_dict[variable][p_id, -1] = self.output_dict[variable][p_id, t_ind] * mass_fraction
                        # Getting the ID of all the new created particles at self.final_t
                        c_id, _ = np.where((base_dict['time'] == base_dict['time'][p_id, t_ind]) &
                                           (base_dict['parent'] == p_id))
                    else:
                        for variable in self.mass_list:
                            self.output_dict[variable][p_id, (t_ind + 1):] = self.output_dict[variable][p_id, t_ind] * mass_fraction
                        # For particle_mass_sink, we account for the loss of mass due to the P_sink term after the
                        # fragmentation event
                        sink_correction = np.power(1 - settings.P_SINK, np.arange(0, self.time_step_number - (t_ind + 1)))
                        self.output_dict['particle_mass_sink'][p_id, (t_ind + 1):] *= sink_correction
                        # Getting the ID of all the new created particles, where we need to add the +1 to the time index
                        # since the new particles are only technically present in the next time step
                        c_id, _ = np.where((base_dict['time'] == base_dict['time'][p_id, t_ind + 1]) & (base_dict['parent'] == p_id))
                    # Removing all particles that were already involved in a previous split so that we only consider
                    # newly created particles
                    for prospect_p in c_id:
                        if prospect_p in previous_split:
                            c_id = np.delete(c_id, np.where(c_id == prospect_p))
                    # Looping through the newly created particles, where the first is skipped as it is the parent
                    if c_id.size > 0:
                        for index_id in range(0, c_id.size):
                            new_particle_mass = utils.mass_per_size_class(k=index_id, f=self.f)
                            if c_id[index_id] not in previous_split:
                                if p_id == 55:
                                    print('c_id[index_id] = {}'.format(c_id[index_id]))
                                    print('size_class = {}'.format(base_dict['size_class'][c_id[index_id], 0]))
                                for variable in self.mass_list:
                                    self.output_dict[variable][c_id[index_id], :] = self.output_dict[variable][p_id, t_ind] * new_particle_mass
                                    if p_id == 55:
                                        print('self.output_dict[variable][c_id[index_id], 0] = {}'.format(self.output_dict[variable][c_id[index_id], 0]))
                                previous_split.append(c_id[index_id])
                            # Accounting again for mass loss
                            self.output_dict['particle_mass_sink'][c_id[index_id], :] *= mass_remainder[c_id[index_id], :]

        for p_id in [255, 256]:
            utils.print_statement('{} {} parent={} {}'.format(p_id,
                                                                 base_dict['size_class'][p_id, 0],
                                                                 base_dict['parent'][p_id, 0],
                                                                 np.unique(self.output_dict['particle_mass'][p_id, :])[::-1]), to_print=True)

        # for p_id in range(self.particle_number):
        #     utils.print_statement('{} {} parent={} {}'.format(p_id,
        #                                                          base_dict['size_class'][p_id, 0],
        #                                                          base_dict['parent'][p_id, 0],
        #                                                          self.output_dict['particle_mass'][p_id, 0], to_print=True))


        # Calculating the particle number from the particle masses
        mass_to_number = np.power(2, self.DN * base_dict['size_class'])
        for number, mass in zip(self.number_list, self.mass_list):
            self.output_dict[number] = np.multiply(self.output_dict[mass], mass_to_number)

        # Masking the output
        for variable in self.mass_number_list:
            self.output_dict[variable] = np.ma.masked_array(self.output_dict[variable], mask=base_dict['time'].mask)

        # Saving the output
        utils.save_obj(self.output_file, self.output_dict)
