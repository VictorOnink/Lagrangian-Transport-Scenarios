import settings
import utils
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import advection_scenarios.create_input_files as create_input_files
import visualization.visualization_utils as vUtils
from calendar import monthrange


class FragmentationKaandorpPartial_mass_loss:
    def __init__(self, scenario):
        self.input_file_prefix = create_input_files.create_input_files(prefix=settings.ADVECTION_DATA, repeat_dt=None,
                                                                       grid=np.ones(1), lon=np.ones(1), lat=np.ones(1))
        self.frag_list = [388, 1000, 10000, 35000, 50000]
        self.scenario = scenario
        self.file_dict = self.scenario.return_full_run_directory()
        self.data_direc = utils.get_output_directory() + 'timeseries/FragmentationKaandorpPartial/'
        self.shore_time = settings.SHORE_TIME
        self.size_class_number = settings.SIZE_CLASS_NUMBER
        self.simulation_length = 3

    def run(self):
        # Getting the monthly input
        monthly_input = 0
        for run in range(settings.RUN_RANGE):
            particle_number = np.load(self.input_file_prefix + '{}_run={}.npy'.format('lon', run)).size
            post_dataset = utils.load_obj(self.file_dict['postprocess'][settings.STARTYEAR][settings.STARTMONTH][run][0])
            monthly_input += np.sum(post_dataset['particle_mass'][:particle_number * settings.SIZE_CLASS_NUMBER, 0])

        # The total mass of plastic within the 3 year simulation
        total_input = monthly_input * 12 * 3

        # Convert the mass loss term from loss per timestep (30 seconds) to loss per day
        p_sink = settings.P_SINK * 2 * 60 * 24

        # Next, we load the timeseries that give us the total mass in the simulation at any point in the simulation,
        # first when incorporating the sink term and then without the sink term
        final_mass = dict.fromkeys(self.frag_list)
        for lambda_frag in self.frag_list:
            data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=self.scenario, prefix='timeseries',
                                                                      data_direc=self.data_direc,
                                                                      shore_time=self.shore_time,
                                                                      lambda_frag=lambda_frag,
                                                                      rho=settings.INIT_DENSITY,
                                                                      postprocess=True)
            final_mass[lambda_frag] = {'mass_sink': 0, 'mass': 0}
            for size_class in range(self.size_class_number):
                final_mass[lambda_frag]['mass_sink'] += data_dict['total'][size_class]['particle_mass_sink'][-1]
                final_mass[lambda_frag]['mass'] += data_dict['total'][size_class]['particle_mass'][-1]

        # Calculate the percentage mass losses, first from sinking, then sinking and fragmentation and finally just
        # fragmentation, where we remove the sinking loss from the sinking and fragmentation loss
        frag_sink_mass_loss = dict.fromkeys(self.frag_list)
        frag_mass_loss = dict.fromkeys(self.frag_list)
        for lambda_frag in self.frag_list:
            frag_sink_mass_loss[lambda_frag] = (final_mass[lambda_frag]['mass_sink'] - total_input) / total_input * 100
            frag_mass_loss[lambda_frag] = (final_mass[lambda_frag]['mass'] - total_input) / total_input * 100

        # So, now we print our results
        utils.print_statement("\n", to_print=True)
        utils.print_statement("In 3 years, we release {:.2f} particle mass".format(total_input), to_print=True)
        utils.print_statement("\n", to_print=True)

        for lambda_frag in self.frag_list:
            str_format = lambda_frag, final_mass[lambda_frag]['mass_sink'], frag_sink_mass_loss[lambda_frag]
            statement = "With lambda_frag={}, we are left with {:.2f}, a a {:.2f}% reduction".format(*str_format)
            utils.print_statement(statement, to_print=True)
            statement = "Just fragmentation leads to a {:.2f}% reduction".format(frag_mass_loss[lambda_frag])
            utils.print_statement(statement, to_print=True)






