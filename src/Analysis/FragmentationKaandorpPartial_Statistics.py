import settings
import utils
import pandas as pd
import numpy as np
import visualization.visualization_utils as vUtils


class FragmentationKaandorpPartial_Statistics:
    """
    This class will have a number of basic statistical tests, such as fractions below any given depth
    """
    def __init__(self, scenario, lambda_frag, shore_time=26, rho=920):
        # Scenario specific variables
        self.scenario = scenario
        self.lambda_frag = lambda_frag
        self.rho = rho
        self.shore_time = shore_time
        # Data specific variables
        self.base_data_direc = utils.get_output_directory(server=settings.SERVER)
        # Other useful variables
        self.seasons = ['JFM', 'AMJ', 'JAS', 'OND']
        self.total_number = 85196

    def fraction_below_depth(self, reference_depth):
        """
        Calculate the fraction of particles found below a given reference depth at the end of the simulation
        :param reference_depth:
        :param shore:
        :param adrift:
        :return:
        """
        prefix = 'vertical_concentration'
        data_direc = self.base_data_direc + 'concentrations/FragmentationKaandorpPartial/'
        conc_types = ['concentration_mass', 'concentration_number']

        # Loading the profile in question
        data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=self.scenario, prefix=prefix,
                                                                  data_direc=data_direc,
                                                                  shore_time=self.shore_time,
                                                                  lambda_frag=self.lambda_frag,
                                                                  rho=self.rho, postprocess=True)
        depth_bins = data_dict['depth']

        # Getting the last year and month
        data_dict = data_dict[utils.analysis_simulation_year_key(2)][11]

        # For each season, calculating the fraction of particles below reference_depth (note greater depths are higher
        # values!)
        below_dict = dict.fromkeys(range(settings.SIZE_CLASS_NUMBER))
        depth_selection = depth_bins > reference_depth

        for size_class in range(settings.SIZE_CLASS_NUMBER):
            below_dict[size_class] = {}
            for conc in conc_types:
                below_dict[size_class][conc] = np.nansum(data_dict[size_class][conc][depth_selection]) / np.nansum(data_dict[size_class][conc]) * 100

        # Print the percentage of particles below a given depth
        print('For lambda_frag {}, we have the following fractions below {} m'.format(self.lambda_frag, reference_depth))
        for size_class in range(settings.SIZE_CLASS_NUMBER):
            print("\tk = {} mass = {:.2f}%, counts = {:.2f}".format(size_class,
                                                                    below_dict[size_class]['concentration_mass'],
                                                                    below_dict[size_class]['concentration_number']))

    def mass_number_fraction_per_size_class(self):
        prefix = 'vertical_concentration'
        data_direc = self.base_data_direc + 'concentrations/FragmentationKaandorpPartial/'
        conc_types = ['concentration_mass', 'concentration_number']

        # Loading the profile in question
        data_dict = vUtils.FragmentationKaandorpPartial_load_data(scenario=self.scenario, prefix=prefix,
                                                                  data_direc=data_direc,
                                                                  shore_time=self.shore_time,
                                                                  lambda_frag=self.lambda_frag,
                                                                  rho=self.rho, postprocess=True,
                                                                  input=self.input)

        # Getting the last year and month
        data_dict = data_dict[utils.analysis_simulation_year_key(2)][11]

        # Calculating the total counts and mass per size class throughout the entire water column
        total_dict = dict.fromkeys(range(settings.SIZE_CLASS_NUMBER))
        for size_class in range(settings.SIZE_CLASS_NUMBER):
            total_dict[size_class] = {}
            for conc in conc_types:
                total_dict[size_class][conc] = np.nansum(data_dict[size_class][conc])

        # Calculate the total counts and mass
        total_mass, total_counts = 0, 0
        for size_class in range(settings.SIZE_CLASS_NUMBER):
            total_mass += total_dict[size_class]['concentration_mass']
            total_counts += total_dict[size_class]['concentration_number']

        # Print the percentage of particles below a given depth
        print('For lambda_frag {}, we have the following mass/number distribution'.format(self.lambda_frag))
        for size_class in range(settings.SIZE_CLASS_NUMBER):
            print("\tk = {} mass = {:.2f}%, counts = {:.2f}".format(size_class,
                                                                    total_dict[size_class]['concentration_mass'] / total_mass * 100,
                                                                    total_dict[size_class]['concentration_number'] / total_counts * 100))






