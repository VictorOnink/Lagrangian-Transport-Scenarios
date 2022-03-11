import settings
import utils
import pandas as pd
import numpy as np
import visualization.visualization_utils as vUtils


class SizeTransport_Statistics:
    """
    This class will have a number of basic statistical tests, such as fractions below any given depth
    """
    def __init__(self, scenario, size, rho=920, tau=0.0, fixed_resus=False, resus_time=7):
        # Scenario specific variables
        self.scenario = scenario
        self.size = size
        self.rho = rho
        self.tau = tau
        self.fixed_resus = fixed_resus
        self.resus_time = resus_time
        # Data specific variables
        self.base_data_direc = utils.get_output_directory(server=settings.SERVER)
        # Other useful variables
        self.seasons = ['JFM', 'AMJ', 'JAS', 'OND']
        self.total_number = 85196

    def fraction_below_depth(self, reference_depth, shore='all'):
        prefix = 'vertical_concentration'
        data_direc = self.base_data_direc + 'concentrations/SizeTransport/'
        conc_type = {'all': 'concentration', 'offshore': 'concentration_offshore',
                     'nearshore': 'concentration_nearshore'}[shore]

        # Loading the profile in question
        data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=prefix,
                                                   data_direc=data_direc, fixed_resus=self.fixed_resus,
                                                   size=self.size, rho=self.rho, tau=self.tau,
                                                   resus_time=self.resus_time)
        depth_bins = data_dict['depth']

        # Calculating the seasonal averages
        season_dict = {}
        for season_ind, season in enumerate(self.seasons):
            season_dict[season] = (data_dict[utils.analysis_simulation_year_key(0)][season_ind * 3][conc_type] +
                                   data_dict[utils.analysis_simulation_year_key(1)][season_ind * 3 + 1][conc_type] +
                                   data_dict[utils.analysis_simulation_year_key(2)][season_ind * 3 + 2][conc_type]) / 3

        # For each season, calculating the fraction of particles below reference_depth (note greater depths are higher
        # values!)
        below_dict = {}
        depth_selection = depth_bins > reference_depth
        for season in self.seasons:
            below_dict[season] = np.nansum(season_dict[season][depth_selection]) / self.total_number * 100

        # Calculate the average fraction below a given depth
        below_dict['total'] = 0
        for season in self.seasons:
            below_dict['total'] += below_dict[season] / 4

        # Print the percentage of particles below a given depth
        opening = 'For a particle with size {:.3f} and density {}, we have the following fractions below {} m'.format(self.size, self.rho, reference_depth),
        print(opening)
        for key in below_dict.keys():
            print("\t{} {}%".format(key, below_dict[key]))
        print('\n')

