import settings
import utils
import glob
from advection_scenarios import advection_files
import numpy as np
from netCDF4 import Dataset


class CMEMS_mediterranean_mean_MLD:
    def __init__(self):
        self.MLD_variable = 'mlotst'
        self.start_year = 2010
        self.simulation_years = 3
        self.data_dir = utils.get_data_directory(server=settings.SERVER) + 'CMEMS_MED/'
        self.seasons = ['JFM', 'AMJ', 'JAS', 'OND']

    def create_output_dict(self):
        advection_scenario = advection_files.AdvectionFiles(server=settings.SERVER, stokes=settings.STOKES,
                                                            advection_scenario=settings.ADVECTION_DATA,
                                                            repeat_dt=None).file_names
        GRID, GRID_MASK = advection_scenario['GRID'], advection_scenario['GRID'].mask
        LON, LAT = advection_scenario['LON'], advection_scenario['LAT']
        output_dict = {'LON': LON, "LAT": LAT}
        for year in range(self.start_year, self.start_year + self.simulation_years):
            output_dict[year] = {}
            for season in self.seasons:
                output_dict[year][season] = np.zeros(GRID.shape, dtype=float)
                output_dict[year][season][GRID_MASK] = np.nan
        return output_dict

    def calculate_seasonal_mean(self):
        file_name = self.data_dir + '2010-2012_seasonal_mean_MLD'
        if not utils.check_file_exist(File=file_name, without_pkl=True):
            # Create the output dictionary
            output_dict = self.create_output_dict()

            # Loop through the simulation years and months
            for year in range(self.start_year, self.start_year + self.simulation_years):
                for season, months in zip(self.seasons, [(1, 2, 3), (4, 5, 6), (7, 8, 9), (10, 11, 12)]):
                    # Get all files for the particular season in question
                    file_list = []
                    for m in months:
                        file_list += glob.glob(self.data_dir + str(year) + str(m).zfill(2) + '*AMXL*.nc')

                    # Loop through all the files
                    for file in file_list:
                        output_dict[year][season] += Dataset(file).variables[self.MLD_variable][0, :, :]

                    # Divide by the number of files to get the average
                    output_dict[year][season] /= file_list.__len__()

            utils.save_obj(filename=file_name, item=output_dict)

            return output_dict
        else:
            return utils.load_obj(filename=file_name)


