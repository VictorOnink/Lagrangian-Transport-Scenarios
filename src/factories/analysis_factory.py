from Analysis import parcels_to_concentration, parcels_to_timeseries
import os

class AnalysisFactory():
    def create_procedure(file_dict: dict, concentration: bool = False, timeseries: bool = False):
        if concentration:
            os.system('echo "Calculating the concentration"')
            parcels_to_concentration.parcels_to_concentration(file_dict=file_dict)
        if timeseries:
            os.system('echo "Calculating timeseries of beached fractions"')
            parcels_to_timeseries.parcels_to_timeseries(file_dict=file_dict)
