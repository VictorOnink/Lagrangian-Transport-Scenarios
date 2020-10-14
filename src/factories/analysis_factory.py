from Analysis import parcels_to_concentration, parcels_to_timeseries, parcels_to_max_distance
import os

class AnalysisFactory():
    def create_procedure(file_dict: dict, concentration: bool = False, timeseries: bool = False,
                         max_distance: bool = False):
        if concentration:
            os.system('echo "Calculating the concentration"')
            parcels_to_concentration.parcels_to_concentration(file_dict=file_dict)
        if timeseries:
            os.system('echo "Calculating timeseries of beached fractions"')
            parcels_to_timeseries.parcels_to_timeseries(file_dict=file_dict)
        if max_distance:
            os.system('echo "Calculating maximum distance from shore along particle trajectories"')
            parcels_to_max_distance.parcels_to_max_distance(file_dict=file_dict)
