import Analysis
import os

class AnalysisFactory():
    def create_procedure(file_dict: dict, concentration: bool = False, timeseries: bool = False,
                         max_distance: bool = False, vertical_concentration: bool = False):
        if concentration:
            os.system('echo "Calculating the concentration"')
            Analysis.parcels_to_concentration(file_dict=file_dict)
        if vertical_concentration:
            os.system('echo "Calculating the vertical concentration"')
            Analysis.parcels_to_vertical_concentration(file_dict=file_dict)
        if timeseries:
            os.system('echo "Calculating timeseries of beached fractions"')
            Analysis.parcels_to_timeseries(file_dict=file_dict)
        if max_distance:
            os.system('echo "Calculating maximum distance from shore along particle trajectories"')
            Analysis.parcels_to_max_distance(file_dict=file_dict)
