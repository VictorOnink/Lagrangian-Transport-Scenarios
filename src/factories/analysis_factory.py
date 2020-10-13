from Analysis import parcels_to_concentration
import os

class AnalysisFactory():
    def create_procedure(file_dict: dict, concentration: bool = False):
        if concentration == True:
            os.system('echo "Calculating the concentration"')
            parcels_to_concentration.parcels_to_concentration(file_dict=file_dict)
