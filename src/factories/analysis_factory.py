from Analysis import parcels_to_concentration

class AnalysisFactory():
    def create_procedure(cls, file_dict: dict, concentration: bool = False):
        if concentration == True:
            parcels_to_concentration.parcels_to_concentration(file_dict=file_dict)
