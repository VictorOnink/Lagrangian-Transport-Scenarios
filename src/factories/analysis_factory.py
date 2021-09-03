import Analysis
import utils


class AnalysisFactory():
    def create_procedure(file_dict: dict, scenario, concentration: bool = False, timeseries: bool = False,
                         max_distance: bool = False, vertical_concentration: bool = False, timeslicing: bool = False,
                         statistics: bool = False, separation_distance: bool = False, size_spectrum: bool = False):
        if concentration:
            utils.print_statement("Calculating the concentration", to_print=True)
            Analysis.parcels_to_concentration(file_dict=file_dict)
        if vertical_concentration:
            utils.print_statement("Calculating the vertical concentration", to_print=True)
            Analysis.parcels_to_vertical_concentration(file_dict=file_dict)
        if timeseries:
            utils.print_statement("Calculating timeseries of beached fractions", to_print=True)
            Analysis.parcels_to_timeseries_sizebins(file_dict=file_dict)
        if max_distance:
            utils.print_statement("Calculating maximum distance from shore along particle trajectories", to_print=True)
            Analysis.parcels_to_max_distance(file_dict=file_dict)
        if timeslicing:
            utils.print_statement("Creating time slices of the simulation", to_print=True)
            Analysis.parcels_to_timeslicing(file_dict=file_dict)
        if statistics:
            utils.print_statement("Calculating basic particle trajectory statistics", to_print=True)
            Analysis.parcels_to_basicstatistics(file_dict=file_dict)
        if separation_distance:
            utils.print_statement("Computing separation distances", to_print=True)
            Analysis.parcels_to_separation_distance(file_dict=file_dict, scenario=scenario)
        if size_spectrum:
            utils.print_statement("Computing size distribution", to_print=True)
            Analysis.parcels_to_sizespectrum(file_dict=file_dict)
