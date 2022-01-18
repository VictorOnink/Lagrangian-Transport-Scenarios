import Analysis
import utils
import settings
from factories.scenario_factory import ScenarioFactory


class AnalysisFactory:
    def __init__(self):
        self.scenario = ScenarioFactory.create_scenario()
        self.file_dict = self.scenario.return_full_run_directory()

    def run_analysis_procedure(self):
        """
        Depending on the parameters in the src/settings.py file, we run specified analysis procedures on the parcels
        output.
        - settings.CONCENTRATION: Calculate the horizontal particle concentrations
        - settings.VERTICAL_CONCENTRATION: Calculate the vertical particle concentrations
        - settings.LONLAT_CONCENTRATION: Calculate the lon/lat averaged horizontal particle concentrations
        - settings.TIMESERIES: Calculate the number of particles in beach/adrift/seabed reservoirs
        - settings.MAX_DISTANCE: Calculate the maximum distance a particle is removed from shore along its trajectory
        - settings.TIMESLICING: Cut up the parcels output so that there is a slice for each output timestep containing
                                all particle positions at that timepoint.
        - settings.STATISTICS: Calculate a number of basic particle statistics
        - settings.SEPARATION: Calculate the separation distances between particles released at same release point
        - settings.SIZE_SPECTRUM: Calculate the number of particle in each size class.
        :return:
        """
        if settings.CONCENTRATION:
            utils.print_statement("Calculating the concentration", to_print=True)
            Analysis.parcels_to_concentration(file_dict=self.file_dict).run()
        if settings.VERTICAL_CONCENTRATION:
            utils.print_statement("Calculating the vertical concentration", to_print=True)
            Analysis.parcels_to_vertical_concentration(file_dict=self.file_dict).run()
        if settings.LONLAT_CONCENTRATION:
            utils.print_statement("Calculating the lonlat concentration averages", to_print=True)
            Analysis.parcels_to_lonlat_average(file_dict=self.file_dict).run()
        if settings.TIMESERIES:
            utils.print_statement("Calculating timeseries of beached fractions", to_print=True)
            Analysis.parcels_to_timeseries(file_dict=self.file_dict).run()
        if settings.MAX_DISTANCE:
            utils.print_statement("Calculating maximum distance from shore along particle trajectories", to_print=True)
            Analysis.parcels_to_max_distance(file_dict=self.file_dict).run()
        if settings.TIMESLICING:
            utils.print_statement("Creating time slices of the simulation", to_print=True)
            Analysis.parcels_to_timeslicing(file_dict=self.file_dict).run()
        if settings.STATISTICS:
            utils.print_statement("Calculating basic particle trajectory statistics", to_print=True)
            Analysis.parcels_to_basicstatistics(file_dict=self.file_dict).run()
        if settings.SEPARATION:
            utils.print_statement("Computing separation distances", to_print=True)
            utils.print_statement("This is currently not rewritten for 2-stage parallel analysis", to_print=True)
            # Analysis.parcels_to_separation_distance(file_dict=self.file_dict, scenario=self.scenario)
        if settings.SIZE_SPECTRUM:
            utils.print_statement("Computing size distribution", to_print=True)
            Analysis.parcels_to_sizespectrum(file_dict=self.file_dict, scenario=self.scenario).run()
