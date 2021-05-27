import os
import settings
import factories.scenario_factory as scenario_factory
import factories.analysis_factory as analysis_factory
import factories.visualization_factory as visualization_factory

def run():
    if settings.SUBMISSION == 'simulation':
        os.system('echo "Running the simulation"')
        scenario_name = settings.SCENARIO_NAME
        stokes = settings.STOKES
        server = settings.SERVER
        scenario = scenario_factory.ScenarioFactory.create_scenario(scenario_name=scenario_name, stokes=stokes,
                                                                    server=server)
        scenario.run()
    elif settings.SUBMISSION == 'analysis':
        os.system('echo "Running the analysis"')
        scenario_name = settings.SCENARIO_NAME
        stokes = settings.STOKES
        server = settings.SERVER
        scenario = scenario_factory.ScenarioFactory.create_scenario(scenario_name=scenario_name, stokes=stokes,
                                                                    server=server)
        file_dict = scenario.return_full_run_directory()
        analysis_factory.AnalysisFactory.create_procedure(file_dict=file_dict,
                                                          concentration=settings.CONCENTRATION,
                                                          vertical_concentration=settings.VERTICAL_CONCENTRATION,
                                                          timeseries=settings.TIMESERIES,
                                                          max_distance=settings.MAX_DISTANCE,
                                                          timeslicing=settings.TIMESLICING,
                                                          statistics=settings.STATISTICS)
    elif settings.SUBMISSION == 'visualization':
        os.system('echo "Generating all visualizations"')
        scenario_name = settings.SCENARIO_NAME
        stokes = settings.STOKES
        server = settings.SERVER
        visualization_factory.VisualizationFactory(scenario=scenario_name, stokes=stokes, server=server).run()


if __name__ == "__main__":
    run()
