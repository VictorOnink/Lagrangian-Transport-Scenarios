import os
os.system('echo "is this working"')
import settings
import factories.scenario_factory as scenario_factory
import factories.analysis_factory as analysis_factory


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
        analysis_factory.AnalysisFactory.create_procedure(file_dict=file_dict, concentration=settings.CONCENTRATION,
                                                          timeseries=settings.TIMESERIES, max_distance=settings.MAX_DISTANCE)


if __name__ == "__main__":
    run()
