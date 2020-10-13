import settings as set
import factories.scenario_factory as scenario_factory
import factories.analysis_factory as analysis_factory
import os


def run():
    if set.SUBMISSION == 'simulation':
        os.system('echo "Running the simulation"')
        scenario_name = set.SCENARIO_NAME
        stokes = set.STOKES
        server = set.SERVER
        scenario = scenario_factory.ScenarioFactory.create_scenario(scenario_name=scenario_name, stokes=stokes, server=server)
        scenario.run()
    elif set.SUBMISSION == 'analysis':
        os.system('echo "Running the analysis"')
        scenario_name = set.SCENARIO_NAME
        stokes = set.STOKES
        server = set.SERVER
        scenario = scenario_factory.ScenarioFactory.create_scenario(scenario_name=scenario_name, stokes=stokes, server=server)
        file_dict = scenario.return_full_run_directory()
        analysis_factory.AnalysisFactory.create_procedure(file_dict=file_dict, concentration=set.CONCENTRATION)




if __name__ == "__main__":
    run()