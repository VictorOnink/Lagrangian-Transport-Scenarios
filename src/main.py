import src.settings as set
import src.factories.scenario_factory as scenario_factory


def run():
    scenario_name = set.SCENARIO_NAME
    stokes = set.STOKES
    server = set.SERVER
    scenario = scenario_factory.ScenarioFactory.create_scenario(scenario_name=scenario_name, stokes=stokes, server=server)
    scenario.run()


if __name__ == "__main__":
    run()