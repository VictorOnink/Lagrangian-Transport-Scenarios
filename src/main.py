import src.settings as set
from src.factories.scenario_factory import ScenarioFactory


def run():
    scenario_name = set.SCENARIO_NAME
    stokes = set.STOKES
    server = set.SERVER
    scenario = ScenarioFactory.create_scenario(scenario_name=scenario_name, stokes=stokes, server=server)
    fieldset = scenario.create_fieldset()
    pset = scenario.create_particle(fieldset)
    scenario.run(pset)


if __name__ == "__main__":
    run()