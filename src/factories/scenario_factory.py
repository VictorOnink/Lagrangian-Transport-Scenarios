from src.scenarios.base_scenario import BaseScenario
from src.scenarios.coastal_proximity import CoastalProximity


class ScenarioFactory:
    @classmethod
    def create_scenario(cls, scenario_name, server, stokes) -> BaseScenario:
        if scenario_name == "CoastalProximity":
            return CoastalProximity(server=server, stokes=stokes)




