import scenarios.base_scenario as base_scenario
from scenarios.advection_diffusion_only_scenario import AdvectionDiffusionOnly
from scenarios.coastal_proximity import CoastalProximity
from scenarios.stochastic_scenario import Stochastic
from scenarios.shore_dependent_resuspension_scenario import SD_Resuspension
from scenarios.Turrel_Beaching_scenario import Turrell_Resuspension
from scenarios.FragmentationKaandorp import FragmentationKaandorp
from scenarios.FragmentationKaandorpPartial import FragmentationKaandorpPartial
from scenarios.SizeTransport import SizeTransport
from scenarios.BlueCloud import BlueCloud
from settings import SCENARIO_NAME


class ScenarioFactory:
    @staticmethod
    def create_scenario(scenario_name=SCENARIO_NAME) -> base_scenario.BaseScenario:
        if scenario_name == "AdvectionDiffusionOnly":
            return AdvectionDiffusionOnly()
        elif scenario_name == "CoastalProximity":
            return CoastalProximity()
        elif scenario_name == "Stochastic":
            return Stochastic()
        elif scenario_name == 'ShoreDependentResuspension':
            return SD_Resuspension()
        elif scenario_name == 'TurrellResuspension':
            return Turrell_Resuspension()
        elif scenario_name == 'FragmentationKaandorp':
            return FragmentationKaandorp()
        elif scenario_name == 'FragmentationKaandorpPartial':
            return FragmentationKaandorpPartial()
        elif scenario_name == 'SizeTransport':
            return SizeTransport()
        elif scenario_name == 'BlueCloud':
            return BlueCloud()
        else:
            raise ValueError("{} is an invalid model scenario!".format(scenario_name))





