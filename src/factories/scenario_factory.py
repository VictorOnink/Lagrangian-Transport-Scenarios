import scenarios.base_scenario as base_scenario
import scenarios.advection_diffusion_only_scenario as AdvDifOnly
import scenarios.coastal_proximity as prox
import scenarios.stochastic_scenario as stochastic
import scenarios.shore_dependent_resuspension_scenario as SD_resuspension
import scenarios.Turrel_Beaching_scenario as Turrell
import scenarios.FragmentationKaandorp as FragmentationKaandorp
import scenarios.SizeTransport as SizeTransport


class ScenarioFactory:
    @staticmethod
    def create_scenario(scenario_name, server, stokes) -> base_scenario.BaseScenario:
        if scenario_name == "AdvectionDiffusionOnly":
            return AdvDifOnly.AdvectionDiffusionOnly(server=server, stokes=stokes)
        elif scenario_name == "CoastalProximity":
            return prox.CoastalProximity(server=server, stokes=stokes)
        elif scenario_name == "Stochastic":
            return stochastic.Stochastic(server=server, stokes=stokes)
        elif scenario_name == 'ShoreDependentResuspension':
            return SD_resuspension.SD_Resuspension(server=server, stokes=stokes)
        elif scenario_name == 'TurrellResuspension':
            return Turrell.Turrell_Resuspension(server=server, stokes=stokes)
        elif scenario_name == 'FragmentationKaandorp':
            return FragmentationKaandorp.FragmentationKaandorp(server=server, stokes=stokes)
        elif scenario_name == 'SizeTransport':
            return SizeTransport.SizeTransport(server=server, stokes=stokes)
        else:
            raise ValueError("invalid model scenario")





