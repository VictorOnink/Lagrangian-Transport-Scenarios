import utils
import visualization.SizeTransport.SizeTransport_Figures as SizeTransport
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_Figures as FragmentationKaandorpPartial
from factories.scenario_factory import ScenarioFactory
import settings


class VisualizationFactory:
    def __init__(self):
        self.server = settings.SERVER
        self.stokes = settings.STOKES
        self.scenario_name = settings.SCENARIO_NAME
        self.figure_direc = self.get_figure_direc()
        self.scenario = ScenarioFactory().create_scenario()

    def get_figure_direc(self):
        return settings.FIGURE_OUTPUT_SERVER[self.server] + '{}/'.format(self.scenario_name)

    def run(self):
        if self.scenario_name is 'SizeTransport':
            SizeTransport.run(scenario=self.scenario, figure_direc=self.figure_direc)
        if self.scenario_name is 'FragmentationKaandorpPartial':
            FragmentationKaandorpPartial.run(scenario=self.scenario, figure_direc=self.figure_direc)
        else:
            utils.print_statement('No figures for {} have been generated yet'.format(self.scenario_name), to_print=True)
