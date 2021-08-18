import visualization.SizeTransport.SizeTransport_Figures as SizeTransport
import visualization.FragmentationKaandorp.FragmentationKaandorp_Figures as FragmentationKaandorp
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_Figures as FragmentationKaandorpPartial
import factories.scenario_factory as scenario_factory
import settings


class VisualizationFactory:
    def __init__(self, server, stokes, scenario):
        self.server = server
        self.stokes = stokes
        self.scenario_name = scenario
        self.figure_direc = self.get_figure_direc()
        self.scenario = scenario_factory.ScenarioFactory.create_scenario(scenario_name=self.scenario_name,
                                                                         stokes=self.stokes,
                                                                         server=self.server)

    def get_figure_direc(self):
        return settings.FIGURE_OUTPUT_SERVER[self.server] + '{}/'.format(self.scenario_name)

    def run(self):
        if self.scenario_name is 'SizeTransport':
            SizeTransport.run(scenario=self.scenario, figure_direc=self.figure_direc)
        if self.scenario_name is 'FragmentationKaandorp':
            FragmentationKaandorp.run(scenario=self.scenario, figure_direc=self.figure_direc)
        if self.scenario_name is 'FragmentationKaandorpPartial':
            FragmentationKaandorpPartial.run(scenario=self.scenario, figure_direc=self.figure_direc)
