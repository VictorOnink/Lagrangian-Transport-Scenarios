import os
import settings
from factories.scenario_factory import ScenarioFactory
from factories.analysis_factory import AnalysisFactory
from factories.visualization_factory import VisualizationFactory


def run():
    if settings.SUBMISSION == 'simulation':
        os.system('echo "Running the simulation"')
        ScenarioFactory().create_scenario().run()

    elif settings.SUBMISSION == 'analysis':
        os.system('echo "Running the analysis"')
        AnalysisFactory().run_analysis_procedure()

    elif settings.SUBMISSION == 'visualization':
        os.system('echo "Generating all visualizations"')
        VisualizationFactory().run()


if __name__ == "__main__":
    run()
