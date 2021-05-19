import settings
import factories.scenario_factory as scenario_factory
import visualization.SizeTransport.SizeTransport_Figures as SizeTransport
import os


def visualization_main(scenario: str = settings.SCENARIO_NAME):
    server = settings.SERVER

    if scenario == 'SizeTransport':
        # Creating the scenario, so that we can call functions to get scenario-specific parameters such as file names
        # All the simulations are essentially with Stokes drift, so set that parameter here
        stokes = 0
        server = settings.SERVER
        scenario = scenario_factory.ScenarioFactory.create_scenario(scenario_name='SizeTransport', stokes=stokes,
                                                                    server=server)
        # Setting the output direc for all the visualizations
        output_direc = settings.FIGURE_OUTPUT_SERVER[server]
        SizeTransport.run(scenario=scenario)
    else:
        os.system('echo "There is currently no figure generating code for this scenario"')