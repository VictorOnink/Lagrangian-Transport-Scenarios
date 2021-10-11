import settings
import numpy as np
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_Animation as Animation
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_timeseries as timeseries
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_FieldDataComp as FieldDataComp
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_Concentration as concentration
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_vertical_profile as vertical_profile
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_boxmodelcomparison as box_model
from visualization.General.General_input_scenario import General_input_scenario
import Analysis


def run(scenario, figure_direc: str):
    shore_time = 20
    lambda_frag_list = np.array([388, 1000, 10000, 35000, 50000])
    rho = 920

    # Standardizing field data
    # Analysis.FragmentationKaandorpPartial_fielddata(to_overwrite=False)

    # Plotting the input scenario
    # General_input_scenario(scenario=scenario, figure_direc=figure_direc).plot()

    # Animation.FragmentationKaandorpPartial_Animation(scenario=scenario, figure_direc=figure_direc, shore_time=20,
    #                                                  rho=rho, simulation_years=2, ocean_frag=False).animate()

    # timeseries.FragmentationKaandorpPartial_timeseries(scenario=scenario, figure_direc=figure_direc,
    #                                                    shore_time=shore_time, lambda_frag=300, rho=rho,
    #                                                    simulation_length=1, weight='particle_number').plot()
    # timeseries.FragmentationKaandorpPartial_timeseries(scenario=scenario, figure_direc=figure_direc,
    #                                                    shore_time=shore_time, lambda_frag=388, rho=rho,
    #                                                    simulation_length=1, weight='particle_number').plot()
    # timeseries.FragmentationKaandorpPartial_timeseries(scenario=scenario, figure_direc=figure_direc,
    #                                                    shore_time=shore_time, lambda_frag=1000, rho=rho,
    #                                                    simulation_length=1, weight='particle_number').plot()
    # timeseries.FragmentationKaandorpPartial_timeseries(scenario=scenario, figure_direc=figure_direc,
    #                                                    shore_time=shore_time, lambda_frag=10000, rho=rho,
    #                                                    simulation_length=1, weight='particle_number').plot()

    # FieldDataComp.FragmentationKaandorpPartial_FieldDataComp(figure_direc=figure_direc, scenario=scenario,
    #                                                          shore_time=shore_time, lambda_frag_list=lambda_frag_list,
    #                                                          rho=rho, sink=False).plot()
    # FieldDataComp.FragmentationKaandorpPartial_FieldDataComp(figure_direc=figure_direc, scenario=scenario,
    #                                                          shore_time=shore_time, lambda_frag_list=lambda_frag_list,
    #                                                          rho=rho, sink=True).plot()
    box_model.FragmentationKaandorpPartial_boxmodelcomparison(figure_direc=figure_direc, scenario=scenario,
                                                              shore_time=shore_time, lambda_frag=388,
                                                              rho=rho, sink=False, sim_length=2).plot()
    box_model.FragmentationKaandorpPartial_boxmodelcomparison(figure_direc=figure_direc, scenario=scenario,
                                                              shore_time=shore_time, lambda_frag=388,
                                                              rho=rho, sink=True, sim_length=2).plot()

    # vertical_profile.FragmentationKaandorpPartial_vertical_profile(figure_direc=figure_direc, scenario=scenario,
    #                                                                shore_time=shore_time, lambda_frag=38000, rho=rho,
    #                                                                simulation_year=1, weight='particle_number').plot()
    # for beach_state in ['beach']:
    #     for year in [0, 1]:
    #         concentration.FragmentationKaandorpPartial_Concentration(scenario=scenario, figure_direc=figure_direc,
    #                                                                  rho=rho, shore_time=shore_time,
    #                                                                  beach_state=beach_state, simulation_year=year,
    #                                                                  lambda_frag=388, mass=False).plot()
    #         concentration.FragmentationKaandorpPartial_Concentration(scenario=scenario, figure_direc=figure_direc,
    #                                                                  rho=rho, shore_time=shore_time,
    #                                                                  beach_state=beach_state, simulation_year=year,
    #                                                                  lambda_frag=388, mass=True).plot()
