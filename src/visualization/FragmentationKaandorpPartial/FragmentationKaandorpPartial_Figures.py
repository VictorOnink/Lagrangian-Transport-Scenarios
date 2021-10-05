import settings
import numpy as np
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_Animation as Animation
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_timeseries as timeseries
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_FieldDataComp as FieldDataComp
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_Concentration as concentration
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_vertical_profile as vertical_profile
import Analysis


def run(scenario, figure_direc: str):
    shore_time = 20
    lambda_frag_list = np.array([388, 1000, 10000, 35000, 50000])
    rho = 920

    # Standardizing field data
    # Analysis.FragmentationKaandorpPartial_fielddata(to_overwrite=False)

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

    # lambda_frag_list = np.array([388, 1000, 10000, 38000])
    FieldDataComp.FragmentationKaandorpPartial_FieldDataComp(figure_direc=figure_direc, scenario=scenario,
                                                             shore_time=shore_time, lambda_frag_list=lambda_frag_list,
                                                             rho=rho, sink=False).plot()

    # vertical_profile.FragmentationKaandorpPartial_vertical_profile(figure_direc=figure_direc, scenario=scenario,
    #                                                                shore_time=shore_time, lambda_frag=38000, rho=rho,
    #                                                                simulation_year=1, weight='particle_number').plot()
