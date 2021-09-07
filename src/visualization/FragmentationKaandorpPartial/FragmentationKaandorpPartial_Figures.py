import settings
import numpy as np
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_SizeSpectrumTime as SizeSpectrumTime
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_Animation as Animation
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_timeseries as timeseries
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_SizeSpectrumBeach as SizeSpectrumBeach
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_FieldDataComp as FieldDataComp
import Analysis


def run(scenario, figure_direc: str):
    shore_time_list = np.array([1, 20])
    lambda_frag_list = np.array([100, 200, 300, 388, 1000, 10000])
    rho = 920

    # Standardizing field data
    # Analysis.FragmentationKaandorpPartial_fielddata(to_overwrite=False)

    # SizeSpectrumTime.FragmentationKaandorpPartial_SizeSpectrumTime(figure_direc=figure_direc, scenario=scenario,
    #                                                                shore_time=20, lambda_frag_list=lambda_frag_list,
    #                                                                density=rho)

    # Animation.FragmentationKaandorpPartial_Animation(scenario=scenario, figure_direc=figure_direc, shore_time=20,
    #                                                  lambda_frag=388, rho=rho, simulation_years=1, ocean_frag=False)

    # timeseries.FragmentationKaandorpPartial_timeseries(scenario=scenario, figure_direc=figure_direc, shore_time=20,
    #                                                    lambda_frag=300, rho=rho, simulation_length=1)
    # timeseries.FragmentationKaandorpPartial_timeseries(scenario=scenario, figure_direc=figure_direc, shore_time=20,
    #                                                    lambda_frag=388, rho=rho, simulation_length=1)
    # timeseries.FragmentationKaandorpPartial_timeseries(scenario=scenario, figure_direc=figure_direc, shore_time=20,
    #                                                    lambda_frag=1000, rho=rho, simulation_length=1)
    # timeseries.FragmentationKaandorpPartial_timeseries(scenario=scenario, figure_direc=figure_direc, shore_time=20,
    #                                                    lambda_frag=10000, rho=rho, simulation_length=1)


    # SizeSpectrumBeach.FragmentationKaandorpPartial_SizeSpectrumBeach(figure_direc=figure_direc, scenario=scenario,
    #                                                                  shore_time=20, lambda_frag=388,
    #                                                                  density=rho)
    # SizeSpectrumBeach.FragmentationKaandorpPartial_SizeSpectrumBeach(figure_direc=figure_direc, scenario=scenario,
    #                                                                  shore_time=20, lambda_frag=300,
    #                                                                  density=rho)
    # SizeSpectrumBeach.FragmentationKaandorpPartial_SizeSpectrumBeach(figure_direc=figure_direc, scenario=scenario,
    #                                                                  shore_time=20, lambda_frag=10000,
    #                                                                  density=rho)
    # SizeSpectrumBeach.FragmentationKaandorpPartial_SizeSpectrumBeach(figure_direc=figure_direc, scenario=scenario,
    #                                                                  shore_time=20, lambda_frag=1000,
    #                                                                  density=rho)

