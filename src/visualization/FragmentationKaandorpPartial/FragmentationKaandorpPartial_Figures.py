import settings
import numpy as np
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_SizeSpectrumTime as SizeSpectrumTime
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_Animation as Animation
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_timeseries as timeseries


def run(scenario, figure_direc: str):
    shore_time_list = np.array([1, 20])
    lambda_frag_list = np.array([1, 10, 100, 200, 300, 388])
    rho = 920

    # SizeSpectrumTime.FragmentationKaandorpPartial_SizeSpectrumTime(figure_direc=figure_direc, scenario=scenario,
    #                                                                shore_time=20, lambda_frag_list=lambda_frag_list,
    #                                                                density=rho)

    # Animation.FragmentationKaandorpPartial_Animation(scenario=scenario, figure_direc=figure_direc, shore_time=20,
    #                                                  lambda_frag=388, rho=rho)

    for lambda_frag in lambda_frag_list:
        timeseries.FragmentationKaandorpPartial_timeseries(scenario=scenario, figure_direc=figure_direc, shore_time=20,
                                                           lambda_frag=lambda_frag, rho=rho, simulation_length=1)