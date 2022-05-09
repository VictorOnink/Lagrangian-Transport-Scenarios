import settings
import numpy as np
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_Animation as Animation
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_timeseries as timeseries
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_FieldDataComp as FieldDataComp
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_Concentration as concentration
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_vertical_profile as vertical_profile
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_boxmodelcomparison as box_model
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_boxmodel_ocean as box_model_ocean
import visualization.FragmentationKaandorpPartial.FragmentationKaandorpPartial_lonlat_averages as lonlat_average
from visualization.General.General_input_scenario import General_input_scenario
import Analysis
import utils


def run(scenario, figure_direc: str):
    shore_time = 26
    lambda_frag_list = np.array([388, 1000, 10000, 35000, 50000])
    rho = 920

    # Standardizing field data
    # Analysis.FragmentationKaandorpPartial_fielddata(to_overwrite=False)

    # Calculating mass losses
    # Analysis.FragmentationKaandorpPartial_mass_loss(scenario=scenario).run()

    # Plotting the input scenario
    # General_input_scenario(scenario=scenario, figure_direc=figure_direc).plot()

    # Creating an animation showing the horizontal/vertical transport in each size class
    # Animation.FragmentationKaandorpPartial_Animation(scenario=scenario, figure_direc=figure_direc, shore_time=20,
    #                                                  rho=rho, simulation_years=3, ocean_frag=False).animate()

    # Comparing the calculated size distributions from the Lagrangian simulations with field data
    # FieldDataComp.FragmentationKaandorpPartial_FieldDataComp(figure_direc=figure_direc, scenario=scenario,
    #                                                          shore_time=shore_time, lambda_frag_list=lambda_frag_list,
    #                                                          rho=rho, input_list=['LebretonKaandorpInit']).plot()

    # Plot the results from the box model
    # box_model_ocean.FragmentationKaandorpPartial_boxmodel_ocean(figure_direc=figure_direc, size_class_number=12).plot()

    # plot the vertical concentration profiles for all the size classes for each year (seasonally averaged)
    # for year in [0, 1, 2]:
    #     for weight in ['particle_number_sink', 'particle_mass_sink']:
    #         for input in ['LebretonKaandorpInit']:
    #             vertical_profile.FragmentationKaandorpPartial_vertical_profile(figure_direc=figure_direc,
    #                                                                            scenario=scenario, shore_time=shore_time,
    #                                                                            lambda_frag=388, rho=rho, input=input,
    #                                                                            simulation_year=year,
    #                                                                            weight=weight).plot()
    #             vertical_profile.FragmentationKaandorpPartial_vertical_profile(figure_direc=figure_direc,
    #                                                                            scenario=scenario, shore_time=shore_time,
    #                                                                            lambda_frag=388, rho=rho, input=input,
    #                                                                            simulation_year=year,
    #                                                                            weight=weight).plot()

    # Plot the horizontal distribution of the particles in each size class
    for beach_state in ['adrift']:
        for year in [0, 1, 2]:
            for input in ['LebretonKaandorpInit']:
                concentration.FragmentationKaandorpPartial_Concentration(scenario=scenario, figure_direc=figure_direc,
                                                                         rho=rho, shore_time=shore_time,
                                                                         beach_state=beach_state, simulation_year=year,
                                                                         lambda_frag=388, mass=False, input=input).plot()
                concentration.FragmentationKaandorpPartial_Concentration(scenario=scenario, figure_direc=figure_direc,
                                                                         rho=rho, shore_time=shore_time,
                                                                         beach_state=beach_state, simulation_year=year,
                                                                         lambda_frag=388, mass=True, input=input).plot()

    # Calculating some basic statistics about the Lagrangian runs
    # for lambda_frag in lambda_frag_list:
    #     base_statistics = Analysis.FragmentationKaandorpPartial_Statistics(scenario, lambda_frag)
    #     base_statistics.fraction_below_depth(reference_depth=10)
    #     base_statistics.mass_number_fraction_per_size_class()

    utils.print_statement("That is all folks!", to_print=True)

