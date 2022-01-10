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
    Analysis.FragmentationKaandorpPartial_fielddata(to_overwrite=True)

    # Plotting the input scenario
    # General_input_scenario(scenario=scenario, figure_direc=figure_direc).plot()

    # Animation.FragmentationKaandorpPartial_Animation(scenario=scenario, figure_direc=figure_direc, shore_time=20,
    #                                                  rho=rho, simulation_years=3, ocean_frag=False).animate()

    # timeseries.FragmentationKaandorpPartial_timeseries(scenario=scenario, figure_direc=figure_direc,
    #                                                    shore_time=shore_time, lambda_frag=10000, rho=rho,
    #                                                    simulation_length=3, weight='particle_number_sink').plot()
    # timeseries.FragmentationKaandorpPartial_timeseries(scenario=scenario, figure_direc=figure_direc,
    #                                                    shore_time=shore_time, lambda_frag=35000, rho=rho,
    #                                                    simulation_length=3, weight='particle_number_sink').plot()

    FieldDataComp.FragmentationKaandorpPartial_FieldDataComp(figure_direc=figure_direc, scenario=scenario,
                                                             shore_time=shore_time, lambda_frag_list=lambda_frag_list,
                                                             rho=rho, input_list=['LebretonKaandorpInit']).plot()

    box_model.FragmentationKaandorpPartial_boxmodelcomparison(figure_direc=figure_direc, scenario=scenario,
                                                              shore_time=shore_time, lambda_frag=388,
                                                              rho=rho, sim_length=2).plot()

    box_model_ocean.FragmentationKaandorpPartial_boxmodel_ocean(figure_direc=figure_direc, size_class_number=10,
                                                                lambda_frag=388).plot()
    box_model_ocean.FragmentationKaandorpPartial_boxmodel_ocean(figure_direc=figure_direc, size_class_number=10,
                                                                lambda_frag=35000).plot()

    for year in [0, 1, 3]:
        for weight in ['particle_number_sink', 'particle_mass_sink']:
            for input in ['LebretonKaandorpInit']:
                vertical_profile.FragmentationKaandorpPartial_vertical_profile(figure_direc=figure_direc,
                                                                               scenario=scenario, shore_time=shore_time,
                                                                               lambda_frag=388, rho=rho, input=input,
                                                                               simulation_year=year,
                                                                               weight=weight).plot()
                vertical_profile.FragmentationKaandorpPartial_vertical_profile(figure_direc=figure_direc,
                                                                               scenario=scenario, shore_time=shore_time,
                                                                               lambda_frag=388, rho=rho, input=input,
                                                                               simulation_year=year,
                                                                               weight=weight).plot()

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

    for beach_state in ['adrift', 'beach']:
        for year in [0, 1, 2]:
            for input in ['LebretonKaandorpInit']:
                lonlat_average.FragmentationKaandorpPartial_lonlat_averages(scenario=scenario, figure_direc=figure_direc,
                                                                            lambda_frag=35000, time_selection=year,
                                                                            beach_state=beach_state, mass=False,
                                                                            input=input).plot()
                lonlat_average.FragmentationKaandorpPartial_lonlat_averages(scenario=scenario, figure_direc=figure_direc,
                                                                            lambda_frag=35000, time_selection=year,
                                                                            beach_state=beach_state, mass=True,
                                                                            input=input).plot()

    utils.print_statement("That's all folks!", to_print=True)

