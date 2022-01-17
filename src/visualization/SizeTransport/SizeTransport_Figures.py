import settings
import utils
from visualization.SizeTransport.SizeTransport_Animation import SizeTransport_Animation
from visualization.SizeTransport.SizeTransport_relative_concentrations import SizeTransport_relative_concentrations
from visualization.SizeTransport.SizeTransport_beach_timeseries import SizeTransport_beach_timeseries
import visualization.General as General
from visualization.SizeTransport.SizeTransport_CumulativeDistance import SizeTransport_CumulativeDistance
from visualization.SizeTransport.SizeTransport_SeparationDistance import SizeTransport_SeparationDistance
from visualization.SizeTransport.SizeTransport_VerticalProfile import SizeTransport_VerticalProfile
from visualization.SizeTransport.SizeTransport_lonlat_averages import SizeTransport_lonlat_averages
from visualization.SizeTransport.SizeTransport_reservoirs import SizeTransport_reservoirs
from visualization.SizeTransport.SizeTransport_rho_concentrations import SizeTransport_rho_concentrations
from visualization.SizeTransport.SizeTransport_full_concentrations import SizeTransport_full_concentrations
import os
import numpy as np


def run(scenario, figure_direc: str):
    """
    So, this is the function where I call all the functions for creating the figures. The figures that I don't want to
    run will just be commented out.
    :param scenario:
    :return:
    """
    # 5000 2500 1250 625 313 156 78 39 20 10 5 2
    size_list = np.array([5000, 1250, 313, 78, 20, 2]) * settings.SIZE_FACTOR

    # Creating a figure of the basin bathymetry
    # General.General_bathymetry(scenario=scenario, figure_direc=figure_direc).plot()

    # Figure of the mean wind speed
    # General.General_average_wind_speed(scenario=scenario, figure_direc=figure_direc).plot()

    # Figure of the input scenario
    # General.General_input_scenario(scenario=scenario, figure_direc=figure_direc).plot()

    # Figure of the seasonal average MLD and wind speed
    # General.General_season_average(scenario=scenario, figure_direc=figure_direc, variable='MLD').plot()
    # General.General_season_average(scenario=scenario, figure_direc=figure_direc, variable='wind').plot()

    # Creating an animation showing how the six different size classes I have simulations for at the moment look like
    # SizeTransport_Animation(scenario=scenario, figure_direc=figure_direc, size_list=size_list, simulation_years=2).animate()

    # Creating figures showing the relative distribution, averaged over the entire simulation and time-snapshots at the
    # end of each simulation year
    # for time_select in [0]:
    #     for rho in [30, 920, 980, 1020]:
    #         SizeTransport_relative_concentrations(scenario=scenario, figure_direc=figure_direc, size_list=size_list,
    #                                               beach_state='adrift', time_selection=time_select, rho=rho).plot()
    #         SizeTransport_relative_concentrations(scenario=scenario, figure_direc=figure_direc, size_list=size_list,
    #                                               beach_state='beach', time_selection=time_select, rho=rho).plot()

    # Plotting the relative distributions for fixed particle sizes, but with different particle densities
    # for time_select in [0]:
    #     for size in np.array([5000, 2500, 1250, 625, 313, 156, 78, 39, 20, 10, 5, 2]) * settings.SIZE_FACTOR:
    #         SizeTransport_rho_concentrations(scenario=scenario, figure_direc=figure_direc, size=size,
    #                                          beach_state='adrift', time_selection=time_select,
    #                                          rho_list=[30, 920, 980, 1020]).plot()

    # Plotting all horizontal concentrations for a given density
    # for rho in [30, 920, 980, 1020]:
    #     for time_select in [0]:
    #         for depth_level in ['surface_1m', 'surface_5m', 'column']:
    #             SizeTransport_full_concentrations(scenario=scenario, figure_direc=figure_direc, beach_state='adrift',
    #                                               time_selection=time_select, rho=rho, depth_level=depth_level).plot()

    size_list = np.array([5000, 2500, 1250, 625, 313, 156, 78, 39, 20, 10, 5, 2]) * settings.SIZE_FACTOR
    # Creating figures of the timeseries of the number of particles that are beached/adrift/seabed
    # SizeTransport_beach_timeseries(scenario=scenario, figure_direc=figure_direc, size_list=size_list,
    #                                simulation_years=3, rho_list=[30, 920, 980, 1020]).plot()

    # Figure showing the beached/adrift fractions of each size class
    SizeTransport_reservoirs(scenario=scenario, figure_direc=figure_direc, size_list=size_list).plot()

    # Cumulative plots for the total distance travelled vertically and horizontally, and the max depth reached
    # SizeTransport_CumulativeDistance(figure_direc=figure_direc, scenario=scenario, size_list=size_list).plot()

    # Plotting the month average vertical profile
    # for rho in [[920, 980], [30, 1020]]:
    #     SizeTransport_VerticalProfile(scenario=scenario, figure_direc=figure_direc, size_list=size_list,
    #                                   time_selection=0, rho_list=rho).plot()
    #     SizeTransport_VerticalProfile(scenario=scenario, figure_direc=figure_direc, size_list=size_list,
    #                                   time_selection=1, rho_list=rho).plot()
    #     SizeTransport_VerticalProfile(scenario=scenario, figure_direc=figure_direc, size_list=size_list,
    #                                   time_selection=2, rho_list=rho).plot()

    # Plotting the separation distance
    # for size_selection in size_list:
    #     SizeTransport_SeparationDistance(scenario=scenario, figure_direc=figure_direc, size_selection=size_selection,
    #                                      size_list=size_list).plot()

    # Plotting the lon lat concentration averages
    # for time in [0]:
    #     SizeTransport_lonlat_averages(scenario=scenario, figure_direc=figure_direc, size_list=size_list,
    #                                   time_selection=time, beach_state='beach').plot()
    #     SizeTransport_lonlat_averages(scenario=scenario, figure_direc=figure_direc, size_list=size_list,
    #                                   time_selection=time, beach_state='adrift').plot()

    utils.print_statement("That is all folks!", to_print=True)
