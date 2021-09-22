import settings
from visualization.SizeTransport.SizeTransport_Animation import SizeTransport_Animation
from visualization.SizeTransport.SizeTransport_relative_concentrations import SizeTransport_relative_concentrations
from visualization.SizeTransport.SizeTransport_beach_timeseries import SizeTransport_beach_timeseries
import visualization.General as General
import visualization.SizeTransport.SizeTransport_CumulativeDistance as SizeTransport_CumulativeDistance
import visualization.SizeTransport.SizeTransport_SeparationDistance as SizeTransport_SeparationDistance
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
    rho_list = np.ones(size_list.shape, dtype=int) * 920
    tau_list = np.ones(size_list.shape, dtype=float) * 0.0

    # Creating a figure of the basin bathymetry
    # General.General_bathymetry(scenario=scenario, figure_direc=figure_direc)

    # Figure of the mean wind speed
    # General.General_average_wind_speed(scenario=scenario, figure_direc=figure_direc)

    # Figure of the input scenario
    # General.General_input_scenario(scenario=scenario, figure_direc=figure_direc)

    # Figure of the seasonal average MLD and wind speed
    # General.General_season_average(scenario=scenario, figure_direc=figure_direc, variable='MLD')
    # General.General_season_average(scenario=scenario, figure_direc=figure_direc, variable='wind')

    # Creating an animation showing how the six different size classes I have simulations for at the moment look like
    # SizeTransport_Animation(scenario=scenario, figure_direc=figure_direc, size_list=size_list, simulation_years=2).animate()

    # Creating figures showing the relative distribution, averaged over the entire simulation and time-snapshots at the
    # end of each simulation year
    # time_select = 2
    # SizeTransport_relative_concentrations(scenario=scenario, figure_direc=figure_direc, size_list=size_list,
    #                                       beach_state='adrift', time_selection=time_select).plot()
    # SizeTransport_relative_concentrations(scenario=scenario, figure_direc=figure_direc, size_list=size_list,
    #                                       beach_state='beach', time_selection=time_select).plot()

    size_list = np.array([5000, 2500, 1250, 625, 313, 156, 78, 39, 20, 10, 5, 2]) * settings.SIZE_FACTOR
    # Creating figures of the timeseries of the number of particles that are beached/adrift/seabed
    print(type(SizeTransport_beach_timeseries))
    timeseries_figure = SizeTransport_beach_timeseries(scenario=scenario, figure_direc=figure_direc, size_list=size_list,
                                   simulation_years=2)
    timeseries_figure.plot()

    # Cumulative plots for the total distance travelled vertically and horizontally, and the max depth reached
    # SizeTransport_CumulativeDistance.SizeTransport_CumulativeDistance(figure_direc=figure_direc,
    #                                                                   scenario=scenario,
    #                                                                   size_list=size_list,
    #                                                                   rho_list=rho_list,
    #                                                                   tau_list=tau_list)

    # Plotting the separation distance
    # for size_selection in size_list:
    #     SizeTransport_SeparationDistance.SizeTransport_SeparationDistance(scenario=scenario, figure_direc=figure_direc,
    #                                                                       size_selection=size_selection,
    #                                                                       rho_selection=920, tau_selection=0.14,
    #                                                                       size_list=size_list)
    pass
