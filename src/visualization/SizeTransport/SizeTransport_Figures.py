import settings
import visualization.SizeTransport.SizeTransport_Animation as SizeTransport_Animation
import visualization.SizeTransport.SizeTransport_relative_concentrations as SizeTransport_relative_concentrations
import visualization.SizeTransport.SizeTransport_beach_timeseries as SizeTransport_beach_timeseries
import visualization.SizeTransport.SizeTransport_SeaFloorDepthDistribution as SizeTransport_SeaFloorDepthDistribution
import visualization.General as General
import visualization.SizeTransport.SizeTransport_CumulativeDistance as SizeTransport_CumulativeDistance
import visualization.SizeTransport.SizeTransport_SeparationDistance as SizeTransport_SeparationDistance
import visualization.SizeTransport.SizeTransport_TauConcentration as SizeTransport_TauConcentration
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
    # SizeTransport_Animation.SizeTransport_Animation(figure_direc=figure_direc, scenario=scenario, size_list=size_list,
    #                                                 rho_list=rho_list, tau_list=tau_list)

    # Creating figures showing the relative distribution, averaged over the entire simulation and time-snapshots at the
    # end of each simulation year
    time_select = 2
    # SizeTransport_relative_concentrations.SizeTransport_relative_concentrations(figure_direc=figure_direc,
    #                                                                             scenario=scenario,
    #                                                                             size_list=size_list,
    #                                                                             rho_list=rho_list,
    #                                                                             tau_list=tau_list,
    #                                                                             time_selection=time_select,
    #                                                                             beach_state='adrift')
    # SizeTransport_relative_concentrations.SizeTransport_relative_concentrations(figure_direc=figure_direc,
    #                                                                             scenario=scenario,
    #                                                                             size_list=size_list,
    #                                                                             rho_list=rho_list,
    #                                                                             tau_list=tau_list,
    #                                                                             time_selection=time_select,
    #                                                                             beach_state='beach')

    size_list = np.array([5000, 2500, 1250, 625, 313, 156, 78, 39, 20, 10, 5, 2]) * settings.SIZE_FACTOR
    rho_list = np.ones(size_list.shape, dtype=int) * 920
    tau_list = np.array([0])

    # Comparing the relative distributions with different tau values
    # SizeTransport_TauConcentration.SizeTransport_TauConcentration(scenario=scenario, figure_direc=figure_direc,
    #                                                               size_selection=5000 * settings.SIZE_FACTOR,
    #                                                               tau_list=tau_list, beach_state='adrift')
    # SizeTransport_TauConcentration.SizeTransport_TauConcentration(scenario=scenario, figure_direc=figure_direc,
    #                                                               size_selection=625 * settings.SIZE_FACTOR,
    #                                                               tau_list=tau_list, beach_state='adrift')
    # SizeTransport_TauConcentration.SizeTransport_TauConcentration(scenario=scenario, figure_direc=figure_direc,
    #                                                               size_selection=20 * settings.SIZE_FACTOR,
    #                                                               tau_list=tau_list, beach_state='adrift')
    # SizeTransport_TauConcentration.SizeTransport_TauConcentration(scenario=scenario, figure_direc=figure_direc,
    #                                                               size_selection=2 * settings.SIZE_FACTOR,
    #                                                               tau_list=tau_list, beach_state='adrift')

    # Creating figures of the timeseries of the number of particles that are beached/adrift/seabed
    SizeTransport_beach_timeseries.SizeTransport_beach_timeseries(figure_direc=figure_direc, scenario=scenario,
                                                                  size_list=size_list, rho_list=rho_list,
                                                                  tau_list=tau_list, simulation_years=3,
                                                                  tau_comp=False, without_seabed=False)

    # A histogram indicating at which depths particles end up beaching
    # SizeTransport_SeaFloorDepthDistribution.SizeTransport_SeaFloorDepthDistribution(figure_direc=figure_direc,
    #                                                                                 scenario=scenario,
    #                                                                                 size_list=size_list,
    #                                                                                 rho_list=rho_list, histogram=True)
    # SizeTransport_SeaFloorDepthDistribution.SizeTransport_SeaFloorDepthDistribution(figure_direc=figure_direc,
    #                                                                                 scenario=scenario,
    #                                                                                 size_list=size_list,
    #                                                                                 rho_list=rho_list,
    #                                                                                 cumulative=True)

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
