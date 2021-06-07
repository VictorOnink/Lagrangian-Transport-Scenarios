import settings
import visualization.SizeTransport.SizeTransport_Animation as SizeTransport_Animation
import visualization.SizeTransport.SizeTransport_relative_concentrations as SizeTransport_relative_concentrations
import visualization.SizeTransport.SizeTransport_beach_timeseries as SizeTransport_beach_timeseries
import visualization.SizeTransport.SizeTransport_SeaFloorDepthDistribution as SizeTransport_SeaFloorDepthDistribution
import visualization.General as General
import visualization.SizeTransport.SizeTransport_CumulativeDistance as SizeTransport_CumulativeDistance
import os
import numpy as np


def run(scenario, figure_direc: str):
    """
    So, this is the function where I call all the functions for creating the figures. The figures that I don't want to
    run will just be commented out.
    :param scenario:
    :return:
    """
    size_list = np.array([5000, 500, 50, 10, 5, 1]) * settings.SIZE_FACTOR
    rho_list = np.ones(size_list.shape, dtype=int) * 920

    # Creating a figure of the basin bathymetry
    # General.General_bathymetry(scenario=scenario, figure_direc=figure_direc)

    # Figure of the mean wind speed
    # General.General_average_wind_speed(scenario=scenario, figure_direc=figure_direc)

    # Creating an animation showing how the six different size classes I have simulations for at the moment look like
    # SizeTransport_Animation.SizeTransport_Animation(figure_direc=figure_direc, scenario=scenario, size_list=size_list,
    #                                                 rho_list=rho_list)

    # Creating figures showing the relative distribution, averaged over the entire simulation and time-snapshots at the
    # end of each simulation year
    # time_select = 2
    # SizeTransport_relative_concentrations.SizeTransport_relative_concentrations(figure_direc=figure_direc,
    #                                                                             scenario=scenario, size_list=size_list,
    #                                                                             rho_list=rho_list,
    #                                                                             time_selection=time_select,
    #                                                                             beach_state='afloat')
    # SizeTransport_relative_concentrations.SizeTransport_relative_concentrations(figure_direc=figure_direc,
    #                                                                             scenario=scenario, size_list=size_list,
    #                                                                             rho_list=rho_list,
    #                                                                             time_selection=time_select,
    #                                                                             beach_state='seabed')
    # SizeTransport_relative_concentrations.SizeTransport_relative_concentrations(figure_direc=figure_direc,
    #                                                                             scenario=scenario, size_list=size_list,
    #                                                                             rho_list=rho_list,
    #                                                                             time_selection=time_select,
    #                                                                             beach_state='beach')
    # SizeTransport_relative_concentrations.SizeTransport_relative_concentrations(figure_direc=figure_direc,
    #                                                                             scenario=scenario, size_list=size_list,
    #                                                                             rho_list=rho_list,
    #                                                                             time_selection=time_select,
    #                                                                             difference=True,
    #                                                                             beach_state='afloat')


    size_list = np.array([5000, 1000, 500, 100, 50, 10, 5, 1]) * settings.SIZE_FACTOR
    rho_list = np.ones(size_list.shape, dtype=int) * 920
    #
    # # Creating figures of the timeseries of the number of particles that are beached/afloat/seabed/removed
    SizeTransport_beach_timeseries.SizeTransport_beach_timeseries(figure_direc=figure_direc, scenario=scenario,
                                                                  size_list=size_list, rho_list=rho_list)


    # A histogram indicating at which depths particles end up beaching
    SizeTransport_SeaFloorDepthDistribution.SizeTransport_SeaFloorDepthDistribution(figure_direc=figure_direc,
                                                                                    scenario=scenario,
                                                                                    size_list=size_list,
                                                                                    rho_list=rho_list, histogram=True)
    SizeTransport_SeaFloorDepthDistribution.SizeTransport_SeaFloorDepthDistribution(figure_direc=figure_direc,
                                                                                    scenario=scenario,
                                                                                    size_list=size_list,
                                                                                    rho_list=rho_list, cumulative=True)

    # Cumulative plots for the total distance travelled vertically and horizontally, and the max depth reached
    SizeTransport_CumulativeDistance.SizeTransport_CumulativeDistance(figure_direc=figure_direc,
                                                                      scenario=scenario,
                                                                      size_list=size_list,
                                                                      rho_list=rho_list, )
    pass
