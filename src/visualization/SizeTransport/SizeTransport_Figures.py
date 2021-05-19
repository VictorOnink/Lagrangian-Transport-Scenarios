import settings
import visualization.SizeTransport.SizeTransport_Animation as SizeTransport_Animation
import os

def run(scenario, figure_direc: str):
    """
    So, this is the function where I call all the functions for creating the figures. The figures that I don't want to
    run will just be commented out.
    :param scenario:
    :return:
    """
    # Creating an animation showing how the six different size classes I have simulations for at the moment look like
    SizeTransport_Animation.SizeTransport_Animation(figure_direc=figure_direc)
