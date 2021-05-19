import settings
import visualization.SizeTransport as SizeTransport
import os

def run(scenario, figure_direc: str):
    """
    So, this is the function where I call all the functions for creating the figures. The figures that I don't want to
    run will just be commented out.
    :param scenario:
    :return:
    """
    os.system('echo "start to get the output name"')
    file_name = scenario._file_names(new=True, advection_data='CMEMS_MEDITERRANEAN', shore_time=20, init_size=500e-5,
                                     init_density=920, start_year=2010, input='Lebreton', run=settings.RUN,
                                     restart=settings.RESTART)
    os.system('everything seems to be working, file is {}'.format(file_name))

