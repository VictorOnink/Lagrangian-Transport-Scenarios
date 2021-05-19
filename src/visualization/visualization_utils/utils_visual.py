import settings
import utils


def SizeTransport_load_data(scenario, prefix, data_direc, size,  rho, advection_data='CMEMS_MEDITERRANEAN',
                            shore_time=20, start_year=2010, input='Lebreton'):
    """
    Loading the data we want for SizeTransport analysis output, which will generally just differ in terms of which
    particle size and density
    :param scenario:
    :param prefix:
    :param data_direc:
    :param size:
    :param rho:
    :param advection_data:
    :param shore_time:
    :param start_year:
    :param input:
    :return:
    """
    file_name = scenario._file_names(new=True, advection_data=advection_data, shore_time=shore_time, init_size=size,
                                     init_density=rho, start_year=start_year, input=input, run=settings.RUN,
                                     restart=settings.RESTART)
    return data_direc + utils._analysis_save_file_name(input_file=file_name, prefix=prefix)
