import settings
import utils


def SizeTransport_load_data(scenario, prefix, data_direc, size, rho, tau=settings.SEABED_CRIT, restart=settings.RESTART,
                            advection_data='CMEMS_MEDITERRANEAN', shore_time=26, start_year=2010, input='Lebreton',
                            run=settings.RUN, fixed_resus=settings.FIXED_RESUS, resus_time=settings.RESUS_TIME):
    """
    Loading the data we want for SizeTransport analysis output, which will generally just differ in terms of which
    particle size and rho
    :param resus_time:
    :param fixed_resus:
    :param run:
    :param restart:
    :param tau:
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
    file_name = scenario.file_names(new=True, advection_data=advection_data, shore_time=shore_time, init_size=size,
                                    init_density=rho, start_year=start_year, input=input, run=run,
                                    restart=restart, seabed_crit=tau, fixed_resus=fixed_resus, resus_time=resus_time)
    full_path = data_direc + utils.analysis_save_file_name(input_file=file_name, prefix=prefix)
    return utils.load_obj(full_path)


def FragmentationKaandorp_load_data(scenario, prefix, data_direc, shore_time, lambda_frag, rho,
                                    advection_data='CMEMS_MEDITERRANEAN', dn=settings.DN,
                                    size_class_number=settings.SIZE_CLASS_NUMBER, p_frag=settings.P_FRAG,
                                    start_year=settings.STARTYEAR, input='Lebreton'
                                    ):
    file_name = scenario.file_names(new=True, advection_data=advection_data, shore_time=shore_time, dn=dn,
                                    size_class_number=size_class_number, p_frag=p_frag,
                                    lambda_frag=lambda_frag, density=rho, start_year=start_year, input=input)
    full_path = data_direc + utils.analysis_save_file_name(input_file=file_name, prefix=prefix)
    return utils.load_obj(full_path)


def FragmentationKaandorpPartial_load_data(scenario, prefix, data_direc, shore_time, lambda_frag, rho,
                                           advection_data='CMEMS_MEDITERRANEAN', dn=settings.DN,
                                           size_class_number=settings.SIZE_CLASS_NUMBER, p_frag=settings.P_FRAG,
                                           year=settings.STARTYEAR, input=settings.INPUT, postprocess=False,
                                           ocean_frag=settings.OCEAN_FRAG
                                           ):
    file_name = scenario.file_names(new=True, advection_data=advection_data, shore_time=shore_time, dn=dn,
                                    size_class_number=size_class_number, p_frag=p_frag, ocean_frag=ocean_frag,
                                    lambda_frag=lambda_frag, density=rho, year=year, input=input, postprocess=postprocess)
    full_path = data_direc + utils.analysis_save_file_name(input_file=file_name, prefix=prefix)
    return utils.load_obj(full_path)


def FragmentationKaandorpPartial_fielddata_filename(with_type=True):
    file_name = settings.DATA_INPUT_DIR_SERVERS[settings.SERVER] + 'Field_Data/standardized_field_data'
    if with_type:
        file_name += '.pkl'
    return file_name
