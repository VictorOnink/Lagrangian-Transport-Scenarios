import settings
import utils
import numpy as np
import pandas as pd


def FragmentationKaandorpPartial_fielddata(to_overwrite=False):
    file_name = settings.DATA_INPUT_DIR_SERVERS[settings.SERVER] + 'Field_Data/standardized_field_data'
    if not utils.check_file_exist(file_name + '.pkl') or to_overwrite:
        output_dict = {}
        output_dict = Fok2017_standardization(output_dict=output_dict)
        output_dict = Constant2019_standardization(output_dict=output_dict)
        output_dict = Cozar2015_standardization(output_dict=output_dict)
        output_dict = RuizOrejon_standardization(output_dict=output_dict)
        utils.save_obj(filename=file_name, item=output_dict)
        utils.print_statement('The standardized field data has been saved.', to_print=True)


def Fok2017_standardization(output_dict: dict):
    """
    Data for Fok et al. (2017). We normalize the particle counts and mass by the size of the particle, which is taken
    as the midpoint of the size bin
    https://doi.org/10.1016/j.envpol.2016.09.079
    """
    prefix = 'Fok'
    output_dict[prefix] = {}
    # Setting the bin sizes
    bin_edges = np.append(np.array([0.315]), np.arange(1, 11))
    bin_midpoint = 10 ** (.5 * (np.log10(bin_edges)[1:] + np.log10(bin_edges)[:-1]))

    # The non-normalized counts and masses
    sample_counts = np.array([29.2, 27.8, 11.2, 11.3, 5.5, 4.1, 2.2, 4.5, 2.3, 2.0])
    sample_mass = np.array([2.2, 7.5, 10.1, 14.2, 11.0, 7.6, 10.6, 15.2, 9.7, 12.0])

    # Normalizing the pdf of the counts and masses by the size of the particles
    output_dict[prefix]['pdf_counts'] = np.divide(sample_counts, bin_edges[1:] - bin_edges[:-1])
    output_dict[prefix]['pdf_mass'] = np.divide(sample_mass, bin_edges[1:] - bin_edges[:-1])
    output_dict[prefix]['bin_edges'] = bin_edges
    output_dict[prefix]['bin_midpoint'] = bin_midpoint

    return output_dict


def Constant2019_standardization(output_dict: dict):
    """
    Data for Constant et al. (2019), again with pdf's normalized by the particle size
    https://doi.org/10.1016/j.marpolbul.2019.03.032
    """
    data_direc = settings.DATA_INPUT_DIR_SERVERS[settings.SERVER] + 'Field_Data/'
    bin_edges = np.array([0.063, 0.315, 0.5, 1, 2.5, 5, 8])
    bin_midpoint = .5 * (np.log10(bin_edges)[1:] + np.log10(bin_edges)[:-1])

    prefix = 'Constant1'
    output_dict[prefix] = {}
    data_pd = pd.read_csv(data_direc + 'Constant2019_1.csv', header=None)
    output_dict[prefix]['pdf_counts'] = data_pd.loc[0::2, 1].values / (bin_edges[1:] - bin_edges[:-1])
    output_dict[prefix]['bin_edges'] = bin_edges
    output_dict[prefix]['bin_midpoint'] = bin_midpoint

    prefix = 'Constant2'
    output_dict[prefix] = {}
    data_pd = pd.read_csv(data_direc + 'Constant2019_2.csv', header=None)
    output_dict[prefix]['pdf_counts'] = data_pd.loc[0::2, 1].values / (bin_edges[1:] - bin_edges[:-1])
    output_dict[prefix]['bin_edges'] = bin_edges
    output_dict[prefix]['bin_midpoint'] = bin_midpoint

    return output_dict


def Cozar2015_standardization(output_dict: dict):
    """
    Data for Cozar et al. (2015), again with pdf's normalized by the particle size
    https://doi.org/10.1371/journal.pone.0121762
    """
    data_direc = settings.DATA_INPUT_DIR_SERVERS[settings.SERVER] + 'Field_Data/'
    prefix = 'Cozar'
    output_dict[prefix] = {}

    data_pd = pd.read_excel(data_direc + 'Cozar_MedData_SizeSpectra.xls', sheet_name=1)
    data_pd = data_pd.drop([0, 29]).reset_index(drop=True)

    output_dict[prefix]['bin_edges'] = np.append(data_pd['Lower Size (mm)'].values,
                                                 data_pd['Upper Size (mm)'].values[-1])
    output_dict[prefix]['bin_midpoint'] = (10 ** data_pd['log Nominal Size']).values

    pdf_tmp = (10 ** data_pd['MED Log # mm-1']).values
    area = np.trapz(pdf_tmp, output_dict[prefix]['bin_midpoint'])
    output_dict[prefix]['pdf_counts'] = np.divide(pdf_tmp, area)

    return output_dict


def RuizOrejon_standardization(output_dict: dict):
    """
    Data for Ruiz-Orejon et al. (2018), again with pdf's normalized by the particle size
    https://doi.org/10.1016/j.marpolbul.2018.06.010
    """
    prefix = 'RuizOrejon'
    output_dict[prefix] = {}

    output_dict[prefix]['bin_edges'] = np.array([0.33, 0.4, 0.5, 0.7, 1.3, 2.5, 4.0, 7.9, 20., 50., 2000.])
    output_dict[prefix]['bin_midpoint'] = 10 ** (0.5 * (
            np.log10(output_dict[prefix]['bin_edges'][1:]) + np.log10(output_dict[prefix]['bin_edges'][:-1])))
    output_dict[prefix]['pdf_counts'] = np.array([0.14830720281849377, 0.09912213358752645, 0.1724104775115928,
                                                  0.33945798285129647, 0.18343691233511597, 0.04388754453458474,
                                                  0.021324131252479426, 0.0056738747517556904,
                                                  0.0004946212353677522, 0.0000891319875335298, 5.599541949072085e-7])[
                                        1:]
    return output_dict