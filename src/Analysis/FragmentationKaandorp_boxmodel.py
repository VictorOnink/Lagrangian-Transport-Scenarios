import settings
import utils
import numpy as np
from scipy.special import gamma
from copy import deepcopy


class FragmentationKaandorp_box_model:
    def __init__(self, sim_length, start_size=settings.INIT_SIZE, size_classes=settings.SIZE_CLASS_NUMBER,
                 DN=settings.DN, p_frag=settings.P_FRAG, lambda_f=settings.LAMBDA_FRAG, lambda_b=settings.SHORE_TIME,
                 rho=settings.INIT_DENSITY, lambda_f_ocean=settings.LAMBDA_OCEAN_FRAG, ocean_frag=False,
                 steady_state=False):
        # Data directories
        self.input_direc = settings.DATA_INPUT_DIR_SERVERS[settings.SERVER] + 'box_model/'
        self.output_direc = settings.DATA_OUTPUT_DIR_SERVERS[settings.SERVER] + 'box_model/'
        # Fragmentation parameters
        self.start_size = start_size
        self.size_classes = size_classes
        self.size_array = utils.size_range(size_class_number=self.size_classes)
        self.class_array = np.arange(0, self.size_classes)
        self.DN = DN
        self.rho = rho
        self.sim_length = sim_length
        self.steady_state = steady_state
        # Stokes drift dependence parameters
        self.p_shape = 'sp'
        self.quantile = 0.5
        self.lambda_b = lambda_b
        self.tau_cb = np.ones(self.class_array.shape) * self.lambda_b  # beaching timescale (days)
        self.p_sink = 5.6e-4  # Sink probability for size dependent resuspension
        self.lambda_f = lambda_f
        self.i_f = lambda_f_to_i_f(timescale=self.lambda_f)
        self.lambda_f_ocean = lambda_f_ocean
        self.i_f_ocean = lambda_f_to_i_f(timescale=self.lambda_f_ocean)
        self.ocean_frag = ocean_frag
        self.p_frag = p_frag  # fragmentation probability
        # The transition matrices for mass and particle number, matrices for the fragmentation, and a dictionary to save
        # the model parameters
        self.T_mat_N = None
        self.T_mat_m = None
        self.T_NB_m = np.zeros((self.size_classes, self.size_classes))
        self.T_NB_N = np.zeros((self.size_classes, self.size_classes))
        self.T_NB_mO = np.zeros((self.size_classes, self.size_classes))
        self.T_NB_NO = np.zeros((self.size_classes, self.size_classes))
        self.dict_save = {}
        # Indices for the transition matrices
        self.index_o = np.arange(self.size_classes)
        self.index_c = np.arange(self.size_classes, 2 * self.size_classes)
        self.index_b = np.arange(2 * self.size_classes, 3 * self.size_classes)

    def fragmentation_fractions(self):
        '''
        The fragmentation model by Charalambous (2015)
        k_arr: the size classes
        i_f: the fragmentation 'index' (denoted by f in the paper)
        p: fragmentation probability
        d_N: spatial dimensionality
        '''
        # Beach fragmentation variables
        gamma_ratio = (gamma(self.class_array + self.i_f) / (gamma(self.class_array + 1) * gamma(self.i_f)))
        pmf_m = gamma_ratio * (self.p_frag ** self.class_array) * (1 - self.p_frag) ** self.i_f
        pmf_N = 2 ** (self.DN * self.class_array) * pmf_m
        # Ocean fragmentation variables
        gamma_ratio_ocean = (gamma(self.class_array + self.i_f_ocean) / (gamma(self.class_array + 1) * gamma(self.i_f_ocean)))
        pmf_mO = gamma_ratio_ocean * (self.p_frag ** self.class_array) * (1 - self.p_frag) ** self.i_f_ocean
        pmf_NO = 2 ** (self.DN * self.class_array) * pmf_mO
        return pmf_m, pmf_N, pmf_mO, pmf_NO

    def get_fragmentation_probabilities(self):
        # Calculating the weekly fragmentation fractions
        m_NB_dt, N_NB_dt, m_NB_dt_O, N_NB_dt_O = self.fragmentation_fractions()
        for index_class in range(self.size_classes):
            self.T_NB_m += m_NB_dt[index_class] * np.diag(np.ones(self.size_classes - index_class), -index_class)
            self.T_NB_N += N_NB_dt[index_class] * np.diag(np.ones(self.size_classes - index_class), -index_class)
            self.T_NB_mO += m_NB_dt_O[index_class] * np.diag(np.ones(self.size_classes - index_class), -index_class)
            self.T_NB_NO += N_NB_dt_O[index_class] * np.diag(np.ones(self.size_classes - index_class), -index_class)

    def get_ocean_probabilities(self):
        probability_data = utils.load_obj(self.input_direc + 'Stokes_analysis_-06373046_1590.pkl')
        stokes_data = utils.load_obj(self.input_direc + 'Stokes_influence_factor.pkl')

        # Determining the stokes influence factor for the particle sizes in self.size_array
        SIF = np.interp(self.size_array, np.flip(stokes_data['l']), np.flip(stokes_data[self.p_shape][self.quantile]))

        # Loading the transition probabilities as functions of the SIF
        SIF_sim = np.linspace(0, 1, 6, endpoint=True)
        P_oo_sim = np.array([probability_data['p_oo_dist'][key] for key in probability_data['p_oo_dist'].keys()])[:-1]
        P_oc_sim = np.array([probability_data['p_oc_dist'][key] for key in probability_data['p_oc_dist'].keys()])[:-1]
        P_co_sim = np.array([probability_data['p_co_dist'][key] for key in probability_data['p_co_dist'].keys()])[:-1]
        P_cc_sim = np.array([probability_data['p_cc_dist'][key] for key in probability_data['p_cc_dist'].keys()])[:-1]

        # Interpolating the transition probabilities for the particles sizes in self.size_array
        P_oo = np.interp(SIF, SIF_sim, P_oo_sim)
        P_oc = np.interp(SIF, SIF_sim, P_oc_sim)
        P_cc = np.interp(SIF, SIF_sim, P_cc_sim)
        P_co = np.interp(SIF, SIF_sim, P_co_sim)

        return P_oo, P_oc, P_co, P_cc

    def get_resuspension_probability(self):
        '''
        resuspension time scale from Hinata et al. (2017)
        '''
        data_S = utils.load_obj(self.input_direc + 'Stokes_influence_factor.pkl')
        # get stokes influence factor for the l_arr
        wb = np.interp(self.size_array, np.flip(data_S['l']), np.flip(data_S['wb_' + self.p_shape]))
        tau_bc = 2.6e2 * wb + 7.1
        return tau_bc

    def create_transition_matrix(self):
        P_oo_sim, P_oc_sim, P_co_sim, P_cc_sim = self.get_ocean_probabilities()

        # Setting the ocean parameters
        P_os = self.p_sink
        P_oo = (P_oo_sim / (P_oo_sim + P_oc_sim)) * (1 - P_os)
        P_oc = (P_oc_sim / (P_oo_sim + P_oc_sim)) * (1 - P_os)
        assert (np.abs(1 - (P_oo + P_oc + P_os)) < 1e-10).all()

        # Setting the coastal parameters
        P_cs = self.p_sink
        P_cb = tau_to_lambda(self.tau_cb / 7)  # beaching probability for a week
        P_norm_co = P_co_sim / (P_co_sim + P_cc_sim)  # normalize to one
        P_norm_cc = P_cc_sim / (P_co_sim + P_cc_sim)

        P_norm2_co = P_norm_co * (1 - P_cb)  # equal fraction needs to be substracted which ends up on beach
        P_norm2_cc = P_norm_cc * (1 - P_cb)

        P_co = P_norm2_co * (1 - P_cs)  # sink fractions substracted as well
        P_cc = P_norm2_cc * (1 - P_cs)
        P_cb = P_cb * (1 - P_cs)
        assert (np.abs(1 - (P_cc + P_co + P_cs + P_cb)) < 1e-10).all()

        # Setting the beach parameters
        P_bs = self.p_sink
        P_sim_bc = tau_to_lambda(self.get_resuspension_probability() / 7)
        P_sim_bb = 1 - P_sim_bc
        P_bb = (P_sim_bb / (P_sim_bb + P_sim_bc)) * (1 - P_bs)
        P_bc = (P_sim_bc / (P_sim_bb + P_sim_bc)) * (1 - P_bs)
        assert (np.abs(1 - (P_bb + P_bc + P_bs)) < 1e-10).all()

        # Create a dictionary with matrix probabilities
        P_mat = {}
        for variable, name in zip([P_oo, P_oc, P_co, P_cc, P_cb, P_bc, P_bb],
                                  ['oo', 'oc', 'co', 'cc', 'cb', 'bc', 'bb']):
            P_mat[name] = np.diag(variable)
        self.get_fragmentation_probabilities()
        # Beach fragmentation
        for variable, name in zip([P_bc, P_bb], ['bc', 'bb']):
            P_mat[name + '_N'] = (np.ones(self.T_NB_N.shape) * variable).T * self.T_NB_N
            P_mat[name + '_m'] = (np.ones(self.T_NB_N.shape) * variable).T * self.T_NB_m
        # Ocean
        if self.ocean_frag:
            for variable, name in zip([P_oo, P_oc, P_co, P_cc, P_cb], ['oo', 'oc', 'co', 'cc', 'cb']):
                P_mat[name + '_N'] = (np.ones(self.T_NB_N.shape) * variable).T * self.T_NB_NO
                P_mat[name + '_m'] = (np.ones(self.T_NB_N.shape) * variable).T * self.T_NB_mO

        # Finally, finishing the transition matrix
        self.T_mat_N = np.zeros([3 * self.size_classes, 3 * self.size_classes])
        self.T_mat_m = np.zeros([3 * self.size_classes, 3 * self.size_classes])
        if self.ocean_frag:
            self.T_mat_m[self.index_o[0]:self.index_o[-1] + 1, self.index_o[0]:self.index_o[-1] + 1] = P_mat['oo_m']
            self.T_mat_m[self.index_o[0]:self.index_o[-1] + 1, self.index_c[0]:self.index_c[-1] + 1] = P_mat['co_m']
            self.T_mat_m[self.index_c[0]:self.index_c[-1] + 1, self.index_o[0]:self.index_o[-1] + 1] = P_mat['oc_m']
            self.T_mat_m[self.index_c[0]:self.index_c[-1] + 1, self.index_c[0]:self.index_c[-1] + 1] = P_mat['cc_m']
            self.T_mat_m[self.index_b[0]:self.index_b[-1] + 1, self.index_c[0]:self.index_c[-1] + 1] = P_mat['cb_m']
        else:
            self.T_mat_m[self.index_o[0]:self.index_o[-1] + 1, self.index_o[0]:self.index_o[-1] + 1] = P_mat['oo']
            self.T_mat_m[self.index_o[0]:self.index_o[-1] + 1, self.index_c[0]:self.index_c[-1] + 1] = P_mat['co']
            self.T_mat_m[self.index_c[0]:self.index_c[-1] + 1, self.index_o[0]:self.index_o[-1] + 1] = P_mat['oc']
            self.T_mat_m[self.index_c[0]:self.index_c[-1] + 1, self.index_c[0]:self.index_c[-1] + 1] = P_mat['cc']
            self.T_mat_m[self.index_b[0]:self.index_b[-1] + 1, self.index_c[0]:self.index_c[-1] + 1] = P_mat['cb']
        self.T_mat_m[self.index_c[0]:self.index_c[-1] + 1, self.index_b[0]:self.index_b[-1] + 1] = P_mat['bc_m']
        self.T_mat_m[self.index_b[0]:self.index_b[-1] + 1, self.index_b[0]:self.index_b[-1] + 1] = P_mat['bb_m']

        if self.ocean_frag:
            self.T_mat_N[self.index_o[0]:self.index_o[-1] + 1, self.index_o[0]:self.index_o[-1] + 1] = P_mat['oo_N']
            self.T_mat_N[self.index_o[0]:self.index_o[-1] + 1, self.index_c[0]:self.index_c[-1] + 1] = P_mat['co_N']
            self.T_mat_N[self.index_c[0]:self.index_c[-1] + 1, self.index_o[0]:self.index_o[-1] + 1] = P_mat['oc_N']
            self.T_mat_N[self.index_c[0]:self.index_c[-1] + 1, self.index_c[0]:self.index_c[-1] + 1] = P_mat['cc_N']
            self.T_mat_N[self.index_b[0]:self.index_b[-1] + 1, self.index_c[0]:self.index_c[-1] + 1] = P_mat['cb_N']
        else:
            self.T_mat_N[self.index_o[0]:self.index_o[-1] + 1, self.index_o[0]:self.index_o[-1] + 1] = P_mat['oo']
            self.T_mat_N[self.index_o[0]:self.index_o[-1] + 1, self.index_c[0]:self.index_c[-1] + 1] = P_mat['co']
            self.T_mat_N[self.index_c[0]:self.index_c[-1] + 1, self.index_o[0]:self.index_o[-1] + 1] = P_mat['oc']
            self.T_mat_N[self.index_c[0]:self.index_c[-1] + 1, self.index_c[0]:self.index_c[-1] + 1] = P_mat['cc']
            self.T_mat_N[self.index_b[0]:self.index_b[-1] + 1, self.index_c[0]:self.index_c[-1] + 1] = P_mat['cb']
        self.T_mat_N[self.index_c[0]:self.index_c[-1] + 1, self.index_b[0]:self.index_b[-1] + 1] = P_mat['bc_N']
        self.T_mat_N[self.index_b[0]:self.index_b[-1] + 1, self.index_b[0]:self.index_b[-1] + 1] = P_mat['bb_N']

        # Saving the model parameters
        for name, variable in zip(['oo', 'oc', 'co', 'cc', 'cb', 'bb', 'bc', 'p_sink', 'lambda_frag', 'lambda_b'],
                                  [P_oo, P_oc, P_co, P_cc, P_cb, P_bb, P_bc, self.p_sink, self.lambda_f, self.lambda_b]):
            self.dict_save[name] = variable
        if self.ocean_frag:
            self.dict_save['lambda_frag_ocean'] = self.lambda_f_ocean

        return self.T_mat_m, self.T_mat_N, self.dict_save

    def weekly_input(self):
        mass_0, number_0 = np.zeros(3 * self.size_classes), np.zeros(3 * self.size_classes)
        input_week = 453 * 12 / 365 * 7
        mass_0[self.index_b[0]], number_0[self.index_b[0]] = input_week, input_week
        return mass_0, number_0

    def calculate_transport(self):
        self.T_mat_m, self.T_mat_N, self.dict_save = self.create_transition_matrix()

        # Create an input array, with input in the biggest size class
        mass_0, number_0 = self.weekly_input()
        mass, number = deepcopy(mass_0), deepcopy(number_0)

        # Creating an output dictionary
        output_dict_mass = {'size': self.size_array}
        output_dict_number = {'size': self.size_array}

        for time in range(0, int(52 * self.sim_length)):
            mass = np.dot(self.T_mat_m, mass + mass_0)
            number = np.dot(self.T_mat_N, number + number_0)
            output_dict_mass[time], output_dict_number[time] = {}, {}
            for reservoir, indices in zip(['ocean', 'coastal', 'beach'], [self.index_o, self.index_c, self.index_b]):
                output_dict_mass[time][reservoir] = mass[indices]
                output_dict_number[time][reservoir] = number[indices]
            output_dict_mass[time]['total'] = mass[self.index_o] + mass[self.index_b] + mass[self.index_c]
            output_dict_number[time]['total'] = number[self.index_o] + number[self.index_b] + number[self.index_c]
        output_dict_mass['final_index'], output_dict_number['final_index'] = time, time

        self.dict_save['mass'] = output_dict_mass
        self.dict_save['number'] = output_dict_number

        return self.dict_save

    def load_box_model(self, rerun=False):
        file_name = self.get_file_name()
        if utils.check_file_exist(file_name, without_pkl=True) and not rerun:
            return utils.load_obj(file_name)
        else:
            if self.steady_state:
                output_dict = self.steady_state_distribution()
            else:
                output_dict = self.calculate_transport()
            utils.save_obj(filename=file_name, item=output_dict)
            return output_dict

    def get_file_name(self):
        str_format = self.lambda_b, self.p_frag, self.lambda_f, self.DN, self.size_classes, self.rho
        file_name = self.output_direc + 'box_model_st={}_pfrag={}_lambdafrag={}_dn={}_sizeclasses={}_rho={}'.format(*str_format)
        if self.ocean_frag:
            file_name += 'lambda_f_ocean={}'.format(self.lambda_f_ocean)
        if self.steady_state:
            file_name += '_steady_state'
        return file_name

    def steady_state_distribution(self):
        # Obtain the transition matrix
        self.T_mat_m, self.T_mat_N, self.dict_save = self.create_transition_matrix()
        # Setting the weekly input
        mass_input, number_input = self.weekly_input()
        # Now following the approach at, first for the mass and then the numbers.
        # https://github.com/OceanParcels/ContinuousCascadingFragmentation/blob/main/box_model_Mediterranean.py#L317-L341
        m1 = np.linalg.inv(np.eye(mass_input.size) - self.T_mat_m)
        m2 = np.dot(self.T_mat_m, mass_input)
        mass_steady_state = np.dot(m1, m2)

        n1 = np.linalg.inv(np.eye(number_input.size) - self.T_mat_N)
        n2 = np.dot(self.T_mat_m, number_input)
        number_steady_state = np.dot(n1, n2)

        # Saving the steady state distribution for each reservoir
        self.dict_save['mass'], self.dict_save['number'] = {}, {}
        for reservoir, indices in zip(['ocean', 'coastal', 'beach'], [self.index_o, self.index_c, self.index_b]):
            self.dict_save['mass'][reservoir] = mass_steady_state[indices]
            self.dict_save['number'][reservoir] = number_steady_state[indices]
        return self.dict_save


def tau_to_lambda(tau):
    return 1 - np.exp(-1 / tau)


def lambda_f_to_i_f(timescale):
    # Convert fragmentation timescale (days) to f / week
    return 1 / (timescale / 7)
