import os
from datetime import timedelta
import socket
if socket.gethostname() in ['kuphaven', 'climatestor01']:
    from dotenv import load_dotenv
    load_dotenv()


def load_env_variable(variable, default, variable_type=int):
    """
    Loading environmental variables
    :param variable: name of the variable
    :param default: default value if the environmental variable is not defined
    :param variable_type: form that the variable should take
    :return:
    """
    if variable in os.environ:
        return {int: int(os.environ[variable]), float: float(os.environ[variable]), str: str(os.environ[variable]),
                bool: bool(int(os.environ[variable]))}[variable_type]
    else:
        return default


SUBMISSION = load_env_variable("SUBMISSION", default=None, variable_type=str)
os.system('echo "run="' + SUBMISSION)

########################################################################################################################
#                                                                                                                      #
#                                    Parameters for data retrieval/storage                                             #
#                                                                                                                      #
########################################################################################################################

# DIRECTORIES FOR DATA, INPUTS & OUTPUTS
SERVER: str = load_env_variable("SERVER", default="UBELIX", variable_type=str)
DATA_DIR_SERVERS: dict = {'KUPHAVEN': "/storage/climatestor/Bern3dLPX/onink/alphadata04/lagrangian_sim/",
                          'UBELIX': "/storage/homefs/vo18e689/Data/"}
DATA_INPUT_DIR_SERVERS: dict = {'KUPHAVEN': "/storage/climatestor/Bern3dLPX/onink/alphadata04/lagrangian_sim/Input/",
                                'UBELIX': "/storage/homefs/vo18e689/Data/Input/"}
DATA_OUTPUT_DIR_SERVERS: dict = {'KUPHAVEN': "/storage/climatestor/Bern3dLPX/onink/alphadata04/lagrangian_sim/Output/",
                                 'UBELIX': "/storage/homefs/vo18e689/Data/Output/"}
DATA_WORKSPACE_OUTPUT_DIR: dict = {'KUPHAVEN': None,
                                   'UBELIX': "/storage/workspaces/climate_esm_bgc/bern3d_lpx/onink/"}
FIGURE_OUTPUT_SERVER: dict = {'KUPHAVEN': "/storage/climatestor/Bern3dLPX/onink/alphadata04/lagrangian_sim/Output/Figures/",
                              'UBELIX': "/storage/homefs/vo18e689/Data/Output/Figures/"}
INPUT_DIREC_DICT = {0: DATA_INPUT_DIR_SERVERS[SERVER] + 'Jambeck_Inputs/',
                    1: DATA_INPUT_DIR_SERVERS[SERVER] + 'Lebreton_Inputs/',
                    2: DATA_INPUT_DIR_SERVERS[SERVER] + 'LebretonDivision_Inputs/',
                    3: DATA_INPUT_DIR_SERVERS[SERVER] + 'Lebreton_Kaandorp_init_Inputs/',
                    4: DATA_INPUT_DIR_SERVERS[SERVER] + 'Point_Release/',
                    5: DATA_INPUT_DIR_SERVERS[SERVER] + 'Uniform/'}
SCRATCH_DIR = {'KUPHAVEN': None,
               'UBELIX': "/storage/scratch/users/vo18e689/"}[SERVER]

# STARTING YEAR, MONTH AND DAY OF THE SIMULATION
STARTYEAR: int = load_env_variable("STARTYEAR", default=2010)
STARTMONTH: int = load_env_variable("STARTMONTH", default=1)
STARTDAY: int = load_env_variable("STARTDAY", default=1)

# LENGTH IN YEARS OF THE SIMULATION
SIM_LENGTH: int = load_env_variable("SIMLEN", default=1)

########################################################################################################################
#                                                                                                                      #
#                                         Run/restart specific parameters                                              #
#                                                                                                                      #
########################################################################################################################

# WHICH OF THE SUBRUNS (DIVISION OF PARTICLE INPUTS)
RUN: int = load_env_variable("RUN", default=0)

# RESTART NUMBER OF THE RUN
# 0 = NEW RUN, N>0 INDICATES NTH YEAR FROM SIMULATION START
RESTART: int = load_env_variable("RESTARTNUM", default=0)

# LAGRANGIAN OR POSTPROCESSING RUN
# 0 = RUN LAGRANGIAN SIMULATION, 1 = RUN POST PROCESSING
POST_PROCESS: bool = load_env_variable("POST_PROCESS", default=False, variable_type=bool)

# ENSEMBLE NUMBER
ENSEMBLE: int = load_env_variable("ENSEMBLE", default=1)

########################################################################################################################
#                                                                                                                      #
#                                       Advection scenario specific parameters                                         #
#                                                                                                                      #
########################################################################################################################

# UV ADVECTION DATA
ADVECTION_DICT: dict = {0: 'HYCOM_GLOBAL', 1: 'HYCOM_CARIBBEAN', 2: 'CMEMS_MEDITERRANEAN'}
ADVECTION_DATA: str = ADVECTION_DICT[load_env_variable("ADVECTION_DATA", default=2)]
# STOKES DRIFT: 0 -> STOKES, 1 -> NO STOKES
STOKES: int = load_env_variable("STOKES", default=0)

########################################################################################################################
#                                                                                                                      #
#                      Scenario specific parameters, not considering input and advection scenarios                     #
#                                                                                                                      #
########################################################################################################################
# MODEL SCENARIO SETTINGS
SCENARIO_NAME: str = load_env_variable("SCENARIO", default=None, variable_type=str)
assert SCENARIO_NAME is not None, 'Please pick a scenario!'


# FOR 'CoastalProximity': TIME IN BEACHING ZONE PRIOR TO BEACHING (DAYS)
VICINITY: int = load_env_variable("VICINITY", default=1)

# BEACHING TIMESCALE (DAYS)
SHORE_TIME: int = load_env_variable("SHORETIME", default=26)
# RESUSPENSION TIMESCALE (DAYS)
RESUS_TIME: int = load_env_variable("RESUSTIME", default=69)
# FIXED_RESUS OR SIZE DEPENDENT RESUSPENSION TIMESCALE
FIXED_RESUS: bool = load_env_variable("FIXED_RESUS", default=False, variable_type=bool)

# FOR 'ShoreDependentResuspension': 0 -> SAMARAS ET AL (2015) RESUSPENSION DEPENDENCE,
#                                   1 -> 1:4 RESUSPENSION DEPENDENCE
SHORE_DEP: int = load_env_variable("SHOREDEPEN", default=0)

# FOR 'TurrellResuspension': MINIMUM OFF-SHORE WIND SPEED FOR RESUSPENSION
WMIN: float = load_env_variable("WMIN", default=0.0, variable_type=float)

# FOR 'FragmentationKaandorpPartial': WHETHER WE WANT TO INCLUDE OCEAN FRAGMENTATION OR NOT
OCEAN_FRAG: bool = load_env_variable("OCEAN_FRAG", default=False, variable_type=bool)


########################################################################################################################
#                                                                                                                      #
#                                       Input scenario specific parameters                                             #
#                                                                                                                      #
########################################################################################################################
# THE INPUT SCENARIO
INPUT_NAMES = {0: 'Jambeck', 1: 'Lebreton', 2: 'LebretonDivision', 3: 'LebretonKaandorpInit', 4: 'Point_Release',
               5: 'Uniform'}
INPUT = INPUT_NAMES[load_env_variable("INPUT", default=0)]
# DIRECTORY CONTAINING INITIAL INPUTS FOR RESTART == 0
INPUT_DIREC = INPUT_DIREC_DICT[load_env_variable("INPUT", default=0)]

if INPUT == 'Jambeck':
    # THE NUMBER OF PARTICLES PER RELEASE STEP PER RUN
    INPUT_DIV = 5000
    # NUMBER OF RUNS
    RUN_RANGE: int = 9
    # MAXIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MAX = 10.0
    # MINIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MIN = 0.08
elif INPUT == 'Lebreton':
    # THE NUMBER OF PARTICLES PER RELEASE STEP PER RUN
    INPUT_DIV = 200000
    # NUMBER OF RUNS
    RUN_RANGE: int = 1
    # MAXIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MAX = 0.0125
    # MINIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MIN = 0.00001  # Minimum plastic mass input for a cell in order to be considered for the input
elif INPUT == 'LebretonDivision':
    # THE NUMBER OF PARTICLES PER RELEASE STEP PER RUN
    INPUT_DIV = 50
    # NUMBER OF RUNS
    RUN_RANGE: int = 10
    # MAXIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MAX = 5
    # MINIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MIN = 0.01  # Minimum plastic mass input for a cell in order to be considered for the input
elif INPUT == 'LebretonKaandorpInit':
    # THE NUMBER OF PARTICLES PER RELEASE STEP PER RUN
    INPUT_DIV = 50
    # NUMBER OF RUNS
    RUN_RANGE: int = 10
    # MAXIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MAX = 5
    # MINIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MIN = 0.01  # Minimum plastic mass input for a cell in order to be considered for the input
elif INPUT == 'Point_Release':
    # THE RELEASE SITE OF THE POINT RELEASE
    RELEASE_SITE = load_env_variable("RELEASE_SITE", default=0)
    SITE_LONLAT = {0: (11.0, 42.2)}[RELEASE_SITE]
    # NUMBER OF RELEASE PARTICLES
    PARTICLE_NUMBER = 275
    # THE NUMBER OF PARTICLES PER RELEASE STEP PER RUN
    INPUT_DIV = 1000000
    # NUMBER OF RUNS
    RUN_RANGE: int = 1
    # STARTING LONGITUDE FOR POINT RELEASE (DEGREES EAST)
    INPUT_LON = SITE_LONLAT[0]
    # STARTING LATITUDE FOR POINT RELEASE (DEGREES NORTH)
    INPUT_LAT = SITE_LONLAT[1]
    # MAXIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MAX = 1.0
    # MINIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MIN = 0.0
elif INPUT == 'Uniform':
    # THE NUMBER OF PARTICLES PER RELEASE STEP PER RUN
    INPUT_DIV = 5000
    # Number of runs
    RUN_RANGE: int = 1
    # GRID RESOlUTION ON WHICH PARTICLES ARE RELEASED
    RELEASE_GRID = 0.1
    # MAXIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MAX = 1.0
    # MINIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MIN = 0.0

########################################################################################################################
#                                                                                                                      #
#                                       Analysis specific parameters                                                   #
#                                                                                                                      #
########################################################################################################################
# IF PARALLEL_STEP == 1, THEN WE ARE CALCULATING THE ANALYSIS FOR EACH RUN/RESTART FILE
# IF PARALLEL_STEP == 2, THEN WE COMBINE ALL THE INDIVIDUAL RUN/RESTART ANALYSIS OUTPUT FILES INTO ONE FOR THE WHOLE
#                        MODEL SETUP
PARALLEL_STEP: int = load_env_variable("PARALLEL_STEP", default=0)
# PLASTIC CONCENTRATION
CONCENTRATION = load_env_variable("CONCENTRATION", default=False, variable_type=bool)
# PLASTIC CONCENTRATIONS (MONTHLY CONCENTRATIONS)
MONTHLY_CONCENTRATION = load_env_variable("MONTHLY_CONCENTRATION", default=False, variable_type=bool)
# VERTICAL PLASTIC CONCENTRATION
VERTICAL_CONCENTRATION = load_env_variable("VERTICAL_CONCENTRATION", default=False, variable_type=bool)
# CALCULATE THE VERTICAL PROFILES ON A SPATIAL GRID
SPATIAL_VERTICAL_PROFILES = load_env_variable("SPATIAL_VERTICAL_PROFILES", default=False, variable_type=bool)
# LON LAT AVERAGED CONCENTRATIONS
LONLAT_CONCENTRATION = load_env_variable("LONLAT_CONCENTRATION", default=False, variable_type=bool)
# TIMESERIES OF BEACHED/COASTAL FRACTIONS
TIMESERIES = load_env_variable("TIMESERIES", default=False, variable_type=bool)
# MAX DISTANCE FROM SHORE
MAX_DISTANCE = load_env_variable("MAX_DISTANCE", default=False, variable_type=bool)
# CREATING TIME SLICES OF THE SIMULATION
TIMESLICING = load_env_variable("TIMESLICING", default=False, variable_type=bool)
# CALCULATING BASIC STATISTICS FOR EACH PARTICLE TRAJECTORY
STATISTICS = load_env_variable("STATISTICS", default=False, variable_type=bool)
# CAlCULATING SEPARATION DISTANCES BETWEEN PARTICLES
SEPARATION = load_env_variable("SEPARATION", default=False, variable_type=bool)
# CALCULATING THE PARTICLE SIZE SPECTRUM/DISTRIBUTION
SIZE_SPECTRUM = load_env_variable("SIZE_SPECTRUM", default=False, variable_type=bool)
# CALCULATING BAYESIAN MATRIX KEEPING TRACK OF
BAYESIAN = load_env_variable("BAYESIAN", default=False, variable_type=bool)
# FRACTION OF PLASTIC ENTERING THE OCEAN IN 2010 THAT IS CONSIDERED BUOYANT (BASED ON GEYERS ET AL., 2017)
BUOYANT = 0.54

########################################################################################################################
#                                                                                                                      #
#                                       Model and physical parameters                                                  #
#                                                                                                                      #
########################################################################################################################
# FORWARD OR BACKWARDS IN TIME
BACKWARD = load_env_variable("BACKWARD", default=False, variable_type=bool)
BACKWARD_MULT = {True: -1, False: 1}[BACKWARD]
# MODEL INTEGRATION TIMESTEP
if SCENARIO_NAME == 'BlueCloud':
    DT_MINUTES = 5
else:
    DT_MINUTES = 0.5
TIME_STEP = timedelta(minutes=DT_MINUTES * BACKWARD_MULT)
# MODEL OUTPUT TIMESTEP
OUTPUT_TIME_STEP = timedelta(hours=12)
# PARTICLE RELEASE TIMESTEP
REPEAT_DT_R0 = timedelta(days=31)
REPEAT_DT_ELSE = None
# RNG SEED VALUE
SEED = 'Fixed'
# HORIZONTAL DIFFUSIVITY M^2 / S, FOLLOWING LACERDA ET AL. (2019) AND LIUBERTSEVA ET AL. (2018)
K_HOR = 10
# VERTICAL (DIAPYCNAL) DIFFUSION (M^2/S) BELOW THE MLD (WATERHOUSE ET AL, 2014)
K_Z_BULK = 3e-5
# WIDTH OF THE BEACHING ZONE (KM) WITHIN WHICH BEACHING CAN OCCUR
COAST_D = {'HYCOM_GLOBAL': 10, 'HYCOM_CARIBBEAN': 4, 'CMEMS_MEDITERRANEAN': 6}[ADVECTION_DATA]
# SIZE SCALING FACTOR FOR INIT_SIZE
SIZE_FACTOR = 1E-6

# INITIAL PARTICLE SIZE (M)
INIT_SIZE = load_env_variable("PARTICLE_SIZE", default=5000) * SIZE_FACTOR
# INITIAL DENSITY (KG/M^3): 920 = POLYPROPYLENE, 980 = HIGH DENSITY POLYETHYLENE (BRIGNAC ET AL. 2017)
INIT_DENSITY = load_env_variable("INIT_DENSITY", default=920)
# FRAGMENTATION PROBABILITY
P_FRAG = load_env_variable("P", default=0.4, variable_type=float)
# NUMBER OF SPATIAL DIMENSIONS
DN = load_env_variable("DN", default=2.5, variable_type=float)
# NUMBER OF SIZE CLASSES
SIZE_CLASS_NUMBER = load_env_variable("SIZE_CLASS_NUMBER", default=6)
# FRAGMENTATION TIMESCALE (DAYS)
LAMBDA_FRAG = load_env_variable("LAMBDA_FRAG", default=388)
# OPEN OCEAN FRAGMENTATION TIMESCALE (DAYS)
LAMBDA_OCEAN_FRAG = load_env_variable("LAMBDA_OCEAN_FRAG", default=388)
# CRITICAL SHEAR STRESS FOR RESUSPENSION OF PARTICLES FROM THE SEA BED, HORIZONTAL DIFFUSION AT THE SEA BED
SEABED_CRIT = load_env_variable("SEABED_CRIT", default=0, variable_type=float)

# ACCELERATION DUE TO GRAVITY (M/S^2)
G = 9.81
# THERMAL EXPANSION COEFFICIENT (1/K) AND REFERENCE TEMPERATURE (K)
A_T, T_R = 2e-4, 25
# HALINE CONTRACTION COEFFICIENT (1/PSU) AND REFERENCE SALINITY (PSU)
B_S, S_R = 8e-4, 35
# DENSITY OF AIR (KG/M^3)
RHO_A = 1.22
# VON KARMAN COEFFICIENT
VK = 0.4
# WAVE AGE FOR FULLY DEVELOPED SEA (KUKULKA ET AL., 2012), BASED ON EITHER WIND SPEED OR AIR FRICTIONAL VELOCITY
BETA, BETA_STAR = 1.21, 35
# STABILITY FUNCTION IN MONIN-OBUKOV BOUNDARY LAYER THEORY (BOUFADEL ET AL. 2020)
PHI = 0.9
# HORIZONTAL DIFFUSION AT THE SEA BED
SEABED_KH = 0.2
# MICROPLASTIC REMOVAL RATE (KAANDORP ET AL., 2020), CONVERT FROM 5.3 X 10**-3 WEEK**-1 TO PER TIMESTEP
P_SINK = 5.1e-3 / (604800 / TIME_STEP.total_seconds())

########################################################################################################################
#                                                                                                                      #
#                                                     LOG                                                              #
#                                                                                                                      #
########################################################################################################################
os.system('echo "This is the {} for the {} scenario"'.format(SUBMISSION, SCENARIO_NAME))
os.system('echo "The input scenario is {}"'.format(INPUT))
if SUBMISSION == 'simulation':
    os.system('echo "The starting year is {}-{}, and this is run {}, restart {}"'.format(STARTYEAR, STARTMONTH, RUN,
                                                                                         RESTART))
    os.system('echo "We are using {} for the plastic advection"'.format(ADVECTION_DATA))
if SUBMISSION == 'analysis':
    os.system('echo "We are running analysis parallelization step {}"'.format(PARALLEL_STEP))

if SCENARIO_NAME == 'CoastalProximity':
    os.system('echo "The coastal vicinity time cutoff is {} days "'.format(VICINITY))

elif SCENARIO_NAME == "Stochastic":
    os.system('echo "The beaching timescale is {} days "'.format(SHORE_TIME))
    os.system('echo "The resuspension timescale is {} days "'.format(RESUS_TIME))

elif SCENARIO_NAME == 'ShoreDependentResuspension':
    os.system('echo "The beaching timescale is {} days "'.format(SHORE_TIME))
    os.system('echo "The resuspension timescale is {} days "'.format(RESUS_TIME))
    os.system('echo "The shore dependent resuspension scenario is {}"'.format(SHORE_DEP))

elif SCENARIO_NAME == 'TurrellResuspension':
    os.system('echo "The beaching timescale is {} days "'.format(SHORE_TIME))
    os.system('echo "The minimum offshore wind for resuspension is {} m / s"'.format(WMIN/10.))

elif SCENARIO_NAME == 'SizeTransport':
    os.system('echo "The particle size is {:.3f} mm"'.format(INIT_SIZE * 1e3))
    os.system('echo "The particle density is {} kg/m3 "'.format(INIT_DENSITY))
    os.system('echo "The beaching timescale is {} days "'.format(SHORE_TIME))
    if FIXED_RESUS:
        os.system('echo "The resuspension timescale is fixed at {} days "'.format(RESUS_TIME))
    else:
        os.system('echo "The resuspension timescale is size dependent"')

elif SCENARIO_NAME in ['FragmentationKaandorp', 'FragmentationKaandorpPartial']:
    os.system('echo "The beaching timescale is {} days "'.format(SHORE_TIME))
    os.system('echo "The fragmentation timescale is {} days "'.format(LAMBDA_FRAG))
    os.system('echo "The initial particle size is {} m and a rho of {} kg/m^3"'.format(INIT_SIZE, INIT_DENSITY))
    os.system('echo "This simulation has {} size classes"'.format(SIZE_CLASS_NUMBER))
