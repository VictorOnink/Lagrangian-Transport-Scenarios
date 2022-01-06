import os
from datetime import timedelta
import socket
if socket.gethostname() in ['kuphaven', 'climatestor01']:
    from dotenv import load_dotenv
    load_dotenv()

SUBMISSION = str(os.environ["SUBMISSION"])
os.system('echo "run="' + SUBMISSION)


########################################################################################################################
#                                                                                                                      #
#                                    Parameters for data retrieval/storage                                             #
#                                                                                                                      #
########################################################################################################################

# DIRECTORIES FOR DATA, INPUTS & OUTPUTS
SERVER: int = int(os.environ["SERVER"])
DATA_DIR_SERVERS: dict = {0: "/storage/climatestor/Bern3dLPX/onink/alphadata04/lagrangian_sim/",
                          1: "/storage/homefs/vo18e689/Data/"}
DATA_INPUT_DIR_SERVERS: dict = {0: "/storage/climatestor/Bern3dLPX/onink/alphadata04/lagrangian_sim/Input/",
                                1: "/storage/homefs/vo18e689/Data/Input/"}
DATA_OUTPUT_DIR_SERVERS: dict = {0: "/storage/climatestor/Bern3dLPX/onink/alphadata04/lagrangian_sim/Output/",
                                 1: "/storage/homefs/vo18e689/Data/Output/"}
DATA_WORKSPACE_OUTPUT_DIR: dict = {0: None,
                                   1: "/storage/workspaces/climate_esm_bgc/bern3d_lpx/onink/"}
FIGURE_OUTPUT_SERVER: dict = {0: "/storage/climatestor/Bern3dLPX/onink/alphadata04/lagrangian_sim/Output/Figures/",
                              1: "/storage/homefs/vo18e689/Data/Output/Figures/"}
INPUT_DIREC_DICT = {0: DATA_INPUT_DIR_SERVERS[SERVER] + 'Jambeck_Inputs/',
                    1: DATA_INPUT_DIR_SERVERS[SERVER] + 'Lebreton_Inputs/',
                    2: DATA_INPUT_DIR_SERVERS[SERVER] + 'LebretonDivision_Inputs/',
                    3: DATA_INPUT_DIR_SERVERS[SERVER] + 'Lebreton_Kaandorp_init_Inputs/',
                    4: DATA_INPUT_DIR_SERVERS[SERVER] + 'Point_Release/',
                    5: DATA_INPUT_DIR_SERVERS[SERVER] + 'Uniform/'}
SCRATCH_DIR = {0: None,
               1: "/storage/scratch/users/vo18e689/"}[SERVER]

if SUBMISSION in ['simulation', 'analysis']:
    # STARTING YEAR, MONTH AND DAY OF THE SIMULATION
    STARTYEAR: int = int(os.environ['STARTYEAR'])
    STARTMONTH: int = int(os.environ['STARTMONTH'])
    STARTDAY: int = int(os.environ['STARTDAY'])

    # LENGTH IN YEARS OF THE SIMULATION
    SIM_LENGTH: int = int(os.environ['SIMLEN'])

########################################################################################################################
#                                                                                                                      #
#                                         Run/restart specific parameters                                              #
#                                                                                                                      #
########################################################################################################################
# DICTIONARY TO CONVERT 0/1 INDICATORS TO BOOLEAN FALSE/TRUE STATEMENTS
BOOLEAN_DICT = {0: False, 1: True}

if SUBMISSION == 'simulation':
    # WHICH OF THE SUBRUNS (DIVISION OF PARTICLE INPUTS)
    RUN: int = int(os.environ['RUN'])

    # RESTART NUMBER OF THE RUN
    # 0 = NEW RUN, N>0 INDICATES NTH YEAR FROM SIMULATION START
    RESTART: int = int(os.environ['RESTARTNUM'])

# LAGRANGIAN OR POSTPROCESSING RUN
# 0 = RUN LAGRANGIAN SIMULATION, 1 = RUN POST PROCESSING
if SUBMISSION == 'simulation':
    POST_PROCESS: bool = BOOLEAN_DICT[int(os.environ['POST_PROCESS'])]
else:
    POST_PROCESS = False

if SUBMISSION in ['simulation', 'analysis']:
    # ENSEMBLE NUMBER
    ENSEMBLE = int(os.environ['ENSEMBLE'])

########################################################################################################################
#                                                                                                                      #
#                                       Advection scenario specific parameters                                         #
#                                                                                                                      #
########################################################################################################################

# UV ADVECTION DATA
ADVECTION_DICT: dict = {0: 'HYCOM_GLOBAL', 1: 'HYCOM_CARIBBEAN', 2: 'CMEMS_MEDITERRANEAN'}
ADVECTION_DATA: str = ADVECTION_DICT[int(os.environ['ADVECTION_DATA'])]
# STOKES DRIFT: 0 -> STOKES, 1 -> NO STOKES
STOKES: int = int(os.environ['STOKES'])

########################################################################################################################
#                                                                                                                      #
#                      Scenario specific parameters, not considering input and advection scenarios                     #
#                                                                                                                      #
########################################################################################################################
# MODEL SCENARIO SETTINGS
SCENARIO_DICT: dict = {0: 'AdvectionDiffusionOnly', 1: 'CoastalProximity', 2: 'Stochastic',
                       3: 'ShoreDependentResuspension', 4: 'TurrellResuspension', 5: 'SizeTransport',
                       6: 'FragmentationKaandorp', 7: 'FragmentationKaandorpPartial'}

SCENARIO_NUM: int = int(os.environ["SCENARIO"])
SCENARIO_NAME: str = SCENARIO_DICT[SCENARIO_NUM]

if SUBMISSION in ['simulation', 'analysis']:
    # FOR 'CoastalProximity': TIME IN BEACHING ZONE PRIOR TO BEACHING (DAYS)
    VICINITY: int = int(os.environ['VICINITY'])

    # BEACHING TIMESCALE
    SHORE_TIME: int = int(os.environ['SHORETIME'])  # days
    # RESUSPENSION TIMESCALE
    RESUS_TIME: int = int(os.environ['RESUSTIME'])  # days
    # FIXED_RESUS OR SIZE DEPENDENT RESUSPENSION TIMESCALE
    FIXED_RESUS: bool = BOOLEAN_DICT[int(os.environ['FIXED_RESUS'])]

    # FOR 'ShoreDependentResuspension': 0 -> SAMARAS ET AL (2015) RESUSPENSION DEPENDENCE,
    #                                   1 -> 1:4 RESUSPENSION DEPENDENCE
    SHORE_DEP: int = int(os.environ['SHOREDEPEN'])

    # FOR 'TurrellResuspension': MINIMUM OFF-SHORE WIND SPEED FOR RESUSPENSION (x10 TO NOT GET DECIMAL IN OUTPUT FILES)
    WMIN: int = int(os.environ['WMIN'])

    # FOR 'FragmentationKaandorpPartial': WHETHER WE WANT TO INCLUDE OCEAN FRAGMENTATION OR NOT
    OCEAN_FRAG: bool = BOOLEAN_DICT[int(os.environ['OCEAN_FRAG'])]


########################################################################################################################
#                                                                                                                      #
#                                       Input scenario specific parameters                                             #
#                                                                                                                      #
########################################################################################################################
# THE INPUT SCENARIO
INPUT_NAMES = {0: 'Jambeck', 1: 'Lebreton', 2: 'LebretonDivision', 3: 'LebretonKaandorpInit', 4: 'Point_Release',
               5: 'Uniform'}
INPUT = INPUT_NAMES[int(os.environ['INPUT'])]
# DIRECTORY CONTAINING INITIAL INPUTS FOR RESTART == 0
INPUT_DIREC = INPUT_DIREC_DICT[int(os.environ['INPUT'])]

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
    # THE NUMBER OF PARTICLES PER RELEASE STEP PER RUN
    INPUT_DIV = 5000
    # NUMBER OF RUNS
    RUN_RANGE: int = 1
    # STARTING LONGITUDE FOR POINT RELEASE (DEGREES EAST)
    INPUT_LON = 0.000000
    # STARTING LATITUDE FOR POINT RELEASE (DEGREES NORTH)
    INPUT_LAT = 38.604168
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
if SUBMISSION in ['analysis']:
    # IF PARALLEL_STEP == 1, THEN WE ARE CALCULATING THE ANALYSIS FOR EACH RUN/RESTART FILE
    # IF PARALLEL_STEP == 2, THEN WE COMBINE ALL THE INDIVIDUAL RUN/RESTART ANALYSIS OUTPUT FILES INTO ONE FOR THE WHOLE
    #                        MODEL SETUP
    PARALLEL_STEP: int = int(os.environ['PARALLEL_STEP'])
    # WHICH OF THE SUBRUNS (DIVISION OF PARTICLE INPUTS)
    RUN: int = int(os.environ['RUN'])
    # RESTART NUMBER OF THE RUN
    # 0 = NEW RUN, N>0 INDICATES NTH YEAR FROM SIMULATION START
    RESTART: int = int(os.environ['RESTARTNUM'])
    # PLASTIC CONCENTRATION
    CONCENTRATION = BOOLEAN_DICT[int(os.environ['CONCENTRATION'])]
    # VERTICAL PLASTIC CONCENTRATION
    VERTICAL_CONCENTRATION = BOOLEAN_DICT[int(os.environ['VERTICAL_CONCENTRATION'])]
    # LON LAT AVERAGED CONCENTRATIONS
    LONLAT_CONCENTRATION = BOOLEAN_DICT[int(os.environ['LONLAT_CONCENTRATION'])]
    # TIMESERIES OF BEACHED/COASTAL FRACTIONS
    TIMESERIES = BOOLEAN_DICT[int(os.environ['TIMESERIES'])]
    # MAX DISTANCE FROM SHORE
    MAX_DISTANCE = BOOLEAN_DICT[int(os.environ['MAX_DISTANCE'])]
    # CREATING TIME SLICES OF THE SIMULATION
    TIMESLICING = BOOLEAN_DICT[int(os.environ['TIMESLICING'])]
    # CALCULATING BASIC STATISTICS FOR EACH PARTICLE TRAJECTORY
    STATISTICS = BOOLEAN_DICT[int(os.environ['STATISTICS'])]
    # CAlCULATING SEPARATION DISTANCES BETWEEN PARTICLES
    SEPARATION = BOOLEAN_DICT[int(os.environ['SEPARATION'])]
    # CALCULATING THE PARTICLE SIZE SPECTRUM/DISTRIBUTION
    SIZE_SPECTRUM = BOOLEAN_DICT[int(os.environ['SIZE_SPECTRUM'])]
    # FRACTION OF PLASTIC ENTERING THE OCEAN IN 2010 THAT IS CONSIDERED BUOYANT (BASED ON GEYERS ET AL., 2017)
    BUOYANT = 0.54
else:
    PARALLEL_STEP: int = 1

########################################################################################################################
#                                                                                                                      #
#                                       Model and physical parameters                                                  #
#                                                                                                                      #
########################################################################################################################
# MODEL INTEGRATION TIMESTEP
TIME_STEP = timedelta(minutes=0.5)
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

if SUBMISSION in ['simulation', 'analysis']:
    # INITIAL PARTICLE SIZE (M)
    INIT_SIZE = int(os.environ['PARTICLE_SIZE']) * SIZE_FACTOR
    # INITIAL DENSITY (KG/M^3): 920 = POLYPROPYLENE, 980 = HIGH DENSITY POLYETHYLENE (BRIGNAC ET AL. 2017)
    INIT_DENSITY = int(os.environ['INIT_DENSITY'])
    # FRAGMENTATION PROBABILITY
    P_FRAG = int(os.environ['P']) * 1e-1
    # NUMBER OF SPATIAL DIMENSIONS
    DN = int(os.environ['DN']) * 1e-1
    # NUMBER OF SIZE CLASSES
    SIZE_CLASS_NUMBER = int(os.environ['SIZE_CLASS_NUMBER'])
    # FRAGMENTATION TIMESCALE (DAYS)
    LAMBDA_FRAG = int(os.environ['LAMBDA_FRAG'])
    # OPEN OCEAN FRAGMENTATION TIMESCALE (DAYS)
    LAMBDA_OCEAN_FRAG = int(os.environ['LAMBDA_OCEAN_FRAG'])
    # CRITICAL SHEAR STRESS FOR RESUSPENSION OF PARTICLES FROM THE SEA BED, HORIZONTAL DIFFUSION AT THE SEA BED
    SEABED_CRIT = int(os.environ['SEABED_CRIT']) * 1E-3

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
#                                       Plotting-specific parameters                                                   #
#                                                                                                                      #
########################################################################################################################
if SUBMISSION in ['visualization']:
    # DEFAULTS TO PREVENT ERRORS
    RUN: int = 0
    RESTART: int = 0
    SHORE_TIME: int = 20
    RESUS_TIME: int = 69
    VICINITY: int = 0
    SHORE_DEP: int = 0
    WMIN: int = 0
    INIT_SIZE: float = 5000 * SIZE_FACTOR
    INIT_DENSITY: int = 920
    STARTYEAR: int = 2010
    STARTMONTH: int = 1
    STARTDAY: int = 1
    SEABED_CRIT: float = 0.025
    P_FRAG: float = 4 * 1e-1
    DN: float = 25 * 1e-1
    SIZE_CLASS_NUMBER: int = 6
    LAMBDA_FRAG: int = 388
    ENSEMBLE: int = 1
    OCEAN_FRAG: bool = False
    LAMBDA_OCEAN_FRAG: int = 388
    FIXED_RESUS: bool = False


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

elif SCENARIO_NAME in ['FragmentationKaandorp', 'FragmentationKaandorpPartial']:
    os.system('echo "The beaching timescale is {} days "'.format(SHORE_TIME))
    os.system('echo "The fragmentation timescale is {} days "'.format(LAMBDA_FRAG))
    os.system('echo "The initial particle size is {} m and a rho of {} kg/m^3"'.format(INIT_SIZE, INIT_DENSITY))
    os.system('echo "This simulation has {} size classes"'.format(SIZE_CLASS_NUMBER))
