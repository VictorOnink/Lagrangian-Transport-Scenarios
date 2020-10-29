import os
from dotenv import load_dotenv
from datetime import timedelta
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
DATA_DIR_SERVERS: dict = {0: "/alphadata04/onink/lagrangian_sim/", 1: "/home/ubelix/climate/shared/onink/"}
DATA_INPUT_DIR_SERVERS: dict = {0: "/alphadata04/onink/lagrangian_sim/BeachingSim/Input/",
                                1: "/home/ubelix/climate/shared/onink/Input/"}
DATA_OUTPUT_DIR_SERVERS: dict = {0: "/alphadata04/onink/lagrangian_sim/BeachingSim/Output/",
                                 1: "/home/ubelix/climate/shared/onink/Output/"}
INPUT_DIREC_DICT = {0: DATA_INPUT_DIR_SERVERS[SERVER] + 'Jambeck_Inputs/',
                    1: DATA_INPUT_DIR_SERVERS[SERVER] + 'Lebreton_Inputs/'}


# STARTING YEAR OF THE SIMULATION
START_YEAR: int = int(os.environ['STARTYEAR'])

# LENGTH IN YEARS OF THE SIMULATION
SIM_LENGTH: int = int(os.environ['SIMLEN'])

########################################################################################################################
#                                                                                                                      #
#                                         Run/restart specific parameters                                              #
#                                                                                                                      #
########################################################################################################################
if SUBMISSION == 'simulation':
    # WHICH OF THE SUBRUNS (DIVISION OF PARTICLE INPUTS
    RUN: int = int(os.environ['RUN'])

    # RESTART NUMBER OF THE RUN
    # 0 = NEW RUN, N>0 INDICATES NTH YEAR FROM SIMULATION START
    RESTART: int = int(os.environ['RESTARTNUM'])

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
                       3: 'ShoreDependentResuspension', 4: 'TurrellResuspension', 5: 'FragmentationCozar'}

SCENARIO_NUM: int = int(os.environ["SCENARIO"])
SCENARIO_NAME: str = SCENARIO_DICT[SCENARIO_NUM]

if SCENARIO_NAME == 'CoastalProximity':
    # TIME IN BEACHING ZONE PRIOR TO BEACHING (DAYS)
    VICINITY: int = int(os.environ['VICINITY'])

if SCENARIO_NAME == 'Stochastic' or SCENARIO_NAME == 'ShoreDependentResuspension':
    # BEACHING TIMESCALE
    SHORE_TIME: int = int(os.environ['SHORETIME'])  # days
    # RESUSPENSION TIMESCALE
    RESUS_TIME: int = int(os.environ['RESUSTIME'])  # days

if SCENARIO_NAME == 'ShoreDependentResuspension':
    # 0 -> SAMARAS ET AL (2015) RESUSPENSION DEPENDENCE, 1 -> 1:4 RESUSPENSION DEPENDENCE
    SHORE_DEP: int = int(os.environ['SHOREDEPEN'])

if SCENARIO_NAME == 'TurrellResuspension':
    # BEACHING TIMESCALE
    SHORE_TIME: int = int(os.environ['SHORETIME'])  # days
    # MINIMUM OFF-SHORE WIND SPEED FOR RESUSPENSION (x10 TO NOT GET DECIMAL IN OUTPUT FILES)
    WMIN: int = int(os.environ['WMIN'])

if SCENARIO_NAME == 'FragmentationCozar':
    # BEACHING TIMESCALE
    SHORE_TIME: int = int(os.environ['SHORETIME'])  # days
    # RESUSPENSION TIMESCALE
    RESUS_TIME: int = int(os.environ['RESUSTIME'])  # days

########################################################################################################################
#                                                                                                                      #
#                                       Input scenario specific parameters                                             #
#                                                                                                                      #
########################################################################################################################
# THE INPUT SCENARIO
INPUT_NAMES = {0: 'Jambeck', 1: 'Lebreton', 2: 'Point release'}
INPUT = INPUT_NAMES[int(os.environ['INPUT'])]
# DIRECTORY CONTAINING INITIAL INPUTS FOR RESTART == 0
INPUT_DIREC = INPUT_DIREC_DICT[int(os.environ['INPUT'])]
# THE NUMBER OF PARTICLES PER RELEASE STEP PER RUN
INPUT_DIV = 5000 # The number of particles per release step per run

if INPUT == 'Jambeck':
    # NUMBER OF RUNS
    RUN_RANGE: int = 9
    # MAXIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MAX = 10.0
    # MINIMUM PLASTIC MASS INPUT ASSIGNED TO ONE PARTICLE (TONS)
    INPUT_MIN = 0.08  # Minimum plastic mass input for a cell in order to be considered for the input
elif INPUT == 'Lebreton':
    RUN_RANGE: int = 4
    INPUT_MAX = 10.0  # Maximum plastic mass input assigned to one particle in tons
    INPUT_MIN = 0.0  # Minimum plastic mass input for a cell in order to be considered for the input


########################################################################################################################
#                                                                                                                      #
#                                       Analysis specific parameters                                                   #
#                                                                                                                      #
########################################################################################################################
if SUBMISSION == 'analysis':
    # DEFAULTS TO PREVENT ERRORS
    RUN: int = 0
    RESTART: int = 0
    # SELECTING ANALYSIS TO RUN
    ANALYSIS_DICT = {0: False, 1: True}
    # PLASTIC CONCENTRATION
    CONCENTRATION = ANALYSIS_DICT[int(os.environ['CONCENTRATION'])]
    # TIMESERIES OF BEACHED/COASTAL FRACTIONS
    TIMESERIES = ANALYSIS_DICT[int(os.environ['TIMESERIES'])]
    # MAX DISTANCE FROM SHORE
    MAX_DISTANCE = ANALYSIS_DICT[int(os.environ['MAX_DISTANCE'])]
    # FRACTION OF PLASTIC ENTERING THE OCEAN IN 2010 THAT IS CONSIDERED BUOYANT (BASED ON GEYERS ET AL., 2017)
    BUOYANT = 0.54

########################################################################################################################
#                                                                                                                      #
#                                       Model and physical parameters                                                  #
#                                                                                                                      #
########################################################################################################################
# MODEL INTEGRATION TIMESTEP
TIME_STEP = timedelta(minutes=10)  # integration timestep
# MODEL OUTPUT TIMESTEP
OUTPUT_TIME_STEP = timedelta(hours=24)
# PARTICLE RELEASE TIMESTEP
if SCENARIO_NAME != 'FragmentationCozar':
    REPEAT_DT_R0 = timedelta(days=31) # timestep of releasing particles for restart == 0. Otherwise, REPEAT_DT = None
else:
    REPEAT_DT_R0 = None
REPEAT_DT_ELSE = None
# RNG SEED VALUE
SEED = 'Fixed'
# HORIZONTAL DIFFUSIVITY M^2 / S, FOLLOWING LACERDA ET AL. (2019) AND LIUBERTSEVA ET AL. (2018)
KH_HOR = 10
# WIDTH OF THE BEACHING ZONE (KM) WITHIN WHICH BEACHING CAN OCCUR
if SCENARIO_NAME != 'AdvectionDiffusionOnly':
    COAST_D = 10  # km, the distance from the nearest shoreline that falls under the coastal zone.
if SCENARIO_NAME == 'FragmentationCozar':
    # INITIAL PARTICLE SIZE (m): THE HEIGHT OF A STANDARD COKE BOTTLE
    INIT_SIZE = 0.35 # m, the height of a standard coke bottle
    # INITIAL DENSITY (KG/M^3): 920 = polypropylene
    INIT_DENSITY = 920
########################################################################################################################
#                                                                                                                      #
#                                                     LOG                                                              #
#                                                                                                                      #
########################################################################################################################
os.system('echo "This is the {} for the {} scenario"'.format(SUBMISSION, SCENARIO_NAME))
os.system('echo "The input scenario is {}"'.format(INPUT))
if SUBMISSION == 'simulation':
    os.system('echo "The starting year is {}, and this is run {}, restart {}"'.format(START_YEAR, RUN, RESTART))
    os.system('echo "We are using {} for the plastic advection"'.format(ADVECTION_DATA))

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

elif SCENARIO_NAME == 'FragmentationCozar':
    os.system('echo "The beaching timescale is {} days "'.format(SHORE_TIME))
    os.system('echo "The resuspension timescale is {} days "'.format(RESUS_TIME))
