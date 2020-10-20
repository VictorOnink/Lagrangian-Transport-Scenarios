import os
from dotenv import load_dotenv
from datetime import timedelta

load_dotenv()

SUBMISSION = str(os.environ["SUBMISSION"])
os.system('echo "run="' + SUBMISSION)
# DIRECTORIES
DATA_DIR_SERVERS: dict = {0: "/alphadata04/onink/lagrangian_sim/", 1: "/home/ubelix/climate/shared/onink/"}
DATA_INPUT_DIR_SERVERS: dict = {0: "/alphadata04/onink/lagrangian_sim/BeachingSim/Input/",
                                1: "/home/ubelix/climate/shared/onink/Input/"}
DATA_OUTPUT_DIR_SERVERS: dict = {0: "/alphadata04/onink/lagrangian_sim/BeachingSim/Output/",
                                 1: "/home/ubelix/climate/shared/onink/Output/"}

# SCENARIO SETTINGS
SERVER: int = int(os.environ["SERVER"])
SCENARIO_DICT: dict = {0: 'AdvectionDiffusionOnly', 1: 'CoastalProximity', 2: 'Stochastic',
                       3: 'ShoreDependentResuspension', 4: 'TurrellResuspension'}

SCENARIO_NUM: int = int(os.environ["SCENARIO"])
SCENARIO_NAME: str = SCENARIO_DICT[SCENARIO_NUM]

# for scenario 1, the time a particle needs to be near the coast to be deleted
VICINITY: int = int(os.environ['VICINITY'])

# for scenario 2, the beaching and resuspension timescales
SHORE_TIME: int = int(os.environ['SHORETIME'])  # days
RESUS_TIME: int = int(os.environ['RESUSTIME'])  # days

# for scenario 3, how does the coastline dependence work?
# 0 = more sand is less likely beaching, 1 = more land is more likely beaching
SHORE_DEP: int = int(os.environ['SHOREDEPEN'])

# for scenario 4, what is the minimum off-shore wind speed for resuspension? Divide by 10 to get the actual
# minimum speed in m/s
WMIN: int = int(os.environ['WMIN'])

if SUBMISSION == 'simulation':
    # for multiple sub runs, which one is being run
    RUN: int = int(os.environ['RUN'])

    # restart stage, where 0 = a completely new set up runs, and progressively higher
    # values indicate the year of the simulation at which new run starts
    # e.g. 1 means we start after one year of simulation, 2 after two years, etc.
    RESTART: int = int(os.environ['RESTARTNUM'])

    # Which advection data we want to use.
    ADVECTION_DICT: dict = {0: 'HYCOM_GLOBAL', 1: 'HYCOM_CARIBBEAN'}
    ADVECTION_DATA: str = ADVECTION_DICT[int(os.environ['ADVECTION_DATA'])]

if SUBMISSION == 'analysis':
    # Just default values to not run into errors
    RUN: int = 0
    RESTART: int = 0
    # Now, which analysis do we want to run
    ANALYSIS_DICT = {0: False, 1: True}
    CONCENTRATION = ANALYSIS_DICT[int(os.environ['CONCENTRATION'])]
    TIMESERIES = ANALYSIS_DICT[int(os.environ['TIMESERIES'])]
    MAX_DISTANCE = ANALYSIS_DICT[int(os.environ['MAX_DISTANCE'])]

# starting year of the simulation
START_YEAR: int = int(os.environ['STARTYEAR'])

# How many years does the simulation run?
SIM_LENGTH: int = int(os.environ['SIMLEN'])

# Which input file we use. if input=0, then we use Jambeck
INPUT_NAMES = {0: 'Jambeck', 1: 'Lebreton'}
INPUT_DIREC_DICT = {0: DATA_INPUT_DIR_SERVERS[SERVER] + 'Jambeck_Inputs/',
                    1: DATA_INPUT_DIR_SERVERS[SERVER] + 'Lebreton_Inputs/'}
INPUT = INPUT_NAMES[int(os.environ['INPUT'])]
INPUT_DIREC = INPUT_DIREC_DICT[int(os.environ['INPUT'])]

INPUT_MAX = 1000000.0 # Maximum plastic mass input assigned to one particle in tons
INPUT_MIN = 0.0 # Minimum plastic mass input for a cell in order to be considered for the input
if INPUT == 'Jambeck':
    RUN_RANGE: int = 9
elif INPUT == 'Lebreton':
    RUN_RANGE: int = 4

# Inclusion of Stokes Drift. 0 = included, 1 = not included
STOKES: int = int(os.environ['STOKES'])

# The option to run multiple ensembles, which we of course don't want saved in the
# same folder since they they would overwrite each other...
ENSEMBLE = int(os.environ['ENSEMBLE'])

# Model setup parameters
TIME_STEP = timedelta(minutes=10)  # integration timestep
OUTPUT_TIME_STEP = timedelta(hours=24)  # timestep of saving particle positions to the output file
REPEAT_DT_R0 = timedelta(days=31) # timestep of releasing particles for restart == 0. Otherwise, REPEAT_DT = None
REPEAT_DT_ELSE = None
SEED = 'Fixed'  # Seed for the rng
KH_HOR = 10  # m^2 / s, following Lacerda et al. (2019) and Liubertseva et al. (2018)
COAST_D = 10  # km, the distance from the nearest shoreline that falls under the coastal zone.
BUOYANT = 0.54  # 54% of plastic entering the ocean in 2010 is buoyant (based on Geyers et al., 2017)

# LOG
os.system('echo "scenario = "' + str(SCENARIO_NAME))
if SCENARIO_NAME == 'CoastalProximity':
    os.system('echo "vicinity = "' + str(VICINITY))
elif SCENARIO_NAME == "Stochastic":
    os.system('echo "shoretime = "' + str(SHORE_TIME))
    os.system('echo "resustime = "' + str(RESUS_TIME))
elif SCENARIO_NAME == 'ShoreDependentResuspension':
    os.system('echo "shoretime = "' + str(SHORE_TIME))
    os.system('echo "resustime = "' + str(RESUS_TIME))
    os.system('echo "shoredepen = "' + str(SHORE_DEP))
elif SCENARIO_NAME == 'TurrellResuspension':
    os.system('echo "shoretime = "' + str(SHORE_TIME))
    os.system('echo "Wmin = "' + str(WMIN))
if SUBMISSION == 'simulation':
    os.system('echo "run="' + str(RUN))
    os.system('echo "restart="' + str(RESTART))
os.system('echo "starting year="' + str(START_YEAR))
os.system('echo "input scenario = "' + str(INPUT))
