import os
from dotenv import load_dotenv
from datetime import timedelta

load_dotenv()


# DIRECTORIES
DATA_DIR_SERVERS: dict = {0: "/alphadata04/onink/lagrangian_sim/", 1: "/home/ubelix/climate/shared/onink/"}
DATA_INPUT_DIR_SERVERS: dict = {0: "/alphadata04/onink/lagrangian_sim/BeachingSim/Input/",
                                         1: "/home/ubelix/climate/shared/onink/Input/"}



# SCENARIO SETTINGS
SERVER: int = int(os.environ["SERVER"])
SCENARIO_DICT: dict = {0: 'AdvectionDiffusionOnly', 1: 'CoastalProximity', 2: 'Stochastic', 3: 'ShoreDependentResuspension', 4: 'TurrellResuspension'}

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

# for scenario 4, what is the minimum off-shore wind speed for resuspension?
WMIN: int = int(os.environ['WMIN'])

# for multiple sub runs, which one is being run
RUN: int = int(os.environ['RUN'])

# restart stage, where 0 = a completely new set up runs, and progressively higher
# values indicate the year of the simulation at which new run starts
# e.g. 1 means we start after one year of simulation, 2 after two years, etc.
RESTART: int = int(os.environ['RESTARTNUM'])

# starting year of the simulation
START_YEAR: int = int(os.environ['STARTYEAR'])

# Which input file we use. if input=0, then we use Jambeck
INPUT_NAMES = {0: 'Jambeck', 1: 'Lebreton'}
INPUT = INPUT_NAMES[int(os.environ['INPUT'])]

# Inclusion of Stokes Drift. 0 = included, 1 = not included
STOKES: int = int(os.environ['STOKES'])

# The option to run multiple ensembles, which we of course don't want saved in the
# same folder since they they would overwrite each other...
ENSEMBLE = int(os.environ['ENSEMBLE'])

#The integration timestep of the model
TIME_STEP = timedelta(minutes=10)
OUTPUT_TIME_STEP = timedelta(hours=24)
SEED = 'Fixed'
KH_HOR = 10 #m^2 / s, following Lacerda et al. (2019) and Liubertseva et al. (2018)

# LOG
os.system('echo "scenario = "'+str(SCENARIO_NAME))
if SCENARIO_NAME=='CoastalProximity':
    os.system('echo "vicinity = "' + str(VICINITY))
elif SCENARIO_NAME=="Stochastic":
    os.system('echo "shoretime = "' + str(SHORE_TIME))
    os.system('echo "resustime = "' + str(RESUS_TIME))
elif SCENARIO_NAME=='ShoreDependentResuspension':
    os.system('echo "shoretime = "' + str(SHORE_TIME))
    os.system('echo "resustime = "' + str(RESUS_TIME))
    os.system('echo "shoredepen = "' + str(SHORE_DEP))
elif SCENARIO_NAME=='TurrellResuspension':
    os.system('echo "shoretime = "' + str(SHORE_TIME))
    os.system('echo "Wmin = "' + str(WMIN))
os.system('echo "run="'+str(RUN))
os.system('echo "restart="'+str(RESTART))
os.system('echo "starting year="'+str(START_YEAR))
os.system('echo "input scenario = "'+str(INPUT))
