# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 09:47:40 2019

@author: Victor Onink

Here is one general file that will contain all the different beaching scenarios
"""
#Loading all the components from parcels
from parcels import ErrorCode
import parcels.rng as rng
#Loading the functions I have for the particle classes, advection, beaching and setup
from beachFunc import DeleteParticle,beachAdvDifOnly,beachVicinity,beachStochastic,beachShoreResus,beachTurrellResus
from transportFunc import AntiBeachNudging,AdvectionRK4_floating,BrownianMotion2D_floating,initialInput
from setupFunc import FileNames,restartInitialPosition,CreatePSET,CreateFieldSet
#And various other packages of use
from datetime import timedelta, datetime
import numpy as np
import os

###############################################################################
# General run setup parameters                                                #
###############################################################################
os.system('echo "Loading the general run parameters"')
#0=First Order, 1 = Coastal proximity, 2 = simple beaching/resuspension,
#3 = coasttype dependence, 4= Turrell (2020)
scenario=int(os.environ['SCENARIO'])
os.system('echo "scenario="'+str(scenario))
#for scenario 1, the time a particle needs to be near the coast to be deleted
vicinity=int(os.environ['VICINITY']) #days
os.system('echo "vicinity="'+str(vicinity))
#for scenario 2, the beaching and resuspension timescales
shoreTime,resusTime=int(os.environ['SHORETIME']),int(os.environ['RESUSTIME']) #days, days
#For scenario 3, how does the coastline dependence work? 
#0 = more sand is less likely beaching, 1 = more land is more likely beaching
shoreDepen = int(os.environ['SHOREDEPEN'])
#for scenario 4, what is the minimum off-shore wind speed for resuspension?
Wmin=int(os.environ['WMIN'])
#for multiple sub runs, which one is being run
run=int(os.environ['RUN'])
os.system('echo "run="'+str(run))
#restart stage, where 0 = a completely new set up runs, and progressively higher
#values indicate the year of the simulation at which new run starts
#e.g. 1 means we start after one year of simulation, 2 after two years, etc.
restart=int(os.environ['RESTARTNUM'])
os.system('echo "restart="'+str(restart))
#starting year of the simulation
startyear=int(os.environ['STARTTIME'])
os.system('echo "starting year="'+str(startyear))
#Which input file we use. if input=0, then we use Jambeck
Input=int(os.environ['INPUT'])
os.system('echo "input="'+str(Input))
#Inclusion of Stokes Drift. 0 = included, 1 = not included
stokes = int(os.environ['STOKES'])
#The option to run multiple ensembles, which we of course don't want saved in the
#same folder since they they would overwrite each other...
ensemble = int(os.environ['ENSEMBLE'])

#To save myself a lot of hassle, we will have a server parameter:
# 0 = kup cluster, 1 = ubelix
server = 1


#%%
os.system('echo "Create the fieldset"')
"""
Loading in all the relevant data into the fieldset
"""
if server==0:
    datadirec="/alphadata04/onink/lagrangian_sim/"
    dataInputdirec="/alphadata04/onink/lagrangian_sim/BeachingSim/Input/"
elif server==1:
    datadirec="/home/ubelix/climate/shared/onink/"
    dataInputdirec="/home/ubelix/climate/shared/onink/Input/"

fieldset=CreateFieldSet(server,stokes,scenario)
#%%
os.system('echo "Setting up the particleset"')
###############################################################################
# First we set the name of the outputfile, and the restart file               #
###############################################################################
if server==0:
    rootodirec="/alphadata04/onink/lagrangian_sim/BeachingSim/Output/"
elif server==1:
    rootodirec="/home/ubelix/climate/shared/onink/Output/"

ofile,rfile=FileNames(rootodirec,scenario,ensemble,startyear,Input,restart,run,
                      vicinity,shoreTime,resusTime,shoreDepen,stokes)

###############################################################################
# Setting the time parameters of the runs, which will depend on restarts, and #
# the initial particle positions                                              #
###############################################################################
starttime=datetime(startyear+restart,1,1,0,0)
endtime=datetime(startyear+restart+1,1,1,0,0)
simulationlength=endtime-starttime

###############################################################################
# Loading in the starting coordinates. In case of a restart file we make use  #
# of the restart function we use, which grabs the particle positions from the #
# previous file                                                               #
###############################################################################
if restart==0:
    if Input==0:
        lons=np.load(dataInputdirec+'Jambeck2010/Jam'+str(2010)+'Lons'+str(run)+'.npy')
        lats=np.load(dataInputdirec+'Jambeck2010/Jam'+str(2010)+'Lats'+str(run)+'.npy')
        weights=np.load(dataInputdirec+'Jambeck2010/Jam'+str(2010)+'Weight'+str(run)+'.npy')
    if Input==1:
        lons=np.load(dataInputdirec+'Lebreton2010/Leb'+str(2010)+'Lons'+str(run)+'.npy')
        lats=np.load(dataInputdirec+'Lebreton2010/Leb'+str(2010)+'Lats'+str(run)+'.npy')
        weights=np.load(dataInputdirec+'Lebreton2010/Leb'+str(2010)+'Weight'+str(run)+'.npy')
    beached=np.zeros(len(lons),dtype=np.int32)
    age_par=np.zeros(len(lons),dtype=np.int32)
    if scenario==1:
        prox_par=np.zeros(len(lons),dtype=np.int32)
    repeatStep=timedelta(days=31)
else:
    if scenario==1:
        lons,lats,beached,age_par,prox_par,weights=restartInitialPosition(rfile,restart,scenario)
    else:
        lons,lats,beached,age_par,weights=restartInitialPosition(rfile,restart,scenario)
    repeatStep=None
lons[lons>180]-=360

###############################################################################
# Creating the particle set                                                   #
###############################################################################

#Creating the particle set, with an initial re-release time of 4 weeks if restart==0
if scenario==1:
    pset = CreatePSET(scenario,lons,lats,beached,age_par,weights,starttime,
                      repeatStep,fieldset,prox_par)
else:
    pset = CreatePSET(scenario,lons,lats,beached,age_par,weights,starttime,
                      repeatStep,fieldset)
    
rng.seed(102334155)

#%% 
os.system('echo "Loading the relevant kernels"')

###############################################################################
# The beaching kernel                                                         #
###############################################################################
if scenario==0:
    beachK=pset.Kernel(beachAdvDifOnly)
elif scenario==1:
    beachK=pset.Kernel(beachVicinity)
elif scenario==2:
    beachK=pset.Kernel(beachStochastic)
elif scenario==3:
    beachK=pset.Kernel(beachShoreResus)
elif scenario==4:
    beachK=pset.Kernel(beachTurrellResus)
    # beachK=pset.Kernel(beachStochastic)
        
###############################################################################
# And now the overall kernel                                                  #
###############################################################################
if scenario==0:
    totalKernel=pset.Kernel(initialInput)+pset.Kernel(AdvectionRK4_floating)+pset.Kernel(BrownianMotion2D_floating)+beachK
else:
    totalKernel=pset.Kernel(initialInput)+pset.Kernel(AdvectionRK4_floating)+pset.Kernel(BrownianMotion2D_floating)+pset.Kernel(AntiBeachNudging)+beachK
    
#%%
###############################################################################
# finally the execution                                                       #
###############################################################################
os.system('echo "The actual calculation of the particle trajectories"')
pfile = pset.ParticleFile(name=ofile,
                          outputdt=timedelta(hours=24))

pset.execute(totalKernel,
              runtime=timedelta(days=simulationlength.days),
              dt=timedelta(minutes=10),
              recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
              output_file=pfile
              )

pfile.export()
