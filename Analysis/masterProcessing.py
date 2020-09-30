#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:53:54 2020

@author: victoronink
"""
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
