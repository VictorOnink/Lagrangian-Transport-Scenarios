#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 15:04:13 2020

@author: victoronink

Different implementations of beaching
"""
import math 
from parcels import rng as random
import os

scenario=int(os.environ['SCENARIO'])
shoreDepen = int(os.environ['SHOREDEPEN'])
shoreTime,resusTime=int(os.environ['SHORETIME']),int(os.environ['RESUSTIME']) #days, days

###############################################################################
# The delete particle Kernel                                                  #
###############################################################################
def DeleteParticle(particle, fieldset, time):
    #This delete particle format from Philippe Delandmeter
    #https://github.com/OceanParcels/Parcelsv2.0PaperNorthSeaScripts/blob/master/northsea_mp_kernels.py
    print("Particle [%d] lost !! (%g %g %g %g)" % (particle.id, particle.lon, particle.lat, particle.depth, particle.time))
    particle.delete()

###############################################################################
# The beaching kernels                                                        #
###############################################################################
    
def beachAdvDifOnly(particle,fieldset,time):
    #A particle is considered beached if it is within a land cell
    if math.floor(fieldset.landID[time,particle.depth,particle.lat,particle.lon])==1:
        particle.beach=1
    #Update the age of the particle
    particle.age+=particle.dt
        
##############################################################################
##############################################################################
##############################################################################

def beachVicinity(particle,fieldset,time):
    if particle.beach==0:        
        dist=fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
        #If a particle is within 10 km of the shore
        if dist<10:
            particle.prox+=particle.dt
        else:
            particle.prox=0.
        if particle.prox>86400*fieldset.vic:
            particle.beach=1
    #Update the age of the particle
    particle.age+=particle.dt

##############################################################################
##############################################################################
##############################################################################
def beachStochastic(particle,fieldset,time):
    if particle.beach==0:
        dist=fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
        if dist<10: 
            if random.random()>fieldset.p_beach:
                particle.beach=1
    #Now the part where we build in the resuspension
    elif particle.beach==1:
        if random.random()>fieldset.p_resus:
            particle.beach=0
    #Update the age of the particle
    particle.age+=particle.dt

##############################################################################
##############################################################################
##############################################################################
    
def beachShoreResus(particle,fieldset,time):
    if particle.beach==0:
        dist=fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
        if dist<10:
            if random.random()>fieldset.p_beach:
                particle.beach=1
    #Next the resuspension part
    elif particle.beach==1:
        p_resus=fieldset.p_resus[time, particle.depth, particle.lat, particle.lon]
        if random.random()>p_resus:
            particle.beach=0
    #Update the age of the particle
    particle.age+=particle.dt
    
##############################################################################
##############################################################################
##############################################################################

def ProbShore(shoreDepen,scenario,sandy):
    """
    Here we calculate a field for the resuspension timescale, depending on
    which dependence scenario we want to examine. The function is written to
    be flexible, additional scenarios can be added easily and depending on
    future work I can also compute and return spatially varying beaching
    timescales.

    Parameters
    ----------
    shoreDepen : int
        Indicator of which shore dependence scenario we want to run.
    scenario : int
        Which scenario are we running in general.
    sandy : array
        Degree of sandiness for each HYCOM cell.

    Returns
    -------
    p_resus, which is an array with the same size as sandy.

    """
    if scenario==3:
        if shoreDepen==0:
            p_r=resusTime*(0.75+0.25*sandy)
        if shoreDepen==1:
            p_r=resusTime*(0.25+0.75*sandy)
        return p_r
