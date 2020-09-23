#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 15:04:13 2020

@author: victoronink

Different implementations of beaching
"""
import math 
from parcels import rng as random


        
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
            if particle.prox>86400*particle.vic:
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
            beach_prob=math.exp(-particle.dt/(particle.LamBeach*86400.))
            if random.random(0.,1.)>beach_prob:
                particle.beach=1
    #Now the part where we build in the resuspension
    elif particle.beach==1:
        resus_prob=math.exp(-particle.dt/(particle.LamResus*86400.))
        if random.random(0.,1.)>resus_prob:
            particle.beach=0
    #Update the age of the particle
    particle.age+=particle.dt

##############################################################################
##############################################################################
##############################################################################
    
def beachShoreResus(particle,fieldset,time):
    if particle.beach==0:
        dist=fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
        #s=sandiness
        s=fieldset.coasttype[time, particle.depth, particle.lat, particle.lon]
        if dist<10:
            # beach_prob=math.exp(-particle.dt/(particle.coastPar*(0.75+0.25*s)*86400.))
            beach_prob=math.exp(-particle.dt/(particle.LamBeach*86400.))
            if random.random(0,1.)>beach_prob:
                particle.beach=1
    #Next the resuspension part
    elif particle.beach==1:
        s=fieldset.coasttype[time, particle.depth, particle.lat, particle.lon]
        resus_prob=math.exp(-particle.dt/(particle.LamResus*(0.75+0.25*s)*86400.))
        if random.random(0,1.)>resus_prob:
            particle.beach=0
    #Update the age of the particle
    particle.age+=particle.dt
