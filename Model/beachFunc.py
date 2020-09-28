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
    print("Particle [%d]" % (particle.id))
    if particle.beach==0:
        dist=fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
        if dist<10:
            if random.random()>fieldset.p_beach:
                particle.beach=1
    #Next the resuspension part
    elif particle.beach==1:
        if random.random()>fieldset.p_resus[time, particle.depth, particle.lat, particle.lon]:
            particle.beach=0
    #Update the age of the particle
    particle.age+=particle.dt
    
##############################################################################
##############################################################################
##############################################################################

def beachTurrellResus(particle,fieldset,time):
    """
    Beaching is implemented the same way as in beachStochastic and beachShoreResus.
    
    Resuspension is based on Turrell 2018 & 2020. Resuspension is possible when
    water levels are at the same level as that of the beached particle. Then,
    only when the offshore wind component is greater than the threshold Wmin
    will the particle actually be resuspended
    """
    # print("Particle [%d]" % (particle.id))
    # if particle.beach==0:
    #     dist=fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
    #     if dist<10:
    #         if random.random()>fieldset.p_beach:
    #             particle.beach=1
    # #Next the resuspension part
    # elif particle.beach==1:
    #     if random.random()>fieldset.p_resus[time, particle.depth, particle.lat, particle.lon]:
    #         particle.beach=0
    #Beaching
    if particle.beach==0:
        dist=fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
        if dist<10:
            if random.random()>fieldset.p_beach:
                particle.beach=1
    #             # particle.depth=fieldset.eta[t,d,la,lo]
    #Resuspension
    # elif particle.beach==1:
    #     sea_elev=fieldset.eta[t,d,la,lo]
    #     #If particles are beached above sea level, then they will remain there
    #     if particle.depth<sea_elev:
    #         #particles will get pushed up to the present water level
    #         particle.depth=sea_elev
    #         #Now, we need to get the offshore wind component
    #         bU,bV=fieldset.borU[t,d,la,lo]*-1*1852*60*math.cos(la*math.pi/180),fieldset.borV[t,d,la,lo]*-1*1852*60
    #         wU,wV=fieldset.U[t,d,la,lo]*1852*60*math.cos(la*math.pi/180),fieldset.V[t,d,la,lo]*1852*60
    #         #magnitude of the b and w vectors, and then dot product between them
    #         mB,mW=math.sqrt(bU**2+bV**2),math.sqrt(wU**2+wV**2)
    #         dot=bU*wU+bV*wV
    #         #Angle between b and w
    #         alpha=math.acos(dot/(mB*mW))
    #         #If the angle between thexs wind and border current is <90 degrees,
    #         #then the wind is directed offshore
    #         if alpha<math.pi/2:
    #             #If the offshore component of wind is greater than Wmin, then
    #             #the particle gets resuspended
    #             if mW*math.cos(alpha)<fieldset.Wmin:
    #                 particle.beach=0
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
            p_resus=resusTime*(0.75+0.25*sandy)
        if shoreDepen==1:
            p_resus=resusTime*(0.25+0.75*sandy)
        return p_resus
