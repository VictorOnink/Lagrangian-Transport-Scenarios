#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 15:23:16 2020

@author: victoronink

Plastic transport kernels for floating plastic
"""
import math 
from parcels import rng as random

##############################################################################
##############################################################################
##############################################################################

def AntiBeachNudging(particle,fieldset,time):
    """    
    If a particle is within 0.5km of the nearest coastline (determined by sampling
    the distance2shore field), then it gets nudged away back out to sea. I have
    fields for border currents, so that the particle gets nudged in the right 
    direction with a speed of 1 - 1.414 m s^{-1}. 
    
    With dt=10 minutes a particle gets displaced by 600 - 848 m back out to sea.
    """
    d1=particle.depth
    if fieldset.distance2shore[time,d1,particle.lat,particle.lon]<0.5:
        borUab,borVab=fieldset.borU[time, d1, particle.lat, particle.lon],fieldset.borV[time, d1, particle.lat, particle.lon]
        particle.lon-=borUab*particle.dt
        particle.lat-=borVab*particle.dt

##############################################################################
##############################################################################
##############################################################################

def AdvectionRK4_floating(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration.

    Function needs to be converted to Kernel object before execution
        
    A particle only moves if it has not beached (rather obviously)
    """
    if particle.beach==0:
        particle.distance=fieldset.distance2shore[time,particle.depth,particle.lat,particle.lon]
        d2=particle.depth
        if particle.lon>180:
            particle.lon-=360
        if particle.lon<-180:
            particle.lon+=360
        (u1, v1) = fieldset.UV[time, d2, particle.lat, particle.lon]
        (uS1, vS1) = fieldset.Ust[time, d2, particle.lat, particle.lon],fieldset.Vst[time, d2, particle.lat, particle.lon]
        # lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
        lon1, lat1 = (particle.lon + (u1+uS1)*.5*particle.dt, particle.lat + (v1+vS1)*.5*particle.dt)

        if lon1>180:
            lon1-=360
        if lon1<-180:
            lon1+=360
        (u2, v2) = fieldset.UV[time + .5 * particle.dt, d2, lat1, lon1]
        (uS2, vS2) = fieldset.Ust[time + .5 * particle.dt, d2, lat1, lon1],fieldset.Vst[time + .5 * particle.dt, d2, lat1, lon1]
        # lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
        lon2, lat2 = (particle.lon + (u2+uS2)*.5*particle.dt, particle.lat + (v2+vS2)*.5*particle.dt)

        if lon2>180:
            lon2-=360
        if lon2<-180:
            lon2+=360
        (u3, v3) = fieldset.UV[time + .5 * particle.dt, d2, lat2, lon2]
        (uS3, vS3) = fieldset.Ust[time + .5 * particle.dt, d2, lat2, lon2],fieldset.Vst[time + .5 * particle.dt, d2, lat2, lon2]
        # lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
        lon3, lat3 = (particle.lon + (u3+uS3)*particle.dt, particle.lat + (v3+vS3)*particle.dt)
        
        if lon3>180:
            lon3-=360
        if lon3<-180:
            lon3+=360
        (u4, v4) = fieldset.UV[time + particle.dt, d2, lat3, lon3]
        (uS4, vS4) = fieldset.Ust[time + particle.dt, d2, lat3, lon3],fieldset.Vst[time + particle.dt, d2, lat3, lon3]
        
        # particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lon += ((u1+uS1) + 2*(u2+uS2) + 2*(u3+uS3) + (u4+uS4)) / 6. * particle.dt
        if particle.lon>180:
            particle.lon-=360
        if particle.lon<-180:
            particle.lon+=360
        # particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle.lat += ((v1+vS1) + 2*(v2+vS2) + 2*(v3+vS3) + (v4+vS4)) / 6. * particle.dt

##############################################################################
##############################################################################
##############################################################################

def BrownianMotion2D_floating(particle, fieldset, time):
    """Kernel for simple Brownian particle diffusion in zonal and meridional direction.
    Assumes that fieldset has fields Kh_zonal and Kh_meridional
    we don't want particles to jump on land and thereby beach"""
    if particle.beach==0:
        # Wiener increment with zero mean and std of sqrt(dt)
        dWx = random.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
        dWy = random.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
    
        bx = math.sqrt(2 * fieldset.Kh_zonal[time, particle.depth, particle.lat, particle.lon])
        by = math.sqrt(2 * fieldset.Kh_meridional[time, particle.depth, particle.lat, particle.lon])
    
        particle.lon += bx * dWx
        particle.lat += by * dWy


##############################################################################
##############################################################################
##############################################################################

def initialInput(particle,fieldset,time):
    """
    Since we have many instances that particles start at the very same position,
    when a particle is first added to the simulation it will get a random kick
    that moves it slightly away from the initial position, and so with multiple
    particles we will see a sort of star pattern around the central position.
    However, the particle shouldn't be put on land though, so a particle will
    have 100000 attempts to be placed in the simulation within a cell that is not
    land. If it doesn't work after 100000 attempts, then the particle just ends up
    starting from the unchanged position. We also set it so that the initial
    position of the particle is always closer to land than the original position
    
    Note: Tests show that at most we get at most around 7 or 8 particles per 
    release getting placed immediately on land. this varies a bit
    
    """
    if particle.age==0:
        #The low latitudes/equatorial regions have larger grid sizes
        if math.fabs(particle.lat)<40.:
            check=0
            distCur=fieldset.distance2shore[time,particle.depth,particle.lat,particle.lon]
            while check<100000:
                potLat=particle.lat+random.uniform(-0.08, 0.08)
                potLon=particle.lon+random.uniform(-0.08, 0.08)
                potLand=math.floor(fieldset.landID[time,particle.depth,particle.lat,particle.lon])
                distPot=fieldset.distance2shore[time,particle.depth,potLat,potLon]
                if potLand==0 and distPot<=distCur:
                    check+=100001
                    particle.lat=potLat
                    particle.lon=potLon
                check+=1
        #Higher latitudes above 40 degrees
        else:
            check=0
            distCur=fieldset.distance2shore[time,particle.depth,particle.lat,particle.lon]
            while check<100000:
                potLat=particle.lat+random.uniform(-0.04, 0.04)
                potLon=particle.lon+random.uniform(-0.04, 0.04)
                potLand=math.floor(fieldset.landID[time,particle.depth,particle.lat,particle.lon])
                distPot=fieldset.distance2shore[time,particle.depth,potLat,potLon]
                if potLand==0 and distPot<=distCur:
                    check+=100001
                    particle.lat=potLat
                    particle.lon=potLon
                check+=1

