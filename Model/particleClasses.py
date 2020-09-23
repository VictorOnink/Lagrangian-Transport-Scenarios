#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 15:09:49 2020

@author: victoronink

Creating all the different particle classes
"""
from parcels import FieldSet, Variable, JITParticle
import numpy as np
from operator import attrgetter
import os

shoreTime,resusTime=int(os.environ['SHORETIME']),int(os.environ['RESUSTIME']) #days, days
vicinity=int(os.environ['VICINITY']) #days

class AdvDifParticle(JITParticle):
    #First we keep track of whether a particle is on a land cell or not
    on_land = Variable('on_land', dtype=np.int32, to_write=False,
                        initial=0)
    #Now the beaching variables
    #0=open ocean, 1=beached
    beach    = Variable('beach', dtype=np.int32,
                        initial=attrgetter('beach'))
    #Finally, I want to keep track of the age of the particle
    age = Variable('age',dtype=np.int32,initial=attrgetter('age'))
    #Weight of the particle in tons
    weights= Variable('weights',dtype=np.float32,initial=attrgetter('weights'))
    #Distance of the particle to the coast
    distance= Variable('distance',dtype=np.float32,initial=0)
    
class vicinityParticle(JITParticle):
    #First we keep track of how long a particle has been close to the shore
    prox = Variable('prox', dtype=np.int32, initial=attrgetter('prox'))
    #And we need to define the cutoff point for how long the particle can 
    #be close to shore
    vic = Variable('vic', dtype=np.float32,initial=vicinity, to_write=False)
    #Now the beaching variables
    #0=open ocean, 1=beached
    beach    = Variable('beach', dtype=np.int32,
                        initial=attrgetter('beach'))
    #Finally, I want to keep track of the age of the particle
    age = Variable('age',dtype=np.int32,initial=attrgetter('age'))
    #Weight of the particle in tons
    weights = Variable('weights',dtype=np.float32,initial=attrgetter('weights'))
    #Distance of the particle to the coast
    distance= Variable('distance',dtype=np.float32,initial=0)
    
class StochasticParticle(JITParticle):
    #Now the beaching variables
    #0=open ocean, 1=beached
    beach    = Variable('beach', dtype=np.int32,
                        initial=attrgetter('beach'))
    #Now the setting of the resuspension time and beaching time
    LamResus = Variable('LamResus',dtype=np.float32,
                        initial=resusTime,to_write=False)
    LamBeach = Variable('LamBeach',dtype=np.float32,
                        initial=shoreTime,to_write=False)
    #Finally, I want to keep track of the age of the particle
    age = Variable('age',dtype=np.int32,initial=attrgetter('age'))
    #Weight of the particle in tons
    weights = Variable('weights',dtype=np.float32,initial=attrgetter('weights'))
    #Distance of the particle to the coast
    distance= Variable('distance',dtype=np.float32,initial=0)
