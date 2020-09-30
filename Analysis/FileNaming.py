#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 11:14:04 2020

@author: victoronink

This will be a class that I will use for all my file naming
"""
import os

class parameters:
    """
    This class is for storing all the different run parameters 
    """
    def __init__(self,sce,vic,st,rt,sd,Wmin,inp,year,stoke,ens,simlen):
        #Naming the scenario
        scenarionames = {0:'AdvectionDiffusionOnly',1:'CoastalProximity',2:'Stochastic',
                         3:'ShoreDependentResuspension',4:'TurrellResuspension'}
        self.scenario = scenarionames[sce]
        #Beaching scenario specific parameters
        if self.scenario == 'CoastalProximity':
            #Time near coast before beaching
            self.vicinity = vic
        if self.scenario in ['Stochastic','ShoreDependentResuspension']:
            #Beaching e-folding timescale (days)
            self.shoretime = st
            #Resuspension e-folding timescale (days)
            self.resustime = rt
        if self.scenario == 'ShoreDependentResuspension':
            #Resuspension shore dependence scenario
            self.shoreDepen = sd
        if self.scenario == 'TurrellResuspension':
            #Beaching e-folding timescale (days)
            self.shoretime = st
            #Resuspension windspeed (m/s)
            self.Wmin = Wmin
        #Input scenario specific parameters
        inputnames={0:'Jambeck',1:'Lebreton'}
        self.input = inputnames[inp]
        #Number of runs in each input scenario
        if self.input == 'Jambeck':
            self.run = 8
        elif self.input == 'Lebreton':
            self.run = 3
        #Parameters relevant for all scenarios
        #Starting year of the simulation
        self.startyear = year
        #0 = Stokes included, 1 = Stokes not included
        self.stokes = stoke
        #Ensemble number
        self.ensemble = ens
        #Length of the simulation (years)
        self.simlen = simlen
    
    @classmethod
    def printScenario(self):
        os.system('echo The beaching scenario is '+self.scenario)
        
# class FileNames(parameters):
#     def __init__(self):
#         parameters.__init__(self)
        
#     @classmethod
#     def printScenario(self):
#         os.system('echo The beaching scenario is '+self.scenario)

    
    