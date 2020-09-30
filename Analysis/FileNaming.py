#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 11:14:04 2020

@author: victoronink

This will be a class that I will use for all my file naming
"""
import os

class parameters:
    def __init__(self):
        #Naming the scenario
        scenarionames = {0:'AdvectionDiffusionOnly',1:'CoastalProximity',2:'Stochastic',
                         3:'ShoreDependentResuspension',4:'TurrellResuspension'}
        self.scenario = scenarionames[int(os.environ['SCENARIO'])]
        
        #Beaching scenario specific parameters
        if self.scenario == 'CoastalProximity':
            #Time near coast before beaching
            self.vicinity = int(os.environ['VICINITY'])
        if self.scenario in ['Stochastic','ShoreDependentResuspension']:
            #Beaching e-folding timescale (days)
            self.shoretime = int(os.environ['SHORETIME'])
            #Resuspension e-folding timescale (days)
            self.resustime = int(os.environ['RESUSTIME'])
        if self.scenario == 'ShoreDependentResuspension':
            #Resuspension shore dependence scenario
            self.shoreDepen = int(os.environ['SHOREDEPEN'])
        if self.scenario == 'TurrellResuspension':
            #Beaching e-folding timescale (days)
            self.shoretime = int(os.environ['SHORETIME'])
            #Resuspension windspeed (m/s)
            self.Wmin = int(os.environ['WMIN'])/10
        
        #Input scenario specific parameters
        inputnames={0:'Jambeck',1:'Lebreton'}
        inp=int(os.environ['INPUT'])
        self.input = inputnames[inp]
        #Number of runs in each input scenario
        if self.input == 'Jambeck':
            self.run = 8
        elif self.input == 'Lebreton':
            self.run = 3

        #Parameters relevant for all scenarios
        #Starting year of the simulation
        self.startyear = int(os.environ['STARTYEAR'])
        #0 = Stokes included, 1 = Stokes not included
        self.stokes = int(os.environ['STOKES'])
        #Ensemble number
        self.ensemble = int(os.environ['ENSEMBLE'])
        #Length of the simulation (years)
        self.simlen = int(os.environ['SIMLEN'])
    
    @classmethod
    def printScenario(self):
        os.system('echo The beaching scenario is '+parameters.scenario)
        
# class filenames(parameters):
#     def __init__(self):
#         parameters.__init__(self)