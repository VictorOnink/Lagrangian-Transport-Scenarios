#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 11:14:04 2020

@author: victoronink

This will be a class that I will use for all my file naming
"""
import os

        
class FileNames(object):
    def __init__(self,sce,vic,st,rt,sd,Wmin,inp,year,stoke,ens,simlen,server=1,
                 restart=0,run=None):
        #Naming the scenario
        scenarionames = {0:'AdvectionDiffusionOnly',1:'CoastalProximity',2:'Stochastic',
                         3:'ShoreDependentResuspension',4:'TurrellResuspension'}
        self.sce = sce
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
        self.inp = inp
        #Number of runs in each input scenario
        if run==None:
            if self.input == 'Jambeck':
                self.run = 8
            elif self.input == 'Lebreton':
                self.run = 3
        else:
            self.run=run
        #Restart position of simulation. 0 = new similation
        self.restart=restart
        #Parameters relevant for all scenarios
        #Starting year of the simulation
        self.startyear = year
        #0 = Stokes included, 1 = Stokes not included
        self.stokes = stoke
        #Ensemble number
        self.ensemble = ens
        #Length of the simulation (years)
        self.simlen = simlen
        #Root output branch
        if server==0:
            self.rootdirec="/alphadata04/onink/lagrangian_sim/BeachingSim/Output/"
        elif server==1:
            self.rootdirec="/home/ubelix/climate/shared/onink/Output/"
        
    def printScenario(self):
        os.system('echo The beaching scenario is '+self.scenario)
        
    def parcelsOutput(self):
        """
        A function that returns the filename of the outputfile and the restart file
        for a parcels run. 
        
        
        Parameters
        ----------
        rootodirec : string
            root directory in which everything is saved.
        scenario : string
            which beaching scenario. 0 = Advection & Diffusion Only
                                     1 = Coastal vicinity
                                     2 = Stochastic
                                     3 = Shore dependent resuspension
                                     4 = Turrell (2020) resuspension
        ensemble : int
            Which ensemble member.
        startyear : int
            Start year of the entire simulation.
        Input : int
            Which input scenario. 0 = Jambeck
                                  1 = Lebreton
        restart : int
            Which restart this is
        run : int
            Which particle subset are we dealing with.
        vicinity : int, optional
            Coastal vicinity cutoff for beaching (days). The default is 0.
        shoreTime : int, optional
            Beaching e-folding timescale (days). The default is 10.
        resusTime : int, optional
            Resuspension e-folding timescale (days). The default is 69.
        shoreDepen : int, optional
            Which shore-dependent resuspension scenario are we considering. The default is 0.
        stokes : int, optional
            Do we include stokes drift. The default is 0. 0 = Include Stokes
                                                          1 = Do not include Stokes
    
        Returns
        -------
        ofile, rfile = Output file name, restart file name
    
        """
        def prefixDeterminant(scenario,stokes):
            components={'AdvectionDiffusionOnly':'AdvDifOnly','CoastalProximity':'Prox','Stochastic':'Stochastic',
                        'ShoreDependentResuspension':'SDResus','TurrellResuspension':'TurrellResuspension',
                        0:'',1:'_NS'}
            return components[scenario]+components[stokes]
        prefix=prefixDeterminant(self.scenario,self.stokes)
        if self.scenario=='AdvectionDiffusionOnly':
            odirec=self.rootdirec+"AdvDifOnly_e_"+str(self.ensemble)+"/"
            ofile=odirec+prefix+"_y="+str(self.startyear)+"_I="+str(self.inp)+"_r="+str(self.restart)+"_run="+str(self.run)+".nc"
            rfile=odirec+prefix+"_y="+str(self.startyear)+"_I="+str(self.inp)+"_r="+str(self.restart-1)+"_run="+str(self.run)+".nc"
        elif self.scenario=='CoastalProximity':
            odirec=self.rootdirec+"coastal_v_"+str(self.vicinity)+"_e_"+str(self.ensemble)+"/"
            ofile=odirec+prefix+"_v="+str(self.vicinity)+"_y="+str(self.startyear)+"_I="+str(self.inp)+"_r="+str(self.restart)+"_run="+str(self.run)+".nc"
            rfile=odirec+prefix+"_v="+str(self.vicinity)+"_y="+str(self.startyear)+"_I="+str(self.inp)+"_r="+str(self.restart-1)+"_run="+str(self.run)+".nc"
        elif self.scenario=='Stochastic':
            odirec=self.rootdirec+"st_"+str(self.shoretime)+"_rs_"+str(self.resustime)+"_e_"+str(self.ensemble)+"/"
            ofile=odirec+prefix+"_st="+str(self.shoretime)+"_rt="+str(self.resustime)+"_y="+str(self.startyear)+"_I="+str(self.inp)+"_r="+str(self.restart)+"_run="+str(self.run)+".nc"
            rfile=odirec+prefix+"_st="+str(self.shoretime)+"_rt="+str(self.resustime)+"_y="+str(self.startyear)+"_I="+str(self.inp)+"_r="+str(self.restart-1)+"_run="+str(self.run)+".nc"
        elif self.scenario=='ShoreDependentResuspension':
            odirec=self.rootdirec+"SDResus_"+str(self.shoreDepen)+"/st_"+str(self.shoretime)+"_rt_"+str(self.resustime)+"_e_"+str(self.ensemble)+"/"
            ofile=odirec+prefix+"_dep="+str(self.shoreDepen)+"_st="+str(self.shoretime)+"_rt="+str(self.resustime)+"_y="+str(self.startyear)+"_I="+str(self.inp)+"_r="+str(self.restart)+"_run="+str(self.run)+".nc"
            rfile=odirec+prefix+"_dep="+str(self.shoreDepen)+"_st="+str(self.shoretime)+"_rt="+str(self.resustime)+"_y="+str(self.startyear)+"_I="+str(self.inp)+"_r="+str(self.restart-1)+"_run="+str(self.run)+".nc"
        elif self.scenario=='TurrellResuspension':
            odirec=self.rootdirec+"Turrell/st_"+str(self.shoreTime)+"_W_"+str(self.Wmin)+"_e_"+str(self.ensemble)+"/"
            ofile=odirec+prefix+"_Wmin="+str(self.Wmin)+"_st="+str(self.shoretime)+"_y="+str(self.startyear)+"_I="+str(self.inp)+"_r="+str(self.restart)+"_run="+str(self.run)+".nc"
            rfile=odirec+prefix+"_Wmin="+str(self.Wmin)+"_st="+str(self.shoretime)+"_y="+str(self.startyear)+"_I="+str(self.inp)+"_r="+str(self.restart-1)+"_run="+str(self.run)+".nc"
        os.system("echo "+ofile)
        return ofile, rfile

    
    