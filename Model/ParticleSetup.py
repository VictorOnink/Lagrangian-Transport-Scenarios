#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 18:25:58 2020

@author: victoronink
"""
from netCDF4 import Dataset
import numpy as np
from parcels import ParticleSet
from ParticleClasses import AdvDifParticle,vicinityParticle,StochasticParticle
from ParticleBeaching import beachAdvDifOnly,beachVicinity,beachStochastic,beachShoreResus,beachTurrellResus
from ParticleTransport import AntiBeachNudging,AdvectionRK4_floating,BrownianMotion2D_floating,initialInput
import os

scenarionames={0:'AdvectionDiffusionOnly',1:'CoastalProximity',2:'Stochastic',
               3:'ShoreDependentResuspension',4:'TurrellResuspension'}
scenario=scenarionames[int(os.environ['SCENARIO'])]
vicinity=int(os.environ['VICINITY']) #days
shoreTime,resusTime=int(os.environ['SHORETIME']),int(os.environ['RESUSTIME']) #days, days
shoreDepen = int(os.environ['SHOREDEPEN'])
Wmin=int(os.environ['WMIN'])

def FileNames(rootodirec,scenario,ensemble,startyear,Input,restart,run,Wmin,
              vicinity=0,shoreTime=10,resusTime=69,shoreDepen=0,stokes=0):
    """
    A function that returns the filename of the outputfile and the restart file
    for a parcels run. 
    
    
    Parameters
    ----------
    rootodirec : string
        root directory in which everything is saved.
    scenario : int
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
    prefix=prefixDeterminant(scenario,stokes)
    if scenario=='AdvectionDiffusionOnly':
        odirec=rootodirec+"AdvDifOnly_e_"+str(ensemble)+"/"
        ofile=odirec+prefix+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart)+"_run="+str(run)+".nc"
        rfile=odirec+prefix+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart-1)+"_run="+str(run)+".nc"
    elif scenario=='CoastalProximity':
        odirec=rootodirec+"coastal_v_"+str(vicinity)+"_e_"+str(ensemble)+"/"
        ofile=odirec+prefix+"_v="+str(vicinity)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart)+"_run="+str(run)+".nc"
        rfile=odirec+prefix+"_v="+str(vicinity)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart-1)+"_run="+str(run)+".nc"
    elif scenario=='Stochastic':
        odirec=rootodirec+"st_"+str(shoreTime)+"_rs_"+str(resusTime)+"_e_"+str(ensemble)+"/"
        ofile=odirec+prefix+"_st="+str(shoreTime)+"_rt="+str(resusTime)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart)+"_run="+str(run)+".nc"
        rfile=odirec+prefix+"_st="+str(shoreTime)+"_rt="+str(resusTime)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart-1)+"_run="+str(run)+".nc"
    elif scenario=='ShoreDependentResuspension':
        odirec=rootodirec+"SDResus_"+str(shoreDepen)+"/st_"+str(shoreTime)+"_rt_"+str(resusTime)+"_e_"+str(ensemble)+"/"
        ofile=odirec+prefix+"_dep="+str(shoreDepen)+"_st="+str(shoreTime)+"_rt="+str(resusTime)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart)+"_run="+str(run)+".nc"
        rfile=odirec+prefix+"_dep="+str(shoreDepen)+"_st="+str(shoreTime)+"_rt="+str(resusTime)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart-1)+"_run="+str(run)+".nc"
    elif scenario=='TurrellResuspension':
        odirec=rootodirec+"Turrell/st_"+str(shoreTime)+"_W_"+str(Wmin)+"_e_"+str(ensemble)+"/"
        ofile=odirec+prefix+"_Wmin="+str(Wmin)+"_st="+str(shoreTime)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart)+"_run="+str(run)+".nc"
        rfile=odirec+prefix+"_Wmin="+str(Wmin)+"_st="+str(shoreTime)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart-1)+"_run="+str(run)+".nc"
    os.system("echo "+ofile)
    return ofile, rfile

##############################################################################
##############################################################################
##############################################################################
def removeNan(dataset,variable,lastSelec,finalTime,lastTimeSelec):
    """
    This function inputs a dataset object for the rfile. We then take the last
    non-masked value for each row, which we then use to initialise the new ofile 
    run. However, in some cases a particle has been deleted during the rfile run,
    and while that particle does stay deleted, we do need to incorporate it in
    the ofile run so that we don't misalign the rows. Therefore, in cases where
    a particle is deleted, we return varSelec with 2 for all those particles,
    since particles where particle.beach==2 will not be advected or be resuspended.
    
    Parameters
    ----------
    dataset : netcdf4 dataset object
        the dataset object we get the variable field from.
    variable : string
        name of the variable we wish to examine.
    lastSelec : int
        For each row, the index of the last unmasked cell. If there are no
        masked cells, it indicates the last cell of the .
    finalTime : datetime object
        The last timestep of the previous restart file.
    lastTimeSelec : datetime object
        the time of the last unmasked call, as indicated by the index from 
        lastSelec.

    Returns
    -------
    varSelec : array Nx1
        The restart array for the given variable to start up the ofile run.

    """
    var=np.array(dataset.variables[variable][:])
    varSelec=var[lastSelec[0],lastSelec[1]]
    varSelec[lastTimeSelec!=finalTime]=2
    return varSelec

##############################################################################
##############################################################################
##############################################################################
def restartInitialPosition(rfile,restart,scenario):
    """
    Parameters
    ----------
    rfile : string
        restart file name.
    restart : int
        Which restart this is
    scenario : int
        which beaching scenario. 0 = Advection & Diffusion Only
                                 1 = Coastal vicinity
                                 2 = Stochastic
                                 3 = Shore dependent resuspension

    Returns
    -------
    TYPE
        lon_reset,lat_reset,beach_reset,age_reset,prox_reset,weight_reset
        restart arrays for longitude, latitude, beached status, particle weight,
        particle age, particle proximity to shore (for scenario 1 only)

    """
    dataset=Dataset(rfile)
    time=dataset.variables['time'][:]
    finalTime=time[0,-1]
    lastSelec=np.ma.notmasked_edges(time,axis=1)[1]
    lastTimeSelec=time[lastSelec[0],lastSelec[1]]

    lon_reset=removeNan(dataset,'lon',lastSelec,finalTime,lastTimeSelec)
    lat_reset=removeNan(dataset,'lat',lastSelec,finalTime,lastTimeSelec)
    beach_reset=removeNan(dataset,'beach',lastSelec,finalTime,lastTimeSelec)
    age_reset=removeNan(dataset,'age',lastSelec,finalTime,lastTimeSelec)
    weight_reset=removeNan(dataset,'weights',lastSelec,finalTime,lastTimeSelec)
    if scenario=='CoastalProximity':
        prox_reset=removeNan(dataset,'prox',lastSelec,finalTime,lastTimeSelec)
    #Returning the relevant arrays    
    if scenario=='CoastalProximity':
        return lon_reset,lat_reset,beach_reset,age_reset,prox_reset,weight_reset
    else:
        return lon_reset,lat_reset,beach_reset,age_reset,weight_reset

##############################################################################
##############################################################################
##############################################################################
def CreatePSET(scenario,lons,lats,beached,age_par,weights,starttime,
               repeatStep,fieldset,prox_par=0):
    """
    Creating the

    Parameters
    ----------
    scenario : int
        which beaching scenario. 0 = Advection & Diffusion Only
                                 1 = Coastal vicinity
                                 2 = Stochastic
                                 3 = Shore dependent resuspension
                                 4 = Turrell (2020) resuspension
    lons : array Nx1
        starting longitudes.
    lats : array Nx1
        starting latitude.
    beached : array Nx1
        started beached status. 0 = Afloat
                                1 = Beached
                                2 = Removed from simulation
    age_par : array Nx1
        particle age (s).
    weights : array Nx1
        particle weight (tons).
    prox_par : array Nx1, optional
        time spent in vicinity of coast (s). The default is 0.
    starttime : datetime
        starting date of the simulation.
    repeatStep : timedelta
        Period with which particles are released. None = particles only released
        at the beginning of the simulation
    fieldset : parcels fieldset
        fieldset containing all the relevant fields.

    Returns
    -------
    pset.

    """
    pclass={'AdvectionDiffusionOnly':AdvDifParticle,'CoastalProximity':vicinityParticle,
            'Stochastic':StochasticParticle,'ShoreDependentResuspension':StochasticParticle,
            'TurrellResuspension':StochasticParticle}
    
    if scenario=='CoastalProximity':
        pset = ParticleSet(fieldset=fieldset, pclass=pclass[scenario], 
                           lon=lons, lat=lats,beach=beached,age=age_par,
                           prox=prox_par,weights=weights,
                           time=starttime, repeatdt=repeatStep)
    else:
        pset = ParticleSet(fieldset=fieldset, pclass=pclass[scenario], 
                           lon=lons, lat=lats,beach=beached,age=age_par,
                           weights=weights,
                           time=starttime, repeatdt=repeatStep)
    return pset

##############################################################################
##############################################################################
##############################################################################
def ParticleBehavior(pset,scenario):
    beachK={'AdvectionDiffusionOnly':beachAdvDifOnly,
            'CoastalProximity':beachVicinity,
            'Stochastic':beachStochastic,
            'ShoreDependentResuspension':beachShoreResus,
            'TurrellResuspension':beachTurrellResus}
    baseKernel=pset.Kernel(initialInput)+pset.Kernel(AdvectionRK4_floating)+pset.Kernel(BrownianMotion2D_floating)
    if scenario=='AdvectionDiffusionOnly':
        totalKernel=baseKernel+pset.Kernel(beachK[scenario])
    else:
        totalKernel=baseKernel+pset.Kernel(AntiBeachNudging)+beachK
    return totalKernel

