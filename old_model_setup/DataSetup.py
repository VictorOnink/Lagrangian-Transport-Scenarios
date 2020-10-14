#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 18:25:35 2020

@author: victoronink
"""

import os
from netCDF4 import Dataset
import numpy as np
from parcels import GeographicPolar,Geographic,FieldSet,Field
from ParticleBeaching import ProbShore
from datetime import timedelta
import math
import glob

scenarionames={0:'AdvectionDiffusionOnly',1:'CoastalProximity',2:'Stochastic',
               3:'ShoreDependentResuspension',4:'TurrellResuspension'}
scenario=scenarionames[int(os.environ['SCENARIO'])]
vicinity=int(os.environ['VICINITY']) #days
shoreTime,resusTime=int(os.environ['SHORETIME']),int(os.environ['RESUSTIME']) #days, days
shoreDepen = int(os.environ['SHOREDEPEN'])
Wmin=int(os.environ['WMIN'])

def CreateFieldSet(server,stokes,scenario):
    """
    Creating the fieldset object containing all the fields

    Parameters
    ----------
    server : int
        The server I'm running my simulations on. 0 = alphadata04
                                                  1 = ubelix
    stokes : int
        Do we include stokes drift. The default is 0. 0 = Include Stokes
                                                      1 = Do not include Stokes
    scenario : int
        which beaching scenario. 0 = Advection & Diffusion Only
                                 1 = Coastal vicinity
                                 2 = Stochastic
                                 3 = Shore dependent resuspension
                                 4 = Turrell (2020) resuspension

    Returns
    -------
    fieldset.

    """
    #Defining the folders in which all the data is stored on the different servers
    if server==0:
        datadirec="/alphadata04/onink/lagrangian_sim/"
        dataInputdirec="/alphadata04/onink/lagrangian_sim/BeachingSim/Input/"
    elif server==1:
        datadirec="/home/ubelix/climate/shared/onink/"
        dataInputdirec="/home/ubelix/climate/shared/onink/Input/"
        
    os.system('echo "Creating the main fieldset"')
    #Loading in the surface currents and Stokes drift
    filenames = {'U': datadirec+"HYCOM/HYCOM_Surface*2010*.nc",
                 'V': datadirec+"HYCOM/HYCOM_Surface*2010*.nc",
                 'Ust': datadirec+"WaveWatchIIIstokes/ww3.2010*_uss.nc",
                 'Vst': datadirec+"WaveWatchIIIstokes/ww3.2010*_uss.nc",
                    }
    variables = {'U': 'water_u',
                 'V': 'water_v',
                 'Ust': 'uuss',
                 'Vst': 'vuss',
                    }
    dimensions = {'U':{'time': 'time','depth':'depth','lat': 'lat','lon': 'lon'},
                  'V':{'time': 'time','depth':'depth','lat': 'lat','lon': 'lon'},
                  'Ust':{'lat': 'latitude','lon': 'longitude','time': 'time'},
                  'Vst':{'lat': 'latitude','lon': 'longitude','time': 'time'},
                  }
    
    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions,allow_time_extrapolation=True)
    
    ###########################################################################
    # Adding the Stokes drift to the HYCOM currents                           #
    ###########################################################################
    os.system('echo "Adding Stokes drift"')
    # if stokes==0:
    #     fieldset.Ust.units = GeographicPolar()
    #     fieldset.Vst.units = Geographic()
    #     fieldset=FieldSet(U=fieldset.U+fieldset.Ust,
    #                       V=fieldset.V+fieldset.Vst,
    #                       )
    # else:
    #     fieldset=FieldSet(U=fieldset.U,
    #                       V=fieldset.V
    #                       )
    fieldset.Ust.units = GeographicPolar()
    fieldset.Vst.units = Geographic()
    ###########################################################################
    # Adding the border current, which applies for all scenarios except for 0 #
    ###########################################################################
    if scenario!='AdvectionDiffusionOnly':
        os.system('echo "Adding the border current"')
        datasetBor=Dataset(dataInputdirec+'boundary_velocities_HYCOM.nc')
        borU=datasetBor.variables['MaskUvel'][0,:,:]
        borV=datasetBor.variables['MaskVvel'][0,:,:]
        #Normalizing the border current so that the total current is always 1m/s
        borMag=np.sqrt(np.square(borU)+np.square(borV))
        borMag[borMag==0]=1
        borU=np.divide(borU,borMag)
        borV=np.divide(borV,borMag)
        #Adding the actual field
        lonBor,latBor=datasetBor.variables['lon'][:],datasetBor.variables['lat'][:]
        fieldset.add_field(Field('borU', borU,lon=lonBor,lat=latBor,mesh='spherical'))
        fieldset.add_field(Field('borV', borV,lon=lonBor,lat=latBor,mesh='spherical'))
        #making sure the units are interpreted as 1 m/s
        fieldset.borU.units = GeographicPolar()
        fieldset.borV.units = Geographic()
    ###########################################################################
    # Adding the horizontal diffusion                                         #
    ###########################################################################
    os.system('echo "Adding diffusion"')
    kh=10 #m^2 s^-1, following Lacerda et al. (2019) and Liubertseva et al. (2018)
    dataset=Dataset(datadirec+'HYCOM/HYCOM_Surface_3h_2000-01-01.nc')
    uo=dataset.variables['water_u'][0,0,:,:]
    lat_kh=dataset.variables['lat'][:]
    lon_kh=dataset.variables['lon'][:]
    kh_f=kh*np.ones(uo.shape)
    kh_f[uo.mask==True]=0
    fieldset.add_field(Field('Kh_zonal', kh_f,lon=lon_kh,lat=lat_kh,mesh='spherical'))
    fieldset.add_field(Field('Kh_meridional', kh_f,lon=lon_kh,lat=lat_kh,mesh='spherical'))
    ###########################################################################
    # Adding in the  land cell identifiers                                    #
    ###########################################################################
    os.system('echo "Adding land/water boolean field"')
    landID=np.load(dataInputdirec+'land_cell_identifier.npy')
    fieldset.add_field(Field('landID', landID,lon=lonBor,lat=latBor,mesh='spherical'))
    ###########################################################################
    # Distance to the shore                                                   #
    ###########################################################################
    os.system('echo "Adding distance to shore"')
    datasetCoast=Dataset(dataInputdirec+'distance2coast.nc')
    distance=datasetCoast.variables['distance'][0,:,:]
    lonD,latD=datasetCoast.variables['lon'][:],datasetCoast.variables['lat'][:]
    fieldset.add_field(Field('distance2shore', distance,lon=lonD,lat=latD,mesh='spherical'))
    
    ###########################################################################
    # Adding fieldset constants                                               #
    ###########################################################################
    if scenario=='CoastalProximity':
        #The vicinity timescale
        fieldset.add_constant('vic',vicinity)
    if scenario=='Stochastic':
        # Beaching and resuspension particles are global constants, so now they
        # don't need to be recomputed every timestep
        p_b=math.exp(-timedelta(minutes=10).total_seconds()/(shoreTime*86400.))
        p_r=math.exp(-timedelta(minutes=10).total_seconds()/(resusTime*86400.))
        fieldset.add_constant('p_beach',p_b)
        fieldset.add_constant('p_resus',p_r)
    if scenario=='ShoreDependentResuspension':
        #Here only the beaching probability is a global constant, the resuspension
        #probability will instead be represented using a field
        p_b=math.exp(-timedelta(minutes=10).total_seconds()/(shoreTime*86400.))
        fieldset.add_constant('p_beach',p_b)
    if scenario=='TurrellResuspension':
        #The global constant resuspension probability
        p_b=math.exp(-timedelta(minutes=10).total_seconds()/(shoreTime*86400.))
        fieldset.add_constant('p_beach',p_b)
        #The minimum wind speed for 
        fieldset.add_constant('Wmin',Wmin/10)
        
    ###########################################################################
    # The coastline type                                                      #
    ###########################################################################
    if scenario=='ShoreDependentResuspension':
        os.system('echo "Adding coastline type"')
        s=np.load(dataInputdirec+'coastline_sand_vs_not_sand.npy')
        p_r=np.exp(-timedelta(minutes=10).total_seconds()/(ProbShore(shoreDepen,scenario,s)*86400.))        
        fieldset.add_field(Field('p_resus', p_r,lon=lon_kh,lat=lat_kh,mesh='spherical'))
    ###########################################################################
    # The 10m winds                                                           #
    ###########################################################################
    if scenario=='TurrellResuspension':
        os.system('echo "Adding 10m winds"')
        windfiles = glob.glob(datadirec+"Wind/ERA5-wind10m*.nc")
        windfiles.sort()
        filenames = {'u10': windfiles,
                      'v10': windfiles}
        variables = {'u10': 'u10','v10': 'v10'}
        dimensions = {'time': 'time','lat': 'latitude','lon': 'longitude'}
        #Creating a fieldset for the wind data
        fieldset_wind = FieldSet.from_netcdf(filenames,variables,dimensions,allow_time_extrapolation=True)
        fieldset_wind.u10.units = GeographicPolar()
        fieldset_wind.v10.units = Geographic()
        #Adding the wind fields to the general fieldset
        fieldset.add_field(fieldset_wind.u10)
        fieldset.add_field(fieldset_wind.v10)
    ###########################################################################
    # Sea surface elevation                                                   #
    ###########################################################################
    if scenario=='TurrellResuspension':
        os.system('echo "Adding sea surface elevation"')
        elevfiles = glob.glob(datadirec+"HYCOM/HYCOM_SeaEleve_3h_20*.nc")
        elevfiles.sort()
        filenames = {'eta': elevfiles}
        variables = {'eta': 'surf_el'}
        dimensions = {'time': 'time','lat': 'lat','lon': 'lon'}
        #Creating a fieldset for the wind data
        fieldset_sea = FieldSet.from_netcdf(filenames,variables,dimensions,allow_time_extrapolation=True)
        #Adding the wind fields to the general fieldset
        fieldset.add_field(fieldset_sea.eta)
    ###########################################################################
    # Now the periodic halo for when we go across the 180/-180 degree line    #
    ###########################################################################
    os.system('echo "Finally, the periodic halo"')
    fieldset.add_periodic_halo(zonal=True)

    return fieldset
