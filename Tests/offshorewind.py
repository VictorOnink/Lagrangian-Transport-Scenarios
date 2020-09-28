#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 12:35:09 2020

@author: victoronink

A place to mess around with the off-shore wind determination
"""
from parcels import ParticleSet,GeographicPolar,Geographic,FieldSet,Field,JITParticle,Variable,ErrorCode
import numpy as np
from netCDF4 import Dataset
import math
from datetime import datetime, timedelta
import os

def DeleteParticle(particle, fieldset, time):
    #This delete particle format from Philippe Delandmeter
    #https://github.com/OceanParcels/Parcelsv2.0PaperNorthSeaScripts/blob/master/northsea_mp_kernels.py
    print("Particle [%d] lost !! (%g %g %g %g)" % (particle.id, particle.lon, particle.lat, particle.depth, particle.time))
    particle.delete()

os.system('rm -r out-*')
#%%
datadirec='/Users/victoronink/Desktop/Bern Projects/Plastic Beaching/TestData/'

dataset=Dataset(datadirec+'ERA5-wind10m-06-2014.nc')
u10=dataset.variables['u10'][-1,:,:]
v10=dataset.variables['v10'][-1,:,:]
lonW=dataset.variables['longitude'][:]
lonW[lonW>180]-=360
latW=dataset.variables['latitude'][:]
fieldset = FieldSet(Field('U',u10,lon=lonW,lat=latW,mesh='spherical'),
                    Field('V',v10,lon=lonW,lat=latW,mesh='spherical'))
fieldset.U.units = GeographicPolar()
fieldset.V.units = Geographic()

uW=fieldset.U
vW=fieldset.V


datasetBor=Dataset(datadirec+'boundary_velocities_HYCOM.nc')
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
bU=fieldset.borU
bV=fieldset.borV

#%%
#Some data points for different cases
#case 1: borU and borV >0
point1=(58,4329)
latlon1=(latBor[point1[0]],lonBor[point1[1]])
#case 2: borU and borV <0
point2=(35,195)
latlon2=(latBor[point2[0]],lonBor[point2[1]])
pL={1:latlon1,2:latlon2}

# def dotproduct(v1, v2):
#   return sum((a*b) for a, b in zip(v1, v2))

# def length(v):
#   return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
    magV1=math.sqrt(v1[0]**2+v1[1]**2)
    magV2=math.sqrt(v2[0]**2+v2[1]**2)
    dot=v1[0]*v2[0]+v1[1]*v2[1]
    return math.acos(dot/(magV1*magV2))

for p in pL:
    loc=pL[p]
    print(loc)
    W=(uW.eval(1404086414999552.000000,0,loc[0],loc[1],applyConversion=False), #time,depth,lat,lon
        vW.eval(1404086414999552.000000,0,loc[0],loc[1],applyConversion=False)
        )
    # W2=(uW[1404086414999552,0,loc[0],loc[1]]*1852*60*math.cos(loc[0]*math.pi/180),
    #     vW[1404086414999552,0,loc[0],loc[1]]*1852*60)
    # print(W)
    # print(W2)
    magW=np.linalg.norm(W)
    B=(bU.eval(1404086414999552.000000,0,loc[0],loc[1],applyConversion=False)*-1, #time,depth,lat,lon
        bV.eval(1404086414999552.000000,0,loc[0],loc[1],applyConversion=False)*-1
        )
    B2=(bU[1404086414999552,0,loc[0],loc[1]]*1852*60*math.cos(loc[0]*math.pi/180),
        bV[1404086414999552,0,loc[0],loc[1]]*1852*60)
    # print(B)
    # print(B2)
    magB=np.linalg.norm(B)
    print('wind at point '+str(p)+' is '+str(W)+', magnitude = '+str(magW))
    print('border at point '+str(p)+' is '+str(B)+', magnitude = '+str(magB))
    alpha=angle(B,W)
    print('Angle between wind and border is '+str(alpha))
    print('Offshore wind component is '+str(magW*math.cos(alpha))+' m/s')

#%% So, the whole arrangment seems to work for the data, now we want to get this into a kernel
#First a fieldset with particles at point 1 and point 2.
class angleparticle(JITParticle):
    Angle = Variable('angle', 
                        dtype=np.float32,
                        initial=0)

lonlist=[lonBor[point1[1]],lonBor[point1[1]]]
latlist=[latBor[point1[0]],latBor[point1[0]]]
z=[0,10]
pset = ParticleSet(fieldset=fieldset, pclass=angleparticle, time=datetime(2014,6,30),
                    lon=lonlist,lat=latlist,depth=z)
print(pset)

def windAngle(particle,fieldset,time):
    t,d,la,lo=time, particle.depth, particle.lat, particle.lon
    #Load in the wind and border current components, and convert them both back to m/s
    bU,bV=fieldset.borU[t,d,la,lo]*-1*1852*60*math.cos(la*math.pi/180),fieldset.borV[t,d,la,lo]*-1*1852*60
    wU,wV=fieldset.U[t,d,la,lo]*1852*60*math.cos(la*math.pi/180),fieldset.V[t,d,la,lo]*1852*60
    #magnitude of the b and w vectors, and then dot product between them
    mB,mW=math.sqrt(bU**2+bV**2),math.sqrt(wU**2+wV**2)
    dot=bU*wU+bV*wV
    #Angle between b and w
    alpha=math.acos(dot/(mB*mW))
    #Get the offshore wind component
    wOff=mW*math.cos(alpha)
    particle.angle=wOff
    
    
windsample=pset.Kernel(windAngle)

pfile = pset.ParticleFile(name='test',
                          outputdt=timedelta(minutes=10))

pset.execute(windsample,
              runtime=timedelta(minutes=10),
              dt=timedelta(minutes=10),
              recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
              output_file=pfile
              )
print(pset)

pfile.export()
os.system('rm -r out-*')

