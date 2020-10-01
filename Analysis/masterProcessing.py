#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:53:54 2020

@author: victoronink
"""
import os
from FileNaming import FileNames
#%%
fn=FileNames(sce=int(os.environ['SCENARIO']),vic=int(os.environ['VICINITY']),
             st=int(os.environ['SHORETIME']),rt=int(os.environ['RESUSTIME']),
             sd=int(os.environ['SHOREDEPEN']),Wmin=int(os.environ['WMIN']),
             inp=int(os.environ['INPUT']),year=int(os.environ['STARTYEAR']),
             stoke=int(os.environ['STOKES']),ens=int(os.environ['ENSEMBLE']),
             simlen=int(os.environ['SIMLEN']))
os.system('echo The beaching scenario is '+fn.scenario)
fn.parcelsOutput()
date='01-01-2020'
suffix='_lon'
fn.timeSlice(date)