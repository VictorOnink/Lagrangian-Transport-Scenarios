#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:53:54 2020

@author: victoronink
"""
import os
from FileNaming import parameters
#%%
os.system('echo '+str(globals()))
os.system('echo '+str(int(os.environ['SCENARIO'])))
os.system('echo The beaching scenario is '+parameters.scenario)
parameters.printScenario()
