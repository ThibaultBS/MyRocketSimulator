#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 10:14:09 2023

@author: thibault
"""


from myrocketsimulator import MRSlib, MRSvislib
import pandas as pd
import matplotlib.pyplot as plt


"""
This script performs a five-hour propagation of the CREME spacraft, including
two delta-v maneuvers to raise the spacecraft to a much higher orbit 
(Hohmann transfer).

The values for the delta-v maneuvers are calculated on:
https://www.omnicalculator.com/physics/hohmann-transfer
with the following parameters:
    - Gravitational parameter: 398600.4415 km^3/s^2
    - Radius primary body: 6378.137 km
    - init. orbit radius: 6978.130 km
    - dest. orbit radius: 8500.000 km

Resulting maneuvers
    - Delta-v #1
        - added velocity: 362.84944 m/s
    - Delta-v #2
        - time after delta-v #1: 3387.7685 s
        - added velocity: 345.3667 m/s
        
    
The following MRS features are included:
    - load an external MRS mission file
    - adjust mission data 
        - change mission name
        - modify delta-v maneuvers
            - change time
            - change velocity values
    - expand the data frame of the flown mission (here: orbital elements)
    - use the MRSviewer
        - plot orbital elements
        - 3D plot of trajectory in GCRF
        - manual saving of figures
    - export mission data frame as csv-file
        

The execution time is arund 2 seconds (i7 Quadcore 1.7 GHz). 


"""

# load mission
CREMEmission = MRSlib.MRSmission('./MRSmissions/CREMEdefaultmission.py', checkmission=False)

# set mission specific name
CREMEmission.MD.name = 'MRS example 3 mission'

# adjust delta-v maneuvers times
CREMEmission.MD.missionSegments.loc[1, 'MET'] = 5500 # 1st maneuver
CREMEmission.MD.missionSegments.loc[2, 'MET'] = 5500 + 3387.7685 #2nd maneuver

# adjust delta-v velocity values (in flight direction / VNB frame)
CREMEmission.MD.maneuverSettings.loc[0,'dx'] = 362.84944
CREMEmission.MD.maneuverSettings.loc[1,'dx'] = 345.3667

# adjust mission duratin to five hours
CREMEmission.MD.missionSegments.loc[3, 'MET'] = 3600 * 5 # [s]

# check mission data
CREMEmission.check_MD()
 
# overwrite MRS default config for Earth-GM (default: DE421 values)
# needs to be done after check_MD()
CREMEmission.EARTH_GM = 398600441500000 # EGM96 instead of DE421 value for Earth GM

# run mission
CREMEmission.run_mission()

# add orbital elements (needed for later plot)
CREMEmission.expand_DFname(['EarthOrbElements'])

# export mission data frame
CREMEmission.exportDataframes(folder='./MRSoutput/',missionDFonly=True)

# generate viewer object
CREMEviewer = MRSvislib.MRSviewer(CREMEmission)

# 3D view
figEx3GCRForbit, _ = CREMEviewer.plot_GCRF_orbit()
figEx3GCRForbit.axes[0].view_init(16,25) # rotate view

# save plot
figEx3GCRForbit.savefig('./MRSoutput/MRSexample3_GCRForbit.svg', dpi=600)

# show orbital elements
figEx3OE, _ = CREMEviewer.plot_EarthOE()

# save plot
figEx3OE.savefig('./MRSoutput/MRSexample3_OE.svg', dpi=600)

