#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 10:14:09 2023

@author: thibault
"""


from myrocketsimulator import MRSlib, MRSvislib
import pandas as pd


"""
This script performs a one-day propagation of the CREME spacraft and compares
it with a trajectory created by GMAT. Only Earth's gravity (point mass) is 
used to accelerate the spacecraft, i.e., no perturbating forces are activated.

The following MRS features are included:
    - load an external MRS mission file
    - adjust mission data 
        - change mission name
        - modify mission segments (remove delta-v maneuvers)
    - modify internal MRS variables (GM value for Earth)
    - expand the data frame of the flown mission (here: longitude, latitude, altitude)
    - import external state vectors for comparison (source: GMAT generated csv)
    - use the MRSviewer
        - plot difference between simulation and external state vectors positions
        - plot ground track (Earth)

The execution time is arund 8 seconds (i7 Quadcore 1.7 GHz). 
The plotted error should have a maxium value of ~4e-5 meters (40 micrometers).

The file ./GMATscripts/GMAT_example1_24hours.script is a GMAT script that can be
used to re-generate the GMAT csv report file. Please consider adjusting the path
of the report file.

"""

# load GMAT data for later comparison
GMATfilename = './referenceCSVs/GMAT_example1_24hours.csv'
GMATcsv = pd.read_csv(GMATfilename, sep='\s+', skiprows=0)  
GMATcsv.loc[:,'JD_TBD'] = GMATcsv.loc[:,'CREME.TDBModJulian'] + 2430000.0 # add JD-offset (ModJulian)

# load mission
CREMEmission = MRSlib.MRSmission('./MRSmissions/CREMEdefaultmission.py', checkmission=False)

# set mission specific name
CREMEmission.MD.name = 'MRS example 1 mission'

# drop delta-v maneuvers in mission segments
CREMEmission.MD.missionSegments = CREMEmission.MD.missionSegments.drop([1, 2]).reset_index(drop=True)
  
# check mission data
CREMEmission.check_MD()
 
# overwrite MRS default config for Earth-GM (default: DE421 values)
# needs to be done after check_MD()
CREMEmission.EARTH_GM = 398600441500000 # EGM96 instead of DE421 value for Earth GM

# run mission
CREMEmission.run_mission()

# add LLA values (needed for ground track)
CREMEmission.expand_DFname(['EarthLLA'])

# load external state vectors; multiply by 1000 to convert from [km] to [m]
CREMEmission.load_comparisonStateVec(0, 1000*GMATcsv.iloc[:,1:7].values)

# generate viewer object
CREMEviewer = MRSvislib.MRSviewer(CREMEmission)

# plot GMAT to MRS position difference 
figEx1PosDiff=CREMEviewer.plot_ComparisonPosDiff()

# save plot
figEx1PosDiff.savefig('./MRSoutput/MRSexample1_ComparisonPosDiff.svg', dpi=600)

# plot ground track
figEx1GroundtrackEarth = CREMEviewer.plot_GroundtrackEarth()

# save plot
figEx1GroundtrackEarth.savefig('./MRSoutput/MRSexample1_GroundtrackEarth.svg', dpi=600)