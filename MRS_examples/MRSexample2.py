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
This script performs a ten-hour propagation of the CREME spacraft and compares
it with a trajectory created by GMAT. Additionally, the orbital elements are 
compared one by one. 
Several perturbating forces are used in this mission.

The following MRS features are included:
    - load an external MRS mission file
    - adjust mission data 
        - change mission name
        - modify mission segments (remove delta-v maneuvers)
        - setup force models:
            - use spherical harmonics (n,m=35) for earth
            - use third body gravity (Sun, Moon)
            - use drag
            - use SRP
    - modify internal MRS variables
        - GM value for Earth (actually nt needed, as SH are used)
        - do not use space weather, but fixed default flux values
        - adjust atmospheric density (see Bautze-Scherff, 2023)
    - expand the data frame of the flown mission (here: orbital elements)
    - import external state vectors for comparison (source: GMAT generated csv)
    - use the MRSviewer
        - plot difference between simulation and external state vectors positions
    - direct access to data stored in the data frame (here: orbital elements)
        

The execution time is arund 2.3 minutes (i7 Quadcore 1.7 GHz). 
The plotted error should have a maxium value of ~0.35 meters.

The file ./GMATscripts/GMAT_example2_1hour.script is a GMAT script that can be
used to re-generate the GMAT csv report file. Please consider adjusting the path
of the report file.

"""

# load GMAT data for later comparison
GMATfilename = './referenceCSVs/GMAT_example2_10hours.csv'
GMATcsv = pd.read_csv(GMATfilename, sep='\s+', skiprows=0)  
GMATcsv.loc[:,'JD_TBD'] = GMATcsv.loc[:,'CREME.TDBModJulian'] + 2430000.0 # add JD-offset (ModJulian)

# load mission
CREMEmission = MRSlib.MRSmission('./MRSmissions/CREMEdefaultmission.py', checkmission=False)

# set mission specific name
CREMEmission.MD.name = 'MRS example 2 mission'

# drop delta-v maneuvers in mission segments
CREMEmission.MD.missionSegments = CREMEmission.MD.missionSegments.drop([1, 2]).reset_index(drop=True)

# adjust mission duratin to ten hours
CREMEmission.MD.missionSegments.loc[1, 'MET'] = 3600 * 10 # [s]

# adjust perturbating forces
CREMEmission.MD.forcesSettings.loc[0,'EarthSHn'] = 35 
CREMEmission.MD.forcesSettings.at[0,'planets'] = ['Earth','Sun', 'Moon'] # third bodies
CREMEmission.MD.forcesSettings.loc[0,'drag'] = 1 # use drag
CREMEmission.MD.forcesSettings.loc[0,'SRP'] = 1 # use SRP

# overwrite MRS default config for space weather 
CREMEmission.use_spaceweather = 0 # use fixed values for flux and geomag index (150, 150, 3)
CREMEmission.transToMSISE90 = 1 # modify density (using pre-set transformation) to match MSISE90

# check mission data
CREMEmission.check_MD()
 
# overwrite MRS default config for Earth-GM (default: DE421 values)
# needs to be done after check_MD()
CREMEmission.EARTH_GM = 398600441500000 # EGM96 instead of DE421 value for Earth GM

# run mission
CREMEmission.run_mission()

# add orbital elements (needed for later plot)
CREMEmission.expand_DFname(['EarthOrbElements'])

# load external state vectors; multiply by 1000 to convert from [km] to [m]
CREMEmission.load_comparisonStateVec(0, 1000*GMATcsv.iloc[:,1:7].values)

# generate viewer object
CREMEviewer = MRSvislib.MRSviewer(CREMEmission)

# plot GMAT to MRS position difference 
figEx2PosDiff, _ = CREMEviewer.plot_ComparisonPosDiff()

# save plot
figEx2PosDiff.savefig('./MRSoutput/MRSexample2_ComparisonPosDiff.svg', dpi=600)

# make a plot (using direct access mission data frame) of orbital element 
# differences between GMAT and MRS 

# make plot window
figEx2OEdiff, ax = plt.subplots(2,3, figsize=(9, 6))
plt.subplots_adjust(left=0.1,
                    bottom=0.11, 
                    right=0.945, 
                    top=0.938, 
                    wspace=0.44, 
                    hspace=0.374)


# adjust font sizes
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.labelsize': 12})
plt.rcParams.update({'axes.titlesize': 12})

# plot OE differences

ax[0,0].set_title('Semi-major axis')
ax[0,0].set_xlabel('MET [hours]')
ax[0,0].set_ylabel('SMA [m]')
ax[0,0].grid()
ax[0,0].plot(CREMEmission.missionDF.MET/3600,CREMEmission.missionDF.EarthOESMA - \
             GMATcsv.loc[:,'CREME.Earth.SMA']*1000)

ax[0,1].set_title('Eccentricity')
ax[0,1].set_xlabel('MET [hours]')
ax[0,1].set_ylabel('Eccenctricity [-]')
ax[0,1].grid()
ax[0,1].plot(CREMEmission.missionDF.MET/3600,CREMEmission.missionDF.EarthOEeccentricity - \
             GMATcsv.loc[:,'CREME.Earth.ECC'])

ax[0,2].set_title('Inclination')
ax[0,2].set_xlabel('MET [hours]')
ax[0,2].set_ylabel('Inclination [°]')
ax[0,2].grid()
ax[0,2].plot(CREMEmission.missionDF.MET/3600,CREMEmission.missionDF.EarthOEinclination - \
             GMATcsv.loc[:,'CREME.EarthMJ2000Eq.INC'])

ax[1,0].set_title('Arg. of Perigee')
ax[1,0].set_xlabel('MET [hours]')
ax[1,0].set_ylabel('Arg. of Perigee [°]')
ax[1,0].grid()
ax[1,0].plot(CREMEmission.missionDF.MET/3600,CREMEmission.missionDF.EarthOEargPeriapsis - \
             GMATcsv.loc[:,'CREME.EarthMJ2000Eq.AOP'])

ax[1,1].set_title('  Right Angle of Ascending Node')
ax[1,1].set_xlabel('MET [hours]')
ax[1,1].set_ylabel('RAAN [°]')
ax[1,1].grid()
ax[1,1].plot(CREMEmission.missionDF.MET/3600,CREMEmission.missionDF.EarthOERAAN - \
             GMATcsv.loc[:,'CREME.EarthMJ2000Eq.RAAN'])

ax[1,2].set_title('True Anomaly')
ax[1,2].set_xlabel('MET [hours]')
ax[1,2].set_ylabel('True Anomaly [°]')
ax[1,2].grid()
ax[1,2].plot(CREMEmission.missionDF.MET/3600,CREMEmission.missionDF.EarthOEtrueAnomaly - \
             GMATcsv.loc[:,'CREME.Earth.TA'])


figEx2OEdiff.savefig('./MRSoutput/MRSexample2_ComparisonOrbitalElements.svg', dpi=600)
