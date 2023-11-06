#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script performs shows the simplest way to use MyRocketSimulator:
Say 'Helllo World' to MRS!

The following MRS features are included:
    - Loading the MRS default mission (24 hours of ISS propagation)
    - Running the simulation
    - Use the MRSvilib to show the ground track 

"""

# import MRS libraries
from myrocketsimulator import MRSlib, MRSvislib

# make MRSMission object with the default MRS mission
MRSdemo0mission = MRSlib.MRSmission('defaultMRSmission')

# run mission
MRSdemo0mission.run_mission()

# add latitude and longitude data to the data frame 
MRSdemo0mission.expand_DFname(['EarthLLA'])

# make MRSvislib object with the performed mission
MRSdemo0viewer = MRSvislib.MRSviewer(MRSdemo0mission)

# show ground track
figGtroundTrack = MRSdemo0viewer.plot_GroundtrackEarth()