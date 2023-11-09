#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 21:54:02 2023

@author: thibault
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from importlib_resources import files
import simplekml

PKGNAME = 'myrocketsimulator'

class MRSviewer():
    """
    The MRSviewer class is initialized with the dataframe of an MRSmission. 
    It provides different methods to visualize the results of a flown mission.   
    """
    
    def __init__(self, MRSmissionobject):
        """
        Loads MRSmission object into the MRSviewer object.

        Parameters
        ----------
        MRSmissionobject : MRSmission object
            An MRSmission object that was execuded (i.e. run_mission()) and 
            contains a filled mission dataframe.

        Returns
        -------
        None.

        """
        
        # check if mission data frame exists
        if hasattr(MRSmissionobject, 'missionDF'):
            # load mission dataframe
            print('MRSviewer:\tLoading dataframe of mission '+MRSmissionobject.MD.name+'.')
            self.MD = MRSmissionobject.MD
            self.missionDF = MRSmissionobject.missionDF
            
        else:
            print('MRS:\t\tERROR: no mission data frame available.')
        
        return None
        
        
    def get_DFpointer(self, METrange):
        """
        

        Parameters
        ----------
        METrange : TYPE
            DESCRIPTION.

        Returns
        -------
        DFpointer : TYPE
            DESCRIPTION.

        """
        
        # if no valid METrange provided
        if len(METrange)!=2:
            DFpointer = np.array([0, len(self.missionDF)])
        # valid amount provided
        else:
            # find index  valuesclosest to start and end MET
            DFpointerStart = np.argmin(np.abs(self.missionDF['MET'].values - METrange[0]))
            DFpointerStop  = np.argmin(np.abs(self.missionDF['MET'].values - METrange[1]))
            DFpointer = np.array([DFpointerStart, DFpointerStop])
            
        return DFpointer
    
    
    def plot_GroundtrackEarth(self, METrange=[0]):
        """
        

        Parameters
        ----------
        METrange : TYPE, optional
            DESCRIPTION. The default is [0].

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        # get lat/lon values
        latvalues = self.missionDF['EarthLat'].values
        lonvalues = self.missionDF['EarthLon'].values
        
        # find jumps larger than 180 degrees
        coordelta = lonvalues[1:] - lonvalues[:-1]
        valdelta = (np.abs(coordelta)>180).nonzero()[0]
        # add last value as last point
        valdelta = np.append(valdelta,np.shape(lonvalues)[0]-1)
        
        # set up fig title
        figtitle = self.MD.name + ' Ground Track'
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(10, 5)
        fig.suptitle(figtitle, fontsize=16)
        ax.set_xlabel('longitude [°]')
        ax.set_ylabel('latitude [°]')
        # image source: https://de.wikipedia.org/wiki/Datei:Blue_Marble_2002.png
        backgroundimage = files(PKGNAME+'.data').joinpath('Blue_Marble_2002-2.png').as_posix()
        img = plt.imread(backgroundimage)
        ax.imshow(img, interpolation='bilinear', 
                        origin='upper', extent=[-180,180,-90,90])
        plt.xticks(np.arange(-180, 190, 30.0))
        plt.yticks(np.arange(-90, 100, 30.0))
        plt.subplots_adjust(left=.0, bottom=0.12, right=.99, top=0.90, wspace=0, hspace=0)
        ax.grid()
        
        # loop trough passes
        pointer = 0
        for i in np.nditer(valdelta):
            ax.plot(lonvalues[pointer:i],latvalues[pointer:i],'r', linewidth=1)
            pointer = i+1
            
        return fig
            
            
    def plot_GroundtrackMoon(self, METrange=[0]):
        """
        

        Parameters
        ----------
        METrange : TYPE, optional
            DESCRIPTION. The default is [0].

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        # get lat/lon values
        latvalues = self.missionDF['MoonLat'].values
        lonvalues = self.missionDF['MoonLon'].values
        
        # find jumps larger than 180 degrees
        coordelta = lonvalues[1:] - lonvalues[:-1]
        valdelta = (np.abs(coordelta)>180).nonzero()[0]
        # add last value as last point
        valdelta = np.append(valdelta,np.shape(lonvalues)[0]-1)
        
        # set up fig title
        figtitle = self.MD.name + ' Ground Track'
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(10, 5)
        fig.suptitle(figtitle, fontsize=16)
        ax.set_xlabel('longitude [°]')
        ax.set_ylabel('latitude [°]')
        # image source: https://svs.gsfc.nasa.gov/cgi-bin/details.cgi?aid=4720
        backgroundimage = files(PKGNAME+'.data').joinpath('MoonSurface.jpg').as_posix()
        img = plt.imread(backgroundimage)
        ax.imshow(img, interpolation='bilinear', 
                        origin='upper', extent=[-180,180,-90,90])
        plt.xticks(np.arange(-180, 190, 30.0))
        plt.yticks(np.arange(-90, 100, 30.0))
        plt.subplots_adjust(left=.0, bottom=0.12, right=.99, top=0.90, wspace=0, hspace=0)
        ax.grid()
        
        # loop trough passes
        pointer = 0
        for i in np.nditer(valdelta):
            ax.plot(lonvalues[pointer:i],latvalues[pointer:i],'r', linewidth=1)
            pointer = i+1
        
        return fig
        
    def plot_EarthOE(self, METrange=[0]):
        """
        

        Parameters
        ----------
        METrange : TYPE, optional
            DESCRIPTION. The default is [0].

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        figtitle = self.MD.name + ' Orbital Elements'
        fig, ax = plt.subplots(2,3, figsize=(9, 6))
        fig.suptitle(figtitle, fontsize=16)
       
        
        plt.subplots_adjust(left=0.1,
                            bottom=0.11, 
                            right=0.945, 
                            top=0.89, 
                            wspace=0.44, 
                            hspace=0.374)
        
        # adjust font sizes
        plt.rcParams.update({'font.size': 12})
        plt.rcParams.update({'axes.labelsize': 12})
        plt.rcParams.update({'axes.titlesize': 12})

        
        # SMA
        ax[0,0].set_title('Semi-major axis')
        ax[0,0].set_xlabel('MET [hours]')
        ax[0,0].set_ylabel('SMA [km]')
        ax[0,0].grid()
        ax[0,0].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.EarthOESMA[DFpointer[0]:DFpointer[1]]/1000,
                linewidth=1.5)
        
        # Eccentricity
        ax[0,1].set_title('Eccentricity')
        ax[0,1].set_xlabel('MET [hours]')
        ax[0,1].set_ylabel('Eccenctricity [-]')
        ax[0,1].grid()
        ax[0,1].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.EarthOEeccentricity[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        # Inclination
        ax[0,2].set_title('Inclination')
        ax[0,2].set_xlabel('MET [hours]')
        ax[0,2].set_ylabel('Inclination. [°]')
        ax[0,2].grid()
        ax[0,2].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.EarthOEinclination[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        
        # Argument of Perigee/Periapsis
        ax[1,0].set_title('Arg. of Perigee')
        ax[1,0].set_xlabel('MET [hours]')
        ax[1,0].set_ylabel('Arg. of Perigee [°]')
        ax[1,0].grid()
        ax[1,0].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.EarthOEargPeriapsis[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        # RAAN
        ax[1,1].set_title('Right Angle of Ascending Node')
        ax[1,1].set_xlabel('MET [hours]')
        ax[1,1].set_ylabel('RAAN [°]')
        ax[1,1].grid()
        ax[1,1].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.EarthOERAAN[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        # True Anomaly
        ax[1,2].set_title('True Anomaly')
        ax[1,2].set_xlabel('MET [hours]')
        ax[1,2].set_ylabel('True Anomaly. [°]')
        ax[1,2].grid()
        ax[1,2].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.EarthOEtrueAnomaly[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        return fig
        
    def plot_MoonOE(self, METrange=[0]):
        """
        

        Parameters
        ----------
        METrange : TYPE, optional
            DESCRIPTION. The default is [0].

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        figtitle = self.MD.name + ' Orbital Elements (Moon)'
        fig, ax = plt.subplots(2,3, figsize=(9, 6))
        fig.suptitle(figtitle, fontsize=16)
       
        
        plt.subplots_adjust(left=0.1,
                            bottom=0.11, 
                            right=0.945, 
                            top=0.89, 
                            wspace=0.44, 
                            hspace=0.374)
        
        # adjust font sizes
        plt.rcParams.update({'font.size': 12})
        plt.rcParams.update({'axes.labelsize': 12})
        plt.rcParams.update({'axes.titlesize': 12})
        
        # SMA
        ax[0,0].set_title('Semi-major axis')
        ax[0,0].set_xlabel('MET [hours]')
        ax[0,0].set_ylabel('SMA [km]')
        ax[0,0].grid()
        ax[0,0].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.MoonOESMA[DFpointer[0]:DFpointer[1]]/1000,
                linewidth=1.5)
        
        # Eccentricity
        ax[0,1].set_title('Eccentricity')
        ax[0,1].set_xlabel('MET [hours]')
        ax[0,1].set_ylabel('Eccenctricity [-]')
        ax[0,1].grid()
        ax[0,1].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.MoonOEeccentricity[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        # Inclination
        ax[0,2].set_title('Inclination')
        ax[0,2].set_xlabel('MET [hours]')
        ax[0,2].set_ylabel('Inclination. [°]')
        ax[0,2].grid()
        ax[0,2].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.MoonOEinclination[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        
        # Argument of Perigee/Periapsis
        ax[1,0].set_title('Arg. of Perigee')
        ax[1,0].set_xlabel('MET [hours]')
        ax[1,0].set_ylabel('Arg. of Perigee [°]')
        ax[1,0].grid()
        ax[1,0].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.MoonOEargPeriapsis[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        # RAAN
        ax[1,1].set_title('Right Angle of Ascending Node')
        ax[1,1].set_xlabel('MET [hours]')
        ax[1,1].set_ylabel('RAAN [°]')
        ax[1,1].grid()
        ax[1,1].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.MoonOERAAN[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        # True Anomaly
        ax[1,2].set_title('True Anomaly')
        ax[1,2].set_xlabel('MET [hours]')
        ax[1,2].set_ylabel('True Anomaly. [°]')
        ax[1,2].grid()
        ax[1,2].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.MoonOEtrueAnomaly[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        return fig
        
    def plot_ComparisonPosDiff(self, METrange=[0]):
        """
        

        Parameters
        ----------
        METrange : TYPE, optional
            DESCRIPTION. The default is [0].

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        figtitle = self.MD.name + ' Position Delta to Reference'
        fig, ax = plt.subplots(1,1)
        #fig.set_size_inches(10, 5)
        ax.set_title(figtitle)
        ax.set_xlabel('MET [hours]')
        ax.set_ylabel('Error [m]')
        #plt.subplots_adjust(left=.0, bottom=0.12, right=.99, top=0.90, wspace=0, hspace=0)
        ax.grid()
        ax.plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.compPosDiff[DFpointer[0]:DFpointer[1]],
                linewidth=2)
        
        return fig
        
    def plot_ComparisonStateVecDiff(self, METrange=[0]):
        """
        

        Parameters
        ----------
        METrange : TYPE, optional
            DESCRIPTION. The default is [0].

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        figtitle = self.MD.name + ' State Vector Differences'
        fig, ax = plt.subplots(2,3)
        fig.suptitle(figtitle, fontsize=16)
        #fig.set_size_inches(10, 5)
        
        # x
        ax[0,0].set_title('dx')
        ax[0,0].set_xlabel('MET [hours]')
        ax[0,0].set_ylabel('dx [m]')
        ax[0,0].grid()
        ax[0,0].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.x[DFpointer[0]:DFpointer[1]]-self.missionDF.compX[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        # y
        ax[0,1].set_title('dy')
        ax[0,1].set_xlabel('MET [hours]')
        ax[0,1].set_ylabel('dy [m]')
        ax[0,1].grid()
        ax[0,1].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.y[DFpointer[0]:DFpointer[1]]-self.missionDF.compY[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        # z
        ax[0,2].set_title('dz')
        ax[0,2].set_xlabel('MET [hours]')
        ax[0,2].set_ylabel('dz [m]')
        ax[0,2].grid()
        ax[0,2].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.z[DFpointer[0]:DFpointer[1]]-self.missionDF.compZ[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        
        # vx
        ax[1,0].set_title('dvx')
        ax[1,0].set_xlabel('MET [hours]')
        ax[1,0].set_ylabel('dvx [m/s]')
        ax[1,0].grid()
        ax[1,0].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.vx[DFpointer[0]:DFpointer[1]]-self.missionDF.compVx[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        # vy
        ax[1,1].set_title('dvy')
        ax[1,1].set_xlabel('MET [hours]')
        ax[1,1].set_ylabel('dvy [m/s]')
        ax[1,1].grid()
        ax[1,1].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.vy[DFpointer[0]:DFpointer[1]]-self.missionDF.compVy[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        # vz
        ax[1,2].set_title('dvz')
        ax[1,2].set_xlabel('MET [hours]')
        ax[1,2].set_ylabel('dvz [m/s]')
        ax[1,2].grid()
        ax[1,2].plot(self.missionDF.MET[DFpointer[0]:DFpointer[1]]/3600,
                self.missionDF.vz[DFpointer[0]:DFpointer[1]]-self.missionDF.compVz[DFpointer[0]:DFpointer[1]],
                linewidth=1.5)
        
        return fig
        
        
    def plot_MoonRPS(self, METrange=[0]):
        """
        

        Parameters
        ----------
        METrange : TYPE, optional
            DESCRIPTION. The default is [0].

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        # extract RPS positions
        RPSposition = self.missionDF.loc[DFpointer[0]:DFpointer[1], ['MoonRPSx', 'MoonRPSy', 'MoonRPSz']].values
        
        # set up plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('RPS X')
        ax.set_ylabel('RPS Y')
        ax.set_zlabel('RPS Z')
        figtitle =  self.MD.name + ' trajectory in Rotating-Pulsating System (RPS)'
        ax.set_title(figtitle)
        # plot data 
        ax.plot(RPSposition[:,0], RPSposition[:,1], RPSposition[:,2], label='XYZ trajectory in RPS')
        # fine tune plot
        ax.set_xlim([-0.0, 1.0])
        ax.set_ylim([-0.5, 0.5])
        ax.set_zlim([-0.5, 0.5])
        ax.legend()
        ax.view_init(elev=35., azim=-35)
        
        return fig
        
    def plot_GCRF_XYpos(self, METrange=[0], withMoon=0):
        """
        

        Parameters
        ----------
        METrange : TYPE, optional
            DESCRIPTION. The default is [0].
        withMoon : TYPE, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        figtitle = self.MD.name + ' GCRF X/Y Position'
        fig, ax = plt.subplots(1,1)
        #fig.set_size_inches(10, 5)
        ax.set_title(figtitle)
        ax.set_xlabel('ECI X [1000 km]')
        ax.set_ylabel('ECI Y [1000 km]')
        #plt.subplots_adjust(left=.0, bottom=0.12, right=.99, top=0.90, wspace=0, hspace=0)
        ax.grid()
       
        ax.axis('equal')
        ax.plot(self.missionDF.x[DFpointer[0]:DFpointer[1]]/1000**2,
                self.missionDF.y[DFpointer[0]:DFpointer[1]]/1000**2,
                label='Spacecraft XY (ECI) trajectory', linewidth=2)
        
        if withMoon and 'MoonPosX' in self.missionDF.columns:
            ax.plot(self.missionDF.MoonPosX[DFpointer[0]:DFpointer[1]]/1000**2,
                    self.missionDF.MoonPosY[DFpointer[0]:DFpointer[1]]/1000**2,
                    label='Moon XY (ECI) trajectory', linewidth=2)
            
        # add legend
        ax.legend()
        
        return fig
    
    def plot_GCRF_orbit(self, METrange=[0]):
        """
        

        Parameters
        ----------
        METrange : TYPE, optional
            DESCRIPTION. The default is [0].

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        # find segments jumps
        segvals = self.missionDF.segmentID[DFpointer[0]:DFpointer[1]].values
        dsegvals = segvals[1:] - segvals[:-1]
        segjumps = (dsegvals>0).nonzero()[0] + DFpointer[0] + 1 # +1 to include last value of segment
        
        # add end value for last segment 
        segjumps = np.append(segjumps, DFpointer[1])
            
        # plot orbit/trajectory
        figtitle = self.MD.name + ' GCRF 3D Orbit'
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title(figtitle)
        ax.set_xlabel('X [km]',fontsize=12)
        ax.set_ylabel('Y [km]',fontsize=12)
        ax.set_zlabel('Z [km]',fontsize=12)
        
        # adjust font sizes
        plt.rcParams.update({'font.size': 12})
        plt.rcParams.update({'axes.labelsize': 12})
        plt.rcParams.update({'axes.titlesize': 12})
        
        # loop through segments
        pointer = DFpointer[0]
        for i in range(len(segjumps)):
            ax.plot(self.missionDF.x[pointer:segjumps[i]]/1000, 
                    self.missionDF.y[pointer:segjumps[i]]/1000,
                    self.missionDF.z[pointer:segjumps[i]]/1000, 
                    label='Segment {}'.format(i), linestyle='-')
            pointer = segjumps[i]

        # show legend        
        ax.legend()
        
        # plot Earth 
        EarthR = 6378.137 
        u, v = np.mgrid[0:2 * np.pi:50j, 0:np.pi:50j]
        x = EarthR * np.cos(u) * np.sin(v)
        y = EarthR * np.sin(u) * np.sin(v)
        z = EarthR * np.cos(v)
        ax.plot_surface(x, y, z, color='blue', alpha=0.3)     
        
        return fig
        
        
    def export_Earth_KML(self, METrange=[0], folder='./'):
        """
        Export KML file of trajetory around Earth.

        Parameters
        ----------
        METrange : TYPE, optional
            DESCRIPTION. The default is [0].

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        # make path to file
        filename = folder +' '+ self.MD.name +' Earth_Orbit.kml'
        
        # get coords
        mycoords = self.missionDF.loc[DFpointer[0]:DFpointer[1],
                                      ['EarthLat', 'EarthLon', 'EarthAlt']].to_numpy()

        # invert lat/lon
        mycoords.T[[0,1]] = mycoords.T[[1,0]]
        
        # make tuple
        tuplecoords = tuple(map(tuple, mycoords))
        
        # Create an instance of Kml
        kml = simplekml.Kml(open=1)
        
        linestring = kml.newlinestring(name=self.MD.name)
        
        for i in range(len(tuplecoords)):
            linestring.coords.addcoordinates([tuplecoords[i]])
        
        linestring.altitudemode = simplekml.AltitudeMode.absolute
        linestring.style.linestyle.color = 'ff0000ff'
        linestring.style.linestyle.width = 2
        
        # Save the KML
        kml.save(filename)
        
        return None
        

    def export_Moon_KML(self, METrange=[0], folder='./'):
        """
        Export KML file of trajetory around Moon.

        Parameters
        ----------
        METrange : TYPE, optional
            DESCRIPTION. The default is [0].

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        # make path to file
        filename = folder +' '+ self.MD.name +' Moon_Orbit.kml'
        
        # get coords
        mycoords = self.missionDF.loc[DFpointer[0]:DFpointer[1],
                                      ['MoonLat', 'MoonLon', 'MoonAlt']].to_numpy()

        # invert lat/lon
        mycoords.T[[0,1]] = mycoords.T[[1,0]]
        
        # adjust altitude
        mycoords[:,2] -= 1737.4e3 
        
        # make tuple
        tuplecoords = tuple(map(tuple, mycoords))
        
        # Create an instance of Kml
        kml = simplekml.Kml(open=1)
        
        linestring = kml.newlinestring(name=self.MD.name)
        
        for i in range(len(tuplecoords)):
            linestring.coords.addcoordinates([tuplecoords[i]])
        
        linestring.altitudemode = simplekml.AltitudeMode.absolute
        linestring.style.linestyle.color = 'ff0000ff'
        linestring.style.linestyle.width = 2
        
        # Save the KML
        kml.save(filename)
        
        return None
        
        
        