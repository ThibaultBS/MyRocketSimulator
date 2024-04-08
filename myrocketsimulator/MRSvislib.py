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
from prettytable import PrettyTable
import simplekml
import re

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
            self.SC = MRSmissionobject.SC
            self.missionDF = MRSmissionobject.missionDF
            
            
        else:
            print('MRS:\t\tERROR: no mission data frame available.')
            
        
        if hasattr(MRSmissionobject, 'eventsDF'):
            self.eventsDF = MRSmissionobject.eventsDF
        
        
        return None
        
        
    def get_DFpointer(self, METrange):
        """
        Gets start and end pointers of dataset for given MET range. Internal function.

        Parameters
        ----------
        METrange : Array of two floats
            Start and End MET values.

        Returns
        -------
        DFpointer : array of two ints
            Start and end rows of mission DF.

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
    
    
    def plotDatasets(self,plot_title, plot_xlabel, plot_ylabel1, plot_ylabel2,
                     plot_xdata, plot_ydata1, plot_ydata2, plot_legend1, plot_legend2,
                     bgimage, bgsize, ylim2):
        """
        Plots one or two datasets from the mission data frame. Internal function.
        Users should call the missionDFplotTS() metho.d

        Parameters
        ----------
        plot_title : string
            Plot title.
        plot_xlabel : string
            Label of x axis.
        plot_ylabel1 : string
            Label of left y axis.
        plot_ylabel2 : string
            Label of right y axis.
        plot_xdata : array of floats
            X data.
        plot_ydata1 : array of floats
            Y data for left axis.
        plot_ydata2 : array of floats
            Y data for right axis.
        plot_legend1 : string
            Legend for left axis data.
        plot_legend2 : string
            Legend for right axis data.
        bgimage : string
            Filename/path to background image.
        bgsize : array of floats
            X/Y limits of background image.
        ylim2 : array of two floats
            Limitation of right y axis.

        Returns
        -------
        fig : figure handle
            No description.
        ax : axis handle
            No description.

        """
        
        fig, ax = plt.subplots(1,1)
        ax.set_title(plot_title)
        ax.set_xlabel(plot_xlabel)
        ax.grid()
        
        # plot dataset 1
        p1 = ax.plot(plot_xdata, plot_ydata1, 'b', label=plot_legend1, linewidth=2)
        ax.set_ylabel(plot_ylabel1)
        
        if bgimage != '':
            img = plt.imread(bgimage)
            ax.imshow(img, interpolation='bilinear', 
                        origin='upper', extent=bgsize, aspect='auto')
            ax.set_ylim(bgsize[-2:])
            ax.set_xlim(bgsize[:2])
        
        # plot dataset 2
        if plot_ylabel2 != '':
            ax2 = ax.twinx()
            p2 = ax2.plot(plot_xdata, plot_ydata2, 'r', label=plot_legend2, linewidth=2)
            ax2.set_ylabel(plot_ylabel2)
            
            # legend for two datasets
            ax.legend(handles=p1+p2, loc='best')
            
            # define y range if bgimage is provided
            if bgimage != '' and ylim2 != [0] :
                ax2.set_ylim(ylim2)
            
        else:
            ax.legend(handles=p1, loc='best')
        
        
        plt.subplots_adjust(left=0.15)
        
        return fig, ax
        
        
    def missionDFplotTS(self, dataset1='', METrange=[0], dataset2 = '', 
                        plot_title='', bgimage = '', ylim2=[0], xunit='sec',
                        filename = '', dpi=600):
        """
        Plots one or two datasets from the mission data frame. Y axis labels
        are automatically set in function of chosen dataset. Image can be
        exported.

        Parameters
        ----------
        dataset1 : string
            Name of missionDF dataset to be plotted
        METrange : array of two floats, optional
            Start/End MET values. The default is [0].
        dataset2 : string, optional
            Name of second dataset. The default is ''.
        plot_title : string, optional
            Plot title. The default is ''.
        bgimage : string, optional
            Background image filename/path. The default is ''.
        ylim2 : array of two flaots, optional
            Limitation of right y axis. The default is [0].
        xunit : string, optional
            Unit of X axis. Allowed options: sec, min, hours. The default is 'sec'.
        filename : string, optional
            Where to save the image if required. The default is ''.
        dpi : int, optional
            DPI for saed image. The default is 600.

        Returns
        -------
        fig : figure handle
            No description.
        ax : axis handle
            No description.

        """
      
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        # check validity of first required dataset
        if dataset1=='':
            print('MRSviewer:\tNo dataset name provided. Quitting.')
            return None
        elif not dataset1 in self.missionDF.columns:
            print('MRSviewer:\tNo dataset1 not available in mission DF. Quitting.')
            return None, None
            
        # check valididty of second required dataset
        if dataset2 != '' and not dataset2 in self.missionDF.columns:
            print('MRSviewer:\tNo dataset2 not available in mission DF. Skipping.')
            dataset2 = ''
            
        # set up default title
        if plot_title == '':
            plot_title = self.MD.name + ' ' + dataset1
            if dataset2 != '':
                plot_title += ' & ' + dataset2
            
        # set up default x label
        if xunit == 'sec':
            plot_xlabel = 'seconds [s]'
        elif xunit == 'min':
            plot_xlabel = 'minutes [min]'
        elif xunit == 'hours':
            plot_xlabel = 'hours [h]'
        else:
            print('MRSviewer:\tInvalid unit for x-axis provided. Using seconds.')
            xunit == 'sec'
            plot_xlabel = 'seconds [s]'
        
        
        # set up default y labels
        plot_ylabel1 = dataset1
        plot_ylabel2 = dataset2
        
        # set up default legend labels
        plot_legend1 = dataset1
        plot_legend2 = dataset2
        
        # default x data
        plot_xdata = self.missionDF.loc[DFpointer[0]:DFpointer[1], ['MET']] 
        if xunit == 'min':
            plot_xdata /= 60
        elif xunit == 'hours':
            plot_xdata /= 3600 
        
        
        plot_ydata1 = self.missionDF.loc[DFpointer[0]:DFpointer[1],[dataset1]]
        if dataset2 != '':
            plot_ydata2 = self.missionDF.loc[DFpointer[0]:DFpointer[1],[dataset2]]
        else:
            plot_ydata2 = 0
        
        # check for background image
        if bgimage != '':
            bgsize = [int(s) for s in re.findall(r'\d+', bgimage)]
        
            if len(bgsize)<4:
                print('MRSviewer:\tInvalid background image file name (image size missing). Skipping.')
                bgimage = ''
            else:
                # take only last four values and make an array (was a list)
                bgsize = np.array(bgsize[-4:], dtype=float)
                
                # transform to minutes or hours
                if xunit == 'min':
                    bgsize[:2] /= 60
                elif xunit == 'hours':
                    bgsize[:2] /= 3600 
                
                # check ylim2 provided if needed
                if dataset2 != '':
                    if ylim2 == [0]:
                        print('MRSviewer:\tWarning: no ylim provided for dataset2. No alignment to BG image.')
        else:
            bgsize = 0
     
        # overwrite settings in function of type of data
        if dataset1 == 'EarthFixedVEL':
            plot_ylabel1 = 'SC velocity (Earth-fixed) [m/s]'
        elif dataset1 == 'EarthAlt':
            plot_ylabel1 = 'WGS84 altitude [km]'
            plot_ydata1 /= 1000
        elif dataset1 == 'atmosT':
            plot_ylabel1 = 'Temperature [K]'
        elif dataset1 == 'atmosPa':
            plot_ylabel1 = 'Pressure [kPA]'
            plot_ydata1 /= 1000
        elif dataset1 == 'atmosRho':
            plot_ylabel1 = 'Density [kg/m^3]'
        elif dataset1 == 'atmosM1':
            plot_ylabel1 = 'Mach 1 [m/s]'
        elif dataset1 == 'EarthFixedFPA':
            plot_ylabel1 = 'SC FPA (Earth-fixed) [deg]'
        elif dataset1 == 'EarthFixedHA':
            plot_ylabel1 = 'SC HA (Earth-fixed) [deg]'
        elif dataset1 == 'EarthWGS84FixedFPA':
            plot_ylabel1 = 'SC FPA (Earth-fixed, WGS84) [deg]'
        elif dataset1 == 'EarthGS84FixedHA':
            plot_ylabel1 = 'SC HA (Earth-fixed, WGS84) [deg]'
        elif dataset1 == 'dynPress':
            plot_ylabel1 = 'Dynamic pressure [kPa]'
            plot_ydata1 /= 1000
        elif dataset1 == 'Mach':
            plot_ylabel1 = 'Mach number'
        elif dataset1 == 'AoA':
            plot_ylabel1 = 'Angle of Attack [deg]'
        elif dataset1 == 'SC_active_Thrust':
            plot_ylabel1 = 'Thrust [kN]'
            plot_ydata1 /= 1000
        elif dataset1 == 'SC_active_Mass':
            plot_ylabel1 = 'Mass [1000 kg]'
            plot_ydata1 /= 1000
        elif dataset1 == 'SC_active_acc':
            plot_ylabel1 = 'SC acceleration [m/s]'
        elif dataset1 == 'DragF':
            plot_ylabel1 = 'Drag force [kN]'
            plot_ydata1 /= 1000
        elif dataset1 == 'EarthFPA':
            plot_ylabel1 = 'SC FPA (Earth, inertial) [deg]'
        elif dataset1 == 'EarthHA':
            plot_ylabel1 = 'SC HA (Earth, inertial) [deg]'
        elif dataset1 == 'EarthVEL':
            plot_ylabel1 = 'SC velocity (Earth, inertial) [m/s]'
        elif dataset1 == 'MoonFPA':
            plot_ylabel1 = 'SC FPA (Moon, inertial) [deg]'
        elif dataset1 == 'MoonVEL':
            plot_ylabel1 = 'SC velocity (Moon, inertial) [m/s]'
        elif dataset1 == 'ShadowFunction':
            plot_ylabel1 = 'Sun visibility [%]'
            plot_ydata1 *= 100
        elif dataset1 == 'SRP':
            plot_ylabel1 = 'Solar Radiation Pressure [N]'
        elif dataset1 == 'range':
            plot_ylabel1 = 'Range [km]'
            plot_ydata1 /= 1000
        elif dataset1 == 'Mach':
            plot_ylabel1 = 'Mach number'
        elif dataset1 == 'Mach':
            plot_ylabel1 = 'Mach number'
     
        
        # call plot function 
        fig, ax = self.plotDatasets(plot_title, plot_xlabel, plot_ylabel1, plot_ylabel2,
                         plot_xdata, plot_ydata1, plot_ydata2, plot_legend1, plot_legend2,
                         bgimage, bgsize, ylim2)
        
        # save plot
        if filename != '':
            fig.savefig(filename, dpi=dpi)
        
        
        return fig, ax
    
    def plot_ElevationValues(self, METrange=[0]):
        """
        Plogs various elevation angles.

        Parameters
        ----------
        METrange : array of two floats, optional
            MET start/end range. The default is [0] (no start/end values)

        Returns
        -------
        fig : figure handle
            No description.
        ax : axis handle
            No description.

        """
        
        if not 'gVec_Earth_ENU_abs_elev' in self.missionDF.columns:
            print('MRSviewer:\tNo guidance vectors stored mission DF. Quitting.')
            return None, None
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        plot_title = self.MD.name + ' guidance elevation values'
        plot_xlabel = 'seconds [s]'
        plot_ylabel = 'elevation values [deg]'
        
        fig, ax = plt.subplots(1,1)
        ax.set_title(plot_title)
        ax.set_xlabel(plot_xlabel)
        ax.set_ylabel(plot_ylabel)
        ax.grid()
        
        # same x data for all datasets
        plot_xdata = self.missionDF.loc[DFpointer[0]:DFpointer[1], ['MET']] 
        
        ax.plot(plot_xdata, self.missionDF.loc[DFpointer[0]:DFpointer[1], ['gVec_Earth_ENU_abs_elev']],\
                label='gVec_Earth_ENU_abs_elev', linewidth=2)
        ax.plot(plot_xdata, self.missionDF.loc[DFpointer[0]:DFpointer[1], ['gVec_EFvel_Earth_ENU_delta_elev']],\
                label='gVec_EFvel_Earth_ENU_delta_elev', linewidth=2)
        ax.plot(plot_xdata, self.missionDF.loc[DFpointer[0]:DFpointer[1], ['gVec_SFvel_Earth_ENU_delta_elev']],\
                label='gVec_SFvel_Earth_ENU_delta_elev', linewidth=2)
        ax.plot(plot_xdata, self.missionDF.loc[DFpointer[0]:DFpointer[1], ['gVec_Launch_ENU_abs_pitch']],\
                label='gVec_Launch_ENU_abs_pitch', linewidth=2)
        ax.plot(plot_xdata, self.missionDF.loc[DFpointer[0]:DFpointer[1], ['gVec_GCRF_abs_elev']],\
                label='gVec_GCRF_abs_elev', linewidth=2)
        ax.plot(plot_xdata, self.missionDF.loc[DFpointer[0]:DFpointer[1], ['gVec_VUW_abs_elev']],\
                label='gVec_VUW_abs_elev', linewidth=2)
        ax.plot(plot_xdata, self.missionDF.loc[DFpointer[0]:DFpointer[1], ['gVec_VNB_abs_elev']],\
                label='gVec_VNB_abs_elev', linewidth=2)
        
        ax.legend(loc='best')
    
        #plt.subplots_adjust(left=0.15)
    
        return fig, ax
        
    def print_EventDetails(self, eventList, eventProperties):
        """
        Prints a table of specific values/properties of specific events.

        Parameters
        ----------
        eventList : array of strings
            Array containing the names of events for which properties are shown.
        eventProperties : array of strings
            Array containing the names of eventsDF columns to be shown.

        Returns
        -------
        None.

        """
        
        # check that events are available, otherwise remove them from list
        eventIndexList = np.array([], dtype=np.int8)
        for i, eventName in enumerate(eventList):
            
            # find index of event desc
            eventIndex = self.eventsDF[self.eventsDF['eventType']==eventName].index.values
        
            # if event not found 
            if len(eventIndex)>0:
                eventIndexList = np.append(eventIndexList, int(eventIndex[0]))
       
        # check that properties are available, otherwise remove them from list
        propertiesDeleteIndex = np.array([], dtype=np.int8)
        for i, propertyName in enumerate(eventProperties):
            
            if not propertyName in self.eventsDF.columns:
                propertiesDeleteIndex = np.append(propertiesDeleteIndex, int(i))
        
        # remote invalid properties
        eventProperties = np.delete(eventProperties, propertiesDeleteIndex)
        
        
        # make table and plot events + properties
        tableColNames = np.append(['Event', 'MET [s]'], eventProperties)
        
        t = PrettyTable(list(tableColNames), max_width=25)
        t.float_format = '5.3'
        t.align['Event'] = 'l'   
        #t.hrules = 1
        
        for i, eventIndex in enumerate(eventIndexList):
             
             t.add_row([
                 self.eventsDF.iloc[eventIndex]['eventType'],
                 self.eventsDF.iloc[eventIndex]['MET']
                 ] + list(self.eventsDF.iloc[eventIndex][eventProperties]))
             
        print(t)
        
        
        
        
        
        
    
    def plot_GroundtrackEarth(self, METrange=[0]):
        """
        Ground track on Earth surface.

        Parameters
        ----------
        METrange : array of two floats, optional
            MET start/end range. The default is [0] (no start/end values)

        Returns
        -------
        fig : figure handle
            No description.
        ax : axis handle
            No description.

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
            
        return fig, ax
            
            
    def plot_GroundtrackMoon(self, METrange=[0]):
        """
        Ground track on Moon surface.

        Parameters
        ----------
        METrange : array of two floats, optional
            MET start/end range. The default is [0] (no start/end values)

        Returns
        -------
        fig : figure handle
            No description.
        ax : axis handle
            No description.

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
        
        return fig, ax
        
    def plot_EarthOE(self, METrange=[0]):
        """
        Orbital elements in Earth orbit.

        Parameters
        ----------
        METrange : array of two floats, optional
            MET start/end range. The default is [0] (no start/end values)

        Returns
        -------
        Nfig : figure handle
            No description.
        ax : axis handle
            No description.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        figtitle = self.MD.name + ' Orbital Elements'
        fig, ax = plt.subplots(2,3, figsize=(9, 6))
        fig.suptitle(figtitle, fontsize=16)
       
        
        plt.subplots_adjust(left=0.15,
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
        
        return fig, ax
        
    def plot_MoonOE(self, METrange=[0]):
        """
        Orbital elements in Moon orbit.

        Parameters
        ----------
        METrange : array of two floats, optional
            MET start/end range. The default is [0] (no start/end values)
            
        Returns
        -------
        fig : figure handle
            No description.
        ax : axis handle
            No description.

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
        
        return fig, ax
        
    def plot_ComparisonPosDiff(self, METrange=[0]):
        """
        Absolute distance between simulated trajectory and imported comparison data.

        Parameters
        ----------
        METrange : array of two floats, optional
            MET start/end range. The default is [0] (no start/end values)

        Returns
        -------
        fig : figure handle
            No description.
        ax : axis handle
            No description.

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
        
        return fig, ax
        
    def plot_ComparisonStateVecDiff(self, METrange=[0]):
        """
        Element-wise comparison of state vectors of simulated and imported trajectories.
        

        Parameters
        ----------
        METrange : array of two floats, optional
            MET start/end range. The default is [0] (no start/end values)

        Returns
        -------
        fig : figure handle
            No description.
        ax : axis handle
            No description.

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
        
        return fig, ax
        
        
    def plot_MoonRPS(self, METrange=[0]):
        """
        3D representation of simulated trajectory in Earth-Moon rotating frame.

        Parameters
        ----------
        METrange : array of two floats, optional
            MET start/end range. The default is [0] (no start/end values)

        Returns
        -------
        fig : figure handle
            No description.
        ax : axis handle
            No description.

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
        
        return fig, ax
        
    def plot_GCRF_XYpos(self, METrange=[0], withMoon=0):
        """
        2D plot of GCRF position.

        Parameters
        ----------
        METrange : array of two floats, optional
            MET start/end range. The default is [0] (no start/end values)
        withMoon : TYPE, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        fig : figure handle
            No description.
        ax : axis handle
            No description.

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
        
        return fig, ax
    
    def plot_GCRF_orbit(self, METrange=[0]):
        """
        3D plot of GCRF position including Earth surface.

        Parameters
        ----------
        METrange : array of two floats, optional
            MET start/end range. The default is [0] (no start/end values)

        Returns
        -------
        fig : figure handle
            No description.
        ax : axis handle
            No description.

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
        
        return fig, ax
        
        
    def export_Earth_KML(self, METrange=[0], folder='./'):
        """
        Export KML file of trajetory around Earth.

        Parameters
        ----------
        METrange : array of two floats, optional
            MET start/end range. The default is [0] (no start/end values)
        folder : string, optional
            Folder in which to save the exported KML. Default is current folder (./)

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        # make path to file
        filename = folder + self.MD.name +' Earth_Orbit.kml'
        
        # get coords
        mycoords = self.missionDF.loc[DFpointer[0]:DFpointer[1],
                                      ['EarthLat', 'EarthLon', 'EarthAlt']].to_numpy()
        
        # reduce to 5000 values if needed by interpolation 
        if len(mycoords)>5000:
            reducefactor = int(np.ceil(len(mycoords)/5000))
            mycoords = mycoords[0::reducefactor, :]

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
        METrange : array of two floats, optional
            MET start/end range. The default is [0] (no start/end values)
        folder : string, optional
            Folder in which to save the exported KML. Default is current folder (./)

        Returns
        -------
        None.

        """
        
        # get start/end pointer for mission dataframe
        DFpointer = self.get_DFpointer(METrange)
        
        # make path to file
        filename = folder + self.MD.name +' Moon_Orbit.kml'
        
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
        
    
    