#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 17:01:25 2023

@author: Thibault Bautze-Scherff
"""

import numpy as np
import pandas as pd
import time
import pyshtools as pysh
import datetime
import importlib
import os
import sys
from importlib_resources import files
from scipy.interpolate import CubicSpline, make_interp_spline
from nrlmsise00 import msise_model
import spaceweather as sw
from scipy.integrate import solve_ivp
from skyfield.api import PlanetaryConstants, load, utc, Distance, wgs84, Velocity
from skyfield.positionlib import Geocentric
from skyfield.framelib import itrs, true_equator_and_equinox_of_date
from skyfield.toposlib import ITRSPosition
from skyfield.elementslib import OsculatingElements

from .data import defaultMRSmission
from .MRSrocket import SpaceCraft
from .MRSguidance import Guidance


# general constants
PKGNAME = 'myrocketsimulator'
DEG2RAD = np.pi/180.
RAD2DEG = 180./np.pi
HOURS2DEG = 15.0
DEG2HOURS = 1./15.0
SEC2DAYS = 1./(3600.*24.)
HOURS2RAD = HOURS2DEG * DEG2RAD
RAD2HOURS = RAD2DEG * DEG2HOURS
AU = 149597870700. # astronomical unit [m]
EARTH_ROT_SPEED = 7292115.1467e-11 # [rad/s]
SRP1AU =  4.5344321e-6 #  [N/m^2] solar flux radiation pressure at one AU
R = 287.052874 # specific gas constant dry air

# Radius values for ecplise + tide calculations
SUN_RADIUS = 696340e3 # [m]
EARTH_RADIUS = 6378.1370e3 # [m]
MOON_RADIUS  = 1737.4e3 # [m]

# load Skyfield timescale
ts = load.timescale()

# Moon frames through Skyfield
pc = PlanetaryConstants()
pc.read_text(load('moon_080317.tf'))
pc.read_text(load('pck00008.tpc'))
pc.read_binary(load('moon_pa_de421_1900-2050.bpc'))
MoonFrameME = pc.build_frame_named('MOON_ME_DE421')
MoonFramePA = pc.build_frame_named('MOON_PA_DE421')



class defaultSpacecraft():
    """
    Class for the default spacecraft, used for static propagation.
    Static values are required for drag acceleration
    """

    name = 'MRS default spacecraft'

    class staticValues():
        mass = 500 # [kg]
        dragarea = 1 # [m^2]
        Cd = 0.25
        Cr = 1.8
        SRParea = 1 # [m^2]

class MRSstaticSpacecraft():
    """
    Class for static spacecrafts. Used to load provided data or to make a
    default spacecraft (when no spacecraft data (SCD) was provided).

    """

    def __init__(self, SCD=0):
        """
        Initializes a static spacecraft.

        Parameters
        ----------
        SCD : Spacecraft data Class, optional
            Class that contains spacecraft data, as provided in mission file.
            The default is 0, i.e. no class is provided and default values are
            used.

        Returns
        -------
        None.

        """

        # if spacecraft data is provided
        if SCD:
            self.SCD = SCD
            self.spacecraftname = self.SCD.name
        # otherwise, use default values
        else:
            self.SCD = defaultSpacecraft()
            self.spacecraftname = 'MRS Default Spacecraft'

        # spacecraft mode (active/static)
        self.mode = 'static'

    def get_staticvalues(self):
        # returns previously set values
        return self.SCD.staticValues.mass, \
               self.SCD.staticValues.dragarea, \
               self.SCD.staticValues.Cd, \
               self.SCD.staticValues.Cr, \
               self.SCD.staticValues.SRParea

    def get_ThrustMassOther(self, MET, pa=0.0, verbose='v'):
        # not implemented for static spacecrafts
        # returns static mass with 0 thrust
        return 0, self.SCD.staticValues.mass, np.zeros((1,6))

    def get_DragF(self, MET, vrel=0.0, rho=0.0, mach1=0.0):
        # not implemented for static spacecraft
        # using static values
        staticDrag = .5 * self.get_AreaCd() * vrel**2 * rho
        return staticDrag

    def get_AreaCd(self, MET=0., mach=0.0):
         return self.SCD.staticValues.dragarea * self.SCD.staticValues.Cd

    def set_fixedThrottlePointer(self, MET):
        # not implemented for static spacecrafts
        return None

    def reset_fixedThrottlePointer(self):
        # not implemented for static spacecrafts
        return None

    def get_EventsList(self):
        # make empty list
        self.eventsDF = pd.DataFrame(index = range(0), columns=['MET','eventType'])
        return self.eventsDF


class MRSstaticGuidance():
    """
    Class for static spacecrafts. The returned guidance vector (through 
    get_guidance()) has the same direction as the velocity vector of the 
    provided state vector y.

    """

    def __init__(self):
         """
         Initializes a static guidance object.

         Parameters
         ----------
         None.

         Returns
         -------
         None.

         """
         
         # set name
         self.guidancename = 'No guidance'

         # guidance mode (active/static)
         self.mode = 'static'

    def get_guidance(self, JDnow, y, MET, mode='TrueMET'):
        """

         JDnow : float
            Current Julian Date (TBD); not used.
         y : array of floats
             State vector in GCRF.
         MET : float
             Mission Elapsed Time; not used.
         mode : string, optional
             Mode of guidance determination. The default is 'TrueMET'; not used.

         Returns
         -------
         gVec : array of floats
             Guidance vector in GCRF.

        """

        # guidance direction = velocity direction
        gVec = y[3:]/np.linalg.norm(y[3:])
        
        return gVec



class MRSmission():
    """
    MRSmission is the main class of the MyRocketSimulator software package.
    Based on provided data, it simulate the trajectory of a spaecraft. Later
    on, it can add specific data to the trajectory and export the results.


    Attributes
    ----------
    EarthGravModel : str
        Setting of Earth gravity model (default: EGM96, optional EGM2008)
    EarthGravZonalOnly : int
        0 if using whole spherical harmonics, 1 for zonal harmonics only
    eventCrashed : int
        0 if not crashed, 1 if crashed into Earth
    eventsDF : dataframe
        State vectors and additional information of relevant events
    integrator_atol : float
        Absolute tolerance setting for the integrator (default: 1e-9)
    integrator_max_step : float
        Max. step size setting for the integrator (default: 1000.)
    integrator_rtol : float
        Relative tolerance setting for the integrator (default: 1e-9)
    MD : object
        All provided mission data
    missionDF : dataframe
        Trajectory simulation result and additional data
    SC : object
        Object that represents the loaded spacecraft
    GO: object
        Object that is used to calculate the guidance vector for thrust


    Methods
    -------
    check_MD()
        Sets up runtime variabels & checks the mission data (MD) agains errors
    expand_DF(typelist)
        Adds additiona informatin to the mission dataframe after the simulation
    get_EventsList(eventNames)
        Provides the times and statevectors of defined events
    load_comparisonStateVec(timevec, statevecs)
        Interpolates comparison data and adds it to the missin dataframe
    load_externalMission(timevec, statevecs, name='extMission', atmosModel='')
        Generates a mission dataframe based on provided state vectors
    show_info()
        Displays all user provided data for the missin and some statistics
    return_comparisonStateVec()
        Returns interpolated comparison state vectors
    return_StateVec()
        Returns state vectors from the simulation
    run_mission()
        Performs the simulation and stores the results in the mission dataframe

    """

    def __init__(self, missionname='defaultMRSmission', checkmission=True):
        """
        Initializes the MRS mission object.

        Parameters
        ----------
        missionname : string, optional
             The name of the mission file to be loaded.
             The default is 'defaultMRSmission'.
        checkmission : True/False, optional
            Wether to check the mission data or not.
            The default is True if a mission name is provided

        """


        # default settings, can be changed before checking or running the mission
        self.integrator_atol = 1e-9 # [m]/[m/s]
        self.integrator_rtol = 1e-9
        self.integrator_max_step = 10. # [s]
        self.EarthGravModel = 'EGM96' # alternative: EGM2008
        self.EarthGravZonalOnly = 0
        self.DEephemeris = 'DE421' # alternative: DE440
        self.OE_TEME = 0 # if 0: OE calculated in ICRF, if 1: OE calc. in TEME
        self.fastEphemeris = 0 # if set to 1, planet positions will be fixed 
                               
        # space weather settings
        self.use_spaceweather = 1
        self.f107s = 150
        self.f107  = 150
        self.Ap = 3
        self.transToMSISE90 = 0
        self.RhoGain = 1.07         # exemplary value; used in Bautze-Scherff, 2023
        self.RhoOffset = 3.0495e-14 # exemplary value; used in Bautze-Scherff, 2023


        # load mission if one provided
        if missionname != 'defaultMRSmission':
            # load mission
            self.load_mission(missionname, checkmission)
        else:
            print('MRS:\t\tUsing default mission.')
            self.load_mission('defaultMRSmission', checkmission)


        # overwrite internal default settings if provided in mission profile
        self.overwrite_settings()
        
        # counters; may be removed in later versions
        self.calls_get_acceleration = 0
        self.calls_preload_ObjectPosVel = 0

        return None
    
    def overwrite_settings(self):
        """
        Overwrites internal MRS default settings if provided in missin data.

        Returns
        -------
        None.

        """
        
        if hasattr(self.MD, 'integrator_atol'):
            self.integrator_atol = self.MD.integrator_atol 
        if hasattr(self.MD, 'integrator_rtol'):
            self.integrator_rtol = self.MD.integrator_rtol 
        if hasattr(self.MD, 'integrator_max_step'):
            self.integrator_max_step = self.MD.integrator_max_step 
        if hasattr(self.MD, 'EarthGravModel'):
            self.EarthGravModel = self.MD.EarthGravModel 
        if hasattr(self.MD, 'EarthGravZonalOnly'):
            self.EarthGravZonalOnly = self.MD.EarthGravZonalOnly 
        if hasattr(self.MD, 'DEephemeris'):
            self.DEephemeris = self.MD.DEephemeris 
        if hasattr(self.MD, 'OE_TEME'):
            self.OE_TEME = self.MD.OE_TEME 
        if hasattr(self.MD, 'fastEphemeris'):
            self.fastEphemeris = self.MD.fastEphemeris 
        if hasattr(self.MD, 'use_spaceweather'):
            self.use_spaceweather = self.MD.use_spaceweather 
        if hasattr(self.MD, 'f107s'):
            self.f107s = self.MD.f107s 
        if hasattr(self.MD, 'f107'):
            self.f107 = self.MD.f107 
        if hasattr(self.MD, 'Ap'):
            self.Ap = self.MD.Ap 
        if hasattr(self.MD, 'transToMSISE90'):
            self.transToMSISE90 = self.MD.transToMSISE90 
        if hasattr(self.MD, 'RhoGain'):
            self.RhoGain = self.MD.RhoGain 
        if hasattr(self.MD, 'RhoOffset'):
            self.RhoOffset = self.MD.RhoOffset 
        
        
        

    def load_mission(self, missionname, checkmission=True):
        """
        Loads a mission into the object and checks the validity of the
        mission data.

        Parameters
        ----------
        missionname : string
            Relative path to the missin file (ending with .py)
        checkmission : True/False, optional
             If mission shoudl already be checked. The default is True.

        """


        # save + print missionname
        self.missionname = missionname
        print('MRS:\t\tLoading mission object \''+self.missionname+'\'.')

        # mission variables
        self.MDloaded = 0 # 1 if mission data is loaded
        self.MDvalid = 0 # 1 if mission data is valid


        if missionname != 'defaultMRSmission':
            # load mission data from external file
            try:

                # This code raises an error when changes were made to the
                # mission data file, but it works (i.e. it's really loading
                # the specified file, not just reloading it from somewhere else)
                # Only probelmatic with iPython; solution:
                # execute %autoreload 0 in iPython
                # https://stackoverflow.com/questions/63469147/programmatically-check-for-and-disable-ipython-autoreload-extension
                # https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html

                # https://stackoverflow.com/questions/67631/how-can-i-import-a-module-dynamically-given-the-full-path
                pathtoMRSmissiondata = os.getcwd()+'/'+missionname
                spec = importlib.util.spec_from_file_location('MDmodule', pathtoMRSmissiondata)
                foo = importlib.util.module_from_spec(spec)
                sys.modules['MDmodule'] = foo
                spec.loader.exec_module(foo)
                self.MD = foo.MRSmissionData()

                print('MRS:\t\tLoading MRSmissiondata from file: ')
                print('MRS:\t\t',pathtoMRSmissiondata)

            except ImportError:
                print('MRS:\t\tERROR: Import of specified mission not possible.')
                print('MRS:\t\tLeaving load_mission().')
                return None

            except FileNotFoundError:
                print('MRS:\t\tERROR: File not found:')
                print('MRS:\t\t',pathtoMRSmissiondata)
                print('MRS:\t\tLeaving load_mission().')
                return None

        else:
            self.MD = defaultMRSmission.MRSmissionData


        print('MRS:\t\tMission \''+self.MD.name+'\' loaded.')
        self.MDloaded = 1

        # check mission data validity
        if checkmission:
            # MD not valid (because it returns >0 (=number of errors found))
            if self.check_MD():
                print('MRS:\t\tError: mission data invalid. Leaving load_mission().')
                return None
            else:
                print('MRS:\t\tMission data valid.')
        else:
            print('MRS:\t\tAttention: mission data not checked.')

        return None

    def check_MD(self):
        """
        Checks mission data against errors of all kind (missing, wrong, non
        logical, ...) and sets up internal variables.


        Returns
        -------
        errorfound : float (but could be int)
            Number of found errors in mission data (self.MD)

        """

        # Do not check if already checked, because this may cause errors in
        # calculated tables.
        if self.MDvalid:
            print('MRS:\t\tMission data validity already checked. Leaving check_MD().')
            return 0

        # do not check if MD is not loaded
        if not self.MDloaded:
            print('MRS:\t\tMission data not loaded. Leaving check_MD().')
            return 0

        # welcome message
        print('MRS:\t\tChecking mission data validity.')

        # counter for errors
        errorfound = 0

        # load ephemeris
        self.load_ephemeris()

        # add propagation segments after deltaV segments
        for i in range(len(self.MD.missionSegments)):
            # detect deltaV segments
            if self.MD.missionSegments.type[i]==1:
                # copy propagation row from above
                propagationSegment = \
                    self.MD.missionSegments.loc[[i-1]].reset_index(drop=True)
                # set MET to MET of deltaV maneuver
                propagationSegment.MET = self.MD.missionSegments.MET[i]
                propagationSegment.comment = 'Continue propagation'
                # insert new propagation segment
                self.MD.missionSegments = pd.concat([self.MD.missionSegments.iloc[:i+1],
                                             propagationSegment,
                                             self.MD.missionSegments.iloc[i+1:]]).reset_index(drop=True)

        # early end of simulation
        if self.MD.tend_MET != 0:
            print('MRS:\t\tUsing early ending MET:',self.MD.tend_MET)
            # find mission segments after tend_MET
            indexMissionSegmentsToDrop = \
                self.MD.missionSegments[self.MD.missionSegments.MET>self.MD.tend_MET].index
            # drop if mission ends before end
            if not indexMissionSegmentsToDrop.empty:
                # drop these segments
                self.MD.missionSegments.drop(indexMissionSegmentsToDrop, inplace=True)
                # make new last row
                self.MD.missionSegments.loc[len(self.MD.missionSegments),\
                                            ['MET', 'type', 'configID', 'comment']] = \
                    self.MD.tend_MET, 0, 0, 'Early simulation end (tend_MET)'



        # calc t0_JD if not provided for launchtype 0 (start from state vector)
        # or for launch from pad
        if (self.MD.launchtype==0 and self.MD.t0_JD<=0) or self.MD.launchtype==1:
            # add info to datetime that it's UtC
            self.MD.t0_UTC = self.MD.t0_UTC.replace(tzinfo=utc)
            # updated TDB JD for given UTC time
            self.MD.t0_JD = ts.from_datetime(self.MD.t0_UTC).tdb


        # make sure t0_MET = 0 for launch from pad
        if self.MD.launchtype==1:
            self.MD.t0_MET = 0


        # max degree for Earth gravity SH
        lmaxEarth = np.max(self.MD.forcesSettings['EarthSHn'])
        # max degree for Moon gravity SH
        lmaxMoon = np.max(self.MD.forcesSettings['MoonSHn'])

        # load Earth SH
        if self.EarthGravModel == 'EGM2008':
            self.clmEarth = pysh.datasets.Earth.EGM2008(lmax=lmaxEarth)
        elif self.EarthGravModel == 'EGM96':
            # https://importlib-resources.readthedocs.io/en/latest/using.html
            egm96file = files(PKGNAME+'.data').joinpath('egm96_to360.ascii.txt').as_posix()
            self.clmEarth = pysh.SHCoeffs.from_file(egm96file,
                                               lmax=lmaxEarth,
                                               format='shtools')
            self.clmEarth.gm = 398600441500000.0
            self.clmEarth.r0 = 6.3781363e6
        else:
            print('MRS:\t\tERROR: no valid Earth gravity model selected.')
            errorfound += 1
        
        # back up coefficients
        self.clmEarth.coeffsOrig = self.clmEarth.coeffs * 1.0

        # lood Moon SH
        #clmMoon = pysh.datasets.Moon.GRGM1200B(lmax=lmaxMoon)
        self.clmMoon = pysh.datasets.Moon.GRGM900C(lmax=lmaxMoon)
        # back up coefficients
        self.clmMoon.coeffsOrig = self.clmMoon.coeffs * 1.0

        # keep only zonal SH for Earth if required
        if self.EarthGravZonalOnly:
            self.clmEarth.coeffs[:,:,1:]  = 0


        # load space weather file if needed
        if 'nrlmsise00' in self.MD.forcesSettings['atmosModel'].to_numpy():
            # update space weather if using it
            if self.use_spaceweather:
                try:
                    sw.update_data()
                except ConnectionError:
                    print('MRS:\t\tWARNING: no internet connection to update space weather.')
                    return None
                # load data frame
                self.swDF = sw.sw_daily()

        # events
        self.eventCrashed = 0

        # if launchpad is provided, calc local ENU for time of lifotff (determined by change of propagation mode)
        if isinstance(self.MD.launchsite_LLA, np.ndarray):
            # loop through segments
            for i in range(len(self.MD.missionSegments)):
                # find first occurence of propagationmode = 1
                if self.MD.propaSettings.loc[self.MD.missionSegments.configID[i],'mode'] == 1:
                    # calc JD for given MET
                    self.MD.t0_JD_liftoff = self.MD.t0_JD + self.MD.missionSegments.MET[i] * SEC2DAYS
                    # store lmax of launch segment
                    lmax_launchsegment = self.MD.forcesSettings.EarthSHn[self.MD.propaSettings.loc[self.MD.missionSegments.configID[i],'forcesID']]
                    # do not continue further to look for propmode 1
                    break

            # transform launchsite LLA to GCRF at given time of liftoff
            self.y0_liftoff = self.transform_LLAgeodetic_GCRF(self.MD.t0_JD_liftoff, self.MD.launchsite_LLA)
            # if Earth-SH are used in lauch segment:
            if lmax_launchsegment:
                self.ENU_liftoff = self.get_ENUvec_EarthGravity(self.MD.t0_JD_liftoff, self.y0_liftoff, lmax_launchsegment)
            # else, use simpler ENU calculation (based on ellipsoid of Earth, going through Skyfield)
            else:
                self.ENU_liftoff = self.get_ENUvec_Earth(self.MD.t0_JD_liftoff, self.y0_liftoff, frame='WGS84')

        # no launchsite --> ENU_liftoff = GCRF axes
        else:
            self.ENU_liftoff = np.eye(3)

        #
        # LOAD SPACECRAFT
        #

        # check that spacecraft data is available in mission data
        if hasattr(self.MD, 'spacecraft'):
            # if it got a SC element list, it's a active spacecraft
            if hasattr(self.MD.spacecraft, 'SCelements'):
                print('MRS:\t\tLoading '+self.MD.spacecraft.name+' as active spacecraft.')
                self.SC = SpaceCraft(self.MD.spacecraft)
            # only static data provided (hopefully, otherwhise it crashes later, b/c not checked at the moment)
            else:
                print('MRS:\t\tLoading '+self.MD.spacecraft.name+' as static spacecraft.')
                self.SC = MRSstaticSpacecraft(self.MD.spacecraft)
        # no spacecraft provided
        else:
            # make a default spacecraft
            self.SC = MRSstaticSpacecraft()
            
        #
        # LOAD GUIDANCE
        #
        
        # check that gudiacne data is available in mission data
        if hasattr(self.MD, 'guidanceData'):
            print('MRS:\t\tLoading guidance object '+self.MD.guidanceData.name+'.')
            self.GO = Guidance(self.MD.guidanceData)
            
        # no guidance provided
        else:
            self.GO = MRSstaticGuidance()
            
        # set up ENU liftoff
        self.GO.ENU_liftoff = self.ENU_liftoff * 1.0
        

        #
        # LOOK FOR ERRORS IN MISSION DATA
        #

        # check data related to launch from statevector (launchtype 0)
        if self.MD.launchtype==0:

            # MET offset not within range of mission segments
            if self.MD.t0_MET<self.MD.missionSegments.iloc[0]['MET'] or self.MD.t0_MET>=self.MD.missionSegments.iloc[-1]['MET']:
                print('MRS:\t\tERROR: t0_MET is outside of mission segments.')
                errorfound += 1

        # check if active spacecraft is available if required in mission data
        # force settings
        if self.MD.forcesSettings.activeSC.to_numpy().sum()>0:
            if not self.SC.mode == 'active':
                print('MRS:\t\tERROR: MD forcesettigs require active SC.')
                #errorfound += 1


        # store mission data is valid
        if not errorfound:
            self.MDvalid = 1

        return errorfound


    def load_ephemeris(self):
        """
        Loads the GM values and the DE-files into internal variables.
        The DE-version is specified in self.DEephemeris.


        Returns
        -------
        None.

        """

        if self.DEephemeris == 'DE440':
            planetsEphemeris = load('de440.bsp')
            # gravity constants DE440
            self.SUN_GM   = 132712440041.279419 * 1000**3
            self.VENUS_GM = 324858.592000 * 1000**3
            self.EARTH_GM = 398600.435507 * 1000**3
            self.MOON_GM  = 4902.800118 * 1000**3
            self.MARS_GM  = 42828.375816 * 1000**3
            self.JUPITER_GM = 126712764.100000 * 1000**3

        # default is DE421
        else:
            planetsEphemeris = load('de421.bsp')
            # gravity constants DE421
            self.SUN_GM   = 132712440040.944 * 1000**3
            self.VENUS_GM = 324858.592000 * 1000**3
            self.EARTH_GM = 398600.436233 * 1000**3
            self.MOON_GM  = 4902.800076 * 1000**3
            self.MARS_GM  = 42828.375214 * 1000**3
            self.JUPITER_GM = 126712764.800000 * 1000**3

        self.sun = planetsEphemeris['Sun']
        self.venus = planetsEphemeris['VENUS']
        self.earth = planetsEphemeris['Earth']
        self.moon = planetsEphemeris['Moon']
        self.mars = planetsEphemeris['MARS_BARYCENTER']
        self.jupiter = planetsEphemeris['JUPITER_BARYCENTER']

        return None

    def run_mission(self):
        """
        Executes the mission. Resulting state vectors are saved in the
        mission dataframe (self.missionDF).

        Returns
        -------
        None

        """

        # run only if MD is validated
        if not self.MDvalid:
            print('MRS:\t\tERROR: mission data not validated. Exiting.')
            return None

        # start timer
        tic = time.time()

        # show some mission parameters
        print('MRS:\t\tRunning mission '+self.MD.name+'.')

        # get state vector if starting from it (when launchtype==0)
        if self.MD.launchtype == 0:
            self.statevec = self.MD.y0
            self.JD = self.MD.t0_JD
        # generate empty statevector if starting from lauch pad
        else:
            self.statevec = np.zeros(6)
            self.JD = 0

        # make empty DF to which segment-dataframes will be appended
        self.make_TempDF(int(0))
        self.missionDF = self.TempDF
        del self.TempDF

        # make empty DF to append events
        self.make_TempDF(int(0))
        self.eventsDF = self.TempDF
        del self.TempDF
        # add empty colum (needs to be present because accessed later on)
        self.eventsDF['eventType'] = None

        for i in range(self.MD.missionSegments.shape[0]):

            # check for crash
            if self.eventCrashed:
                # end mission
                break

            # don't show message if in last segment (which is not processed)
            if i != len(self.MD.missionSegments)-1:
                print('MRS:\t\tProcessing mission segment {}.'.format(i))

            ###
            ### PROPAGATE SPACECRAFT
            ###

            # if in propgation segment
            if self.MD.missionSegments['type'][i]==0:

                # save current segment and its properties
                self.segmentType = 0
                self.segmentID = i
                self.propaID   = self.MD.missionSegments['configID'][i]
                self.propaMode = self.MD.propaSettings['mode'][self.propaID]
                self.forcesID  = self.MD.propaSettings['forcesID'][self.propaID]
                self.dragOn    = self.MD.forcesSettings['drag'][self.forcesID]
                self.SRPOn     = self.MD.forcesSettings['SRP'][self.forcesID]
                self.activeSC  = self.MD.forcesSettings['activeSC'][self.forcesID]
                self.planets   = self.MD.forcesSettings['planets'][self.forcesID]
                self.EarthSHn  = self.MD.forcesSettings['EarthSHn'][self.forcesID]
                self.MoonSHn   = self.MD.forcesSettings['MoonSHn'][self.forcesID]
                self.EarthTides = self.MD.forcesSettings['EarthTides'][self.forcesID]
                self.MoonTides  = self.MD.forcesSettings['MoonTides'][self.forcesID]
                
                # make sure objects are preloaded at least once, if needed
                self.deactivate_preload_ObjectPosVel = 0

                ###
                ### MAKE MET (MISSION ELLAPSED TIME) VECTOR
                ###

                # if in last segment, just one single MET value
                if i==len(self.MD.missionSegments)-1:
                    # get MET value
                    self.METvec = np.array([self.MD.missionSegments['MET'][i]])

                # if not in last segment, generate MET vector
                else:
                    # if starting from statevector check t0_MET
                    if self.MD.launchtype == 0:
                        # if reference time of state vector is in next segment
                        if self.MD.t0_MET >= self.MD.missionSegments['MET'][i+1]:
                            # go to next segment
                            continue
                        # if reference time is before current segment start time
                        if self.MD.t0_MET < self.MD.missionSegments['MET'][i]:
                            METsegmentstart = self.MD.missionSegments['MET'][i]
                        # otherwise set starttime to provided t0_MET
                        else:
                            METsegmentstart = self.MD.t0_MET
                    # not starting from statevector, use time from segment
                    else:
                        METsegmentstart = self.MD.missionSegments['MET'][i]

                    # make timevector (MET) for segment
                    # starttime: starttime of segment
                    # endtime: starttime of following segment
                    # step size: taken from propagation settings for this segment
                    self.METvec = np.arange(
                        METsegmentstart,
                        self.MD.missionSegments['MET'][i+1],
                        self.MD.propaSettings['stepsizePropa'][self.propaID])
                    # add starttime of next segment as final value of current segment
                    self.METvec = np.append(self.METvec,self.MD.missionSegments['MET'][i+1])



                ###
                ### PREPARE DATA LOGGING
                ###

                # make vector when to log
                self.LOGvec = np.zeros(len(self.METvec))
                self.LOGvec[0::int(self.MD.propaSettings['downsampleLog'][self.propaID])] = 1

                # make sure not to deactivate the logging of last segment
                # (only one value logged), therefore >1
                if len(self.LOGvec)>1:
                    self.LOGvec[-1] = 0 # last one is not logged, will be taken
                                        # care of in following segment

                # make temporary dataframe for logging
                self.make_TempDF(int(np.sum(self.LOGvec)))

                # write segment properties to dataframe
                self.TempDF['segmentID'] = self.segmentID
                self.TempDF['segmentType'] = self.segmentType
                self.TempDF['configID'] = self.propaID
                self.TempDF['propaMode'] = self.propaMode
                self.TempDF['forcesID'] = self.forcesID
                self.TempDF['activeSC'] = self.activeSC

                ###
                ### LOAD SPACE WEATHER (for date at segment start)
                ###

                # check space weather is needed
                if self.MD.forcesSettings['atmosModel'][self.forcesID] == 'nrlmsise00':

                    # JD_TBD at segment starts
                    JD_SegStart = self.MD.t0_JD + \
                        (self.MD.missionSegments['MET'][self.segmentID] - self.MD.t0_MET) * SEC2DAYS
                    # date string
                    datestring = ts.tdb_jd(JD_SegStart).utc_strftime(format='%Y-%m-%d')

                    # update values if using space weather
                    if self.use_spaceweather:
                        # get indices from spaceweather dataframe
                        indices = self.swDF.loc[datestring]
                        # save values
                        self.f107s = indices['f107_81lst_adj']
                        self.f107  = indices['f107_adj']
                        self.Ap    = indices['Apavg']

                    # save as current atmosModel
                    self.atmosModel = 'nrlmsise00'
                    self.TempDF['atmosModel'] = self.atmosModel

                # no atmos model defined
                else:
                    self.atmosModel = '-'
                    self.TempDF['atmosModel'] = self.atmosModel

                ###
                ### GRAVITY
                ###

                # make sure to reset the Cnm/Snm for Moon gravity in case
                # tides are used (may have altered values in previous segment)
                if self.EarthTides:
                    self.clmEarth.coeffs = self.clmEarth.coeffsOrig * 1.0

                # make sure to reset the Cnm/Snm for Moon gravity in case
                # tides are used (may have altered values in previous segment)
                if self.MoonTides:
                    self.clmMoon.coeffs = self.clmMoon.coeffsOrig * 1.0

                ###
                ### PROPAGATION
                ###

                # in last segment, don't care of propagation mode, just fill
                # in last statevec & exit
                if i==len(self.MD.missionSegments)-1:

                    # write data to TempDF
                    self.TempDF['JD_TBD'] = self.JD
                    self.TempDF['MET'] = self.METvec[0]
                    self.TempDF[['x','y','z','vx','vy','vz']] = self.statevec

                    # overwrite segment ID (otherwise it's a new segment and
                    # leads to bad visualization with the MRSvislib
                    self.TempDF['segmentID'] -= 1

                    # finalize missionDF
                    self.missionDF = pd.concat([self.missionDF, self.TempDF], ignore_index=True)

                    # exit segment loop - mission accomplished!
                    break


                # standing at launch pad
                if self.propaMode==0:

                    # propagate launchsite position in GCRF
                    self.propagate_Mode0()

                # free run mode propagation
                elif self.propaMode==1:

                    # propagate without thrust/guidance control
                    self.propagate_Mode1()

                # step-wise propagation
                elif self.propaMode==2:

                    # propagate SC
                    self.propagate_Mode2()

                # unknown mode
                else:
                    print('MRS:\t\tERROR: Unknown propagation mode. Exiting \
                    segment processing.')
                    # exit segment loop
                    break


                ###
                ### POST-PROPAGATION IN-SEGMENT TODOS
                ###

                # append segmentDF to missionDF
                self.missionDF = pd.concat([self.missionDF, self.TempDF], ignore_index=True)

                


            ###
            ### APPLY DELTA-V MANEUVER
            ###

            elif self.MD.missionSegments['type'][i]==1:

                self.maneuverID = self.MD.missionSegments['configID'][i]
                self.maneuverFrame = self.MD.maneuverSettings['frame'][self.maneuverID]
                self.dv = self.MD.maneuverSettings.loc[self.maneuverID,['dx','dy','dz']].to_numpy().astype('float')
                self.planet = self.MD.maneuverSettings['planet'][self.maneuverID]
                self.args = self.MD.maneuverSettings['args'][self.maneuverID]

                if self.maneuverFrame == 'LVLH':
                    self.statevec = self.apply_deltaV_LVLH(self.JD, self.statevec, self.dv, planet=self.planet)
                elif self.maneuverFrame == 'VNB':
                    self.statevec = self.apply_deltaV_VNB(self.JD, self.statevec, self.dv, planet=self.planet)
                elif self.maneuverFrame == 'GCRF':
                    self.statevec[3:] += self.dv
                else:
                    print('MRS:\t\tERROR: Unknown maneuver type. Exiting segment processing.')
                    # exit segment loop
                    break


            # unknown segment type
            else:
                print('MRS:\t\tERROR: Unknown segment type. Exiting segment processing.')
                # exit segment loop
                break


        ###
        ### POST MISSION TODOS
        ###

        # delete temp vars
        self.del_tempVars()

        # end timer
        self.toc = time.time() - tic

        # print final information
        print('MRS:\t\tMission ended. Processing time: {} seconds.'.format(round(self.toc,3)))

        # save kind of missionDF (1=MRS sim, 2=external data)
        self.missionDFtype = 1
    
        return None

    def del_tempVars(self):
        """
        Internal function.
        Deletes temporary attributes of MRSmission object.

        Returns
        -------
        None

        """

        if hasattr(self, 'TempDF'):
            del self.TempDF
        if hasattr(self, 'LOGvec'):
            del self.LOGvec
        if hasattr(self, 'METvec'):
            del self.METvec
        if hasattr(self, 'propaID'):
            del self.propaID
        if hasattr(self, 'segmentID'):
            del self.segmentID
        if hasattr(self, 'forcesID'):
            del self.forcesID
        if hasattr(self, 'dragOn'):
            del self.dragOn
        if hasattr(self, 'activeSC'):
            del self.activeSC
        if hasattr(self, 'planets'):
            del self.planets
        if hasattr(self, 'planet'):
            del self.planet
        if hasattr(self, 'EarthSHn'):
            del self.EarthSHn
        if hasattr(self, 'MoonSHn'):
            del self.MoonSHn
        if hasattr(self, 'atmosModel'):
            del self.atmosModel
        if hasattr(self, 'args'):
            del self.args
        if hasattr(self, 'maneuverID'):
            del self.maneuverID
        if hasattr(self, 'maneuverFrame'):
            del self.maneuverFrame
        if hasattr(self, 'dv'):
            del self.dv
        if hasattr(self, 'SunPosGCRF'):
            del self.SunPosGCRF
        if hasattr(self, 'VenusPosGCRF'):
            del self.VenusPosGCRF
        if hasattr(self, 'MoonPosGCRF'):
            del self.MoonPosGCRF
        if hasattr(self, 'MoonVelGCRF'):
            del self.MoonVelGCRF
        if hasattr(self, 'MarsPosGCRF'):
            del self.MarsPosGCRF
        if hasattr(self, 'JupiterPosGCRF'):
            del self.JupiterPosGCRF

        return None

    def propagate_Mode0(self):
        """
        Internal function.
        Propagates the state vector of the spacecraft standing on its launch
        pad on earth by taking its coordinates and saving them in GCRF.

        Returns
        -------
        None

        """

        # number of MET times required to be calculated
        numSteps = self.METvec.shape[0]

        # pointer in TempDF
        TempDFpointer = 0

        # loop through steps
        for i in range(numSteps):

            # current JD
            JDnow = self.MD.t0_JD + (self.METvec[i] - self.MD.t0_MET) * SEC2DAYS

            # get current statevec
            statevec = self.transform_LLAgeodetic_GCRF(JDnow, self.MD.launchsite_LLA)

            # call spacecraft-thrust-function to display its output
            _, _, _ = self.SC.get_ThrustMassOther(self.METvec[i], 0, verbose='v')

            # save to logfile if needed
            if self.LOGvec[i]==1:
                self.TempDF.loc[TempDFpointer, ['MET']] = self.METvec[i]
                self.TempDF.loc[TempDFpointer, ['JD_TBD']] = JDnow
                self.TempDF.loc[TempDFpointer, ['x','y','z','vx','vy','vz']] = statevec
                # increase pointer for next row
                TempDFpointer += 1

        # update state final vector
        self.statevec = statevec
        self.JD = self.MD.t0_JD + (self.MD.missionSegments['MET'][self.segmentID+1] - self.MD.t0_MET) * SEC2DAYS

        return None



    def propagate_Mode1(self):
        """
        Internal function.
        Propagates the spaecraft within the given mission segment; auto-step
        size control by the given integrator.
        To be used with passive spacecrafts.
        Not recommended for active spacecrafts with thrust and guidance.
        Resulting state vectors are stored to TempDF and (in run_mission())
        to the missionDF.


        Returns
        -------
        None

        """

        # special function to serve as crash event on Earth for integrator termination
        def EndFlightEvent(t, y):
            return self.get_EarthAlt(t, y)
        # end mission if crash occurs
        EndFlightEvent.terminal  = True
        EndFlightEvent.direction = -1.

        # prepare solve_ivp variables
        fun = self.get_slopes
        t_span = np.array([self.MD.missionSegments['MET'][self.segmentID],
                           self.MD.missionSegments['MET'][self.segmentID+1]])
        y0 = self.statevec
        method = self.MD.propaSettings['method'][self.propaID] # get from dataframe
        dense_output = True
        max_step = self.integrator_max_step
        atol = self.integrator_atol
        rtol = self.integrator_rtol
        first_step = 3

        # perform solve_ivp
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
        sol = solve_ivp(fun, t_span, y0, method, dense_output = dense_output,
                        events=EndFlightEvent, first_step=first_step,
                        max_step=max_step, atol=atol, rtol=rtol)

        # METs at which logging is required; don't consider last MET, as not indexed by LOGvec
        METlog = self.METvec[self.LOGvec==1]

        # update final state vector
        self.statevec = sol.y.T[-1] # last computed state vector
        self.JD = self.MD.t0_JD + (sol.t.T[-1]  - self.MD.t0_MET) * SEC2DAYS

        # check if rocket crashed into Earth
        if sol.t_events[0].size>0:
            print('MRS:\t\tWARNING: Collision with Earth, mission ended.')
            # delete METlog timestampe after crash
            METlog = METlog[METlog<sol.t_events[0][0]]
            # add MET for crash
            METlog = np.append(METlog, sol.t_events[0][0])
            # adjust length of TempDF
            self.TempDF = self.TempDF.head(len(METlog))
            # save crash
            self.eventCrashed = 1
            # append to eventsDF
            self.add_event(self.JD, self.statevec, 'EarthCollision')

        # get dense data for required METs
        sol_statevecs = sol.sol(METlog)

        # save to TempDF
        self.TempDF['MET'] = METlog
        self.TempDF['JD_TBD'] = self.MD.t0_JD + (METlog - self.MD.t0_MET) * SEC2DAYS
        self.TempDF[['x','y','z','vx','vy','vz']] = sol_statevecs.T

        self.sol = sol

        return None

    def propagate_Mode2(self):
        """
        Internal function.
        Propagates the spaecraft within the given mission segment; fixed step
        size (provided in misison data propaSettings)
        To be used with active spacecrafts.
        Not recommended for static spacecrafts.
        Resulting state vectors are stored to TempDF and (in run_mission())
        to the missionDF.


        Returns
        -------
        None

        """

        # save first statevector
        self.TempDF.loc[0,['MET']] = self.METvec[0]
        self.TempDF.loc[0,['JD_TBD']] = self.MD.t0_JD + (self.METvec[0] - self.MD.t0_MET) * SEC2DAYS
        self.TempDF.loc[0,['x','y','z','vx','vy','vz']] = self.statevec

        # special function to serve as crash event on Earth for integrator termination
        def EndFlightEvent(t, y):
            return self.get_EarthAlt(t, y)
        # end mission if crash occurs
        EndFlightEvent.terminal  = True
        EndFlightEvent.direction = -1.

        # prepare fixed solve_ivp variables
        fun = self.get_slopes
        #statevec = self.statevec
        method = self.MD.propaSettings['method'][self.propaID]
        dense_output = True
        max_step = self.integrator_max_step
        atol = self.integrator_atol
        rtol = self.integrator_rtol
        first_step = 3 # TODO: needs to be set to value equal (or less) to step size of mode 2


        # number of MET times required to be calculated
        numSteps = self.METvec.shape[0]

        # pointer in TempDF; skip first row because alrady written
        TempDFpointer = 1

        # loop through steps
        for i in range(numSteps-1):

            # prepare solve_ivp variables
            t_span = np.array([self.METvec[i], self.METvec[i+1]])
            y0 = self.statevec

            # set fixed pointer for spacecraft
            self.SC.set_fixedThrottlePointer(self.METvec[i])
            
            # update guidance
            _ = self.GO.get_guidance(self.MD.t0_JD + (self.METvec[0] - self.MD.t0_MET) * SEC2DAYS,\
                                     y, self.METvec[0], mode='TrueMET')

            # perform solve_ivp
            # https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
            sol = solve_ivp(fun, t_span, y0, method, dense_output = dense_output,
                            events=EndFlightEvent, first_step=first_step,
                            max_step=max_step, atol=atol, rtol=rtol)

            # update state vector
            self.statevec = sol.y.T[-1] # last computed state vector
            self.JD = self.MD.t0_JD + (sol.t.T[-1]  - self.MD.t0_MET) * SEC2DAYS

            # check if rocket crashed into Earth
            if sol.t_events[0].size>0:
                print('MRS:\t\tWARNING: Collision with Earth, mission ended.')
                # save crash
                self.eventCrashed = 1
                # append to eventsDF
                self.add_event(self.JD, self.statevec, 'EarthCollision')
                # shorten TempDF
                self.TempDF = self.TempDF.head(TempDFpointer+1)

            # save to logfile if needed; always one ahead, because the statevector for the following MET is calculated
            if self.LOGvec[i+1]==1 or self.eventCrashed==1:
                self.TempDF.loc[TempDFpointer, ['MET']] = sol.t.T[-1]
                self.TempDF.loc[TempDFpointer, ['JD_TBD']] = self.JD
                self.TempDF.loc[TempDFpointer, ['x','y','z','vx','vy','vz']] = self.statevec
                # increase pointer for next row
                TempDFpointer += 1

            # exit in case of crash
            if self.eventCrashed == 1:
                # exit steps-loop
                break

        # clear spacecraft pointer
        self.SC.reset_fixedThrottlePointer()

        return None

    def get_EarthAlt(self, MET, y):
        """
        Internal function.
        Returns the spacecraft's alititude w.r.t. to the WGS84 Earth.

        Parameters
        ----------
        MET : float
            Mission Elapsed Time
        y : array of floats
            State vector

        Returns
        -------
        float
            Altitude [m]

        """

        # calculate current JD_TBD
        JDnow = self.MD.t0_JD + (MET - self.MD.t0_MET) * SEC2DAYS

        # get geodetic coordinates (lat [°], lon [°], alt [m])
        LLAgeodetic = self.transform_GCRF_LLAgeodetic(JDnow, y[:3])

        return LLAgeodetic[2]


    def get_slopes(self, MET, y):
        """
        Internal function.
        Returns the slopes for given state vector.
        Only to be called by integrator, as get_acceleration is called in mode
        'intermediateMET'.

        Parameters
        ----------
        MET : float
            Mission Elapsed Time
        y : array of floats
            State vector

        Returns
        -------
        Array of floats
            Slopes of state vector.

        """

        return np.hstack((y[3:6], self.get_acceleration(MET, y, mode='intermediateMET')))


    def get_acceleration(self, MET, y, mode='TrueMET'):
        """
        Internal function.
        Calculates the acceleration of the spacecraft (i.e. the slopes of
        its velocity vector) by considering various forces:
            - gravity
            - thrust
            - drag

        Parameters
        ----------
        MET : float
            Mission Elapsed Time
        y : array of floats
            State vector

        Returns
        -------
        acceleration : array of floats
            Acceleration in the GCRF frame.

        """
        
        # increase call counter
        self.calls_get_acceleration += 1

        # calculate current JD_TBD
        JDnow = self.MD.t0_JD + (MET - self.MD.t0_MET) * SEC2DAYS

        # prevent preloading of objects if no object positions are needed
        if (self.planets==['Earth'] or self.planets == []) and \
            self.EarthTides == 0 and \
            self.SRPOn == 0:
                self.deactivate_preload_ObjectPosVel = 1

        # preload pos/vel of relevant celestial bodies 
        # always happens at first call, even if fastEphemeris is set to 1
        if not self.deactivate_preload_ObjectPosVel:
            self.preload_ObjectPosVel(JDnow)
        
        # if fastEphemeris are set to 1, then set flag to not reload ephemeris
        # for this segment
        if self.fastEphemeris == 1:
            self.deactivate_preload_ObjectPosVel = 1
        
        # set preload flag for all functions
        preloadFlag = 1 

        # get acceleration from planets
        accPlanets = self.get_accPlanets(MET, y, JDnow, preloadedObjects=preloadFlag)

        # static values are always loaded
        staticMass, staticDragArea, staticCd, staticCr, staticSRParea\
                = self.SC.get_staticvalues()


        # load values in function of active/static SC + Thrust calculation
        if self.activeSC:

            # get temperature/pressure/density/M1
            atmosvalues = self.get_atmos(JDnow, y[:3])

            # get spacecraft properties
            thrustForce, SCmass, _ = self.SC.get_ThrustMassOther(MET, \
                                                   atmosvalues[1], verbose='v')

            # get guidance
            thrustVector = self.GO.get_guidance(JDnow, y, MET, mode=mode)

            # acceleration by TVC
            accThrust = thrustVector * thrustForce / SCmass
        else:
            accThrust = 0

        # calculate drag acceleration if required
        if self.dragOn:
            # if active spacecraft:
            if self.activeSC:

                # get numbers needed to calc drag
                vrel = self.get_relVelocityEarth(JDnow, y)
                rho = atmosvalues[2]
                mach = atmosvalues[3]

                # calc drag acceleration
                accDrag = self.SC.get_DragF(MET, vrel, rho, mach) / SCmass

            else:
                # dynPress comes as a vector in direction of velocity
                accDrag = - self.get_dynPress(MET, y, JDnow) * \
                         staticDragArea * staticCd  / staticMass
        else:
            accDrag = np.zeros(3)


        # Calculate solar radiation pressure acceleration if required
        # Only static SC values considered at the moment.
        if self.SRPOn:

            # get the visibility factor and the vector pointing from SC to Sun.
            v, SCtoSun, occultingBodies = \
                self.get_sunVisibility(JDnow, y[:3], preloadedObjects=preloadFlag)

            # get SRP force
            SRPforce = self.get_SRPForce(v, staticCr, staticSRParea, SCtoSun)

            # acceleration by SRP
            accSRP = SRPforce / staticMass

        else:
            accSRP = np.zeros(3)


        # sum up accelerations
        acceleration = accPlanets + accDrag + accThrust + accSRP

        return acceleration

    def preload_ObjectPosVel(self, JDnow):
        """
        Internal function.
        Loads positions of required planets (self.planets) into temporary
        variables, as well as the position and velocity of the Moon.
        Used to speed up the overall computation by calling Skyfield planet-
        positions only once.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)

        Returns
        -------
        None.

        """
        
        # count function calls
        self.calls_preload_ObjectPosVel += 1

        # get position & velocity of Earth
        EarthJDnow = self.earth.at(ts.tdb_jd(JDnow))
        EarthPosICRF = EarthJDnow.position.m
        EarthVelICRF = EarthJDnow.velocity.m_per_s

        # loop through planets
        for planet in self.planets:
            if planet=='Sun':
                objectPos = self.sun.at(ts.tdb_jd(JDnow)).position.m
                self.SunPosGCRF = objectPos - EarthPosICRF
            elif planet=='Venus':
                objectPos = self.venus.at(ts.tdb_jd(JDnow)).position.m
                self.VenusPosGCRF = objectPos - EarthPosICRF
            elif planet=='Earth':
                continue
            elif planet=='Moon':
                MoonNow = self.moon.at(ts.tdb_jd(JDnow))
                objectPos = MoonNow.position.m
                objectVel = MoonNow.velocity.m_per_s
                self.MoonPosGCRF = objectPos - EarthPosICRF
                self.MoonVelGCRF = objectVel - EarthVelICRF
            elif planet=='Mars':
                objectPos = self.mars.at(ts.tdb_jd(JDnow)).position.m
                self.MarsPosGCRF = objectPos - EarthPosICRF
            elif planet=='Jupiter':
                objectPos = self.jupiter.at(ts.tdb_jd(JDnow)).position.m
                self.JupiterPosGCRF = objectPos - EarthPosICRF
            else:
                print('MRS:\t\tERROR: invalid planet name provided in forcesSettings.')
                continue

        # check if tides needed
        if self.EarthTides or self.MoonTides:
            useTides = 1
        else:
            useTides = 0

        # load Sun if needed for SRP or Tides
        if not 'Sun' in self.planets and (self.SRPOn or useTides):
            objectPos = self.sun.at(ts.tdb_jd(JDnow)).position.m
            self.SunPosGCRF = objectPos - EarthPosICRF
            
        # load Moon if needed for tides
        if not 'Moon' in self.planets and useTides:
            objectPos = self.moon.at(ts.tdb_jd(JDnow)).position.m
            self.MoonPosGCRF = objectPos - EarthPosICRF

        return None

    def get_SRPForce(self, v, Cr, SRParea, SCtoSun):
        """
        Internal function.
        Calculates force caused by Solar Radiation Pressure (SRP).

        Parameters
        ----------
        v : float
            Visibility of the Sun (0-1)
        Cr : float
            Radiation pressure coefficient
        SRParea : float
            Area affected by SRP (in sqaure meters)
        SCtoSun : array of floats
            Vector (Sun w.r.t. to spacecraft)

        Returns
        -------
        SRPforce : array of floats
            SRP force on spaecraft

        """

        SRPforce = - v * Cr * SRParea * SRP1AU * SCtoSun \
                     * AU**2 / np.linalg.norm(SCtoSun)**3

        return SRPforce

    def get_accPlanets(self, MET, y, JDnow, preloadedObjects=0):
        """
        Internal function.
        Returns the cummulated acceleration caused by planets.
        Selection of relevant planets is made mission data forceSettings.


        Parameters
        ----------
        MET : float
            Mission Elapsed Time
        y : array of floats
            State vector
        JDnow : float
            Current Julian Date (TBD)
        preloadedObjects : int, optional
            Whether to compute the Moon pos/vel while execution (0) or use
            preloaded values (faster.)

        Returns
        -------
        accPlanets : array of floats
            Acceleration by planet gravities in the GCRF frame.

        """

        # get current spacecraft position
        SCpos = y[:3]

        # accumulator for planet acceleration
        accPlanets = np.zeros(3)

        # get position of Earth if not using preloaded GCRF position of objects
        if not preloadedObjects:
            EarthPos = self.earth.at(ts.tdb_jd(JDnow)).position.m

        # loop through planets
        for planet in self.planets:

            # get positions (in ECI) and GM-values of planets
            if planet=='Sun':
                GM = self.SUN_GM
                if preloadedObjects:
                    objectToEarth = self.SunPosGCRF
                else:
                    objectPos = self.sun.at(ts.tdb_jd(JDnow)).position.m

            elif planet=='Venus':
                GM = self.VENUS_GM
                if preloadedObjects:
                    objectToEarth = self.VenusPosGCRF
                else:
                    objectPos = self.venus.at(ts.tdb_jd(JDnow)).position.m

            elif planet=='Earth':
                GM = self.EARTH_GM

            elif planet=='Moon':
                GM = self.MOON_GM
                if preloadedObjects:
                    objectToEarth = self.MoonPosGCRF
                else:
                    objectPos = self.moon.at(ts.tdb_jd(JDnow)).position.m

            elif planet=='Mars':
                GM = self.MARS_GM
                if preloadedObjects:
                    objectToEarth = self.MarsPosGCRF
                else:
                    objectPos = self.mars.at(ts.tdb_jd(JDnow)).position.m

            elif planet=='Jupiter':
                GM = self.JUPITER_GM
                if preloadedObjects:
                    objectToEarth = self.JupiterPosGCRF
                else:
                    objectPos = self.jupiter.at(ts.tdb_jd(JDnow)).position.m
            else:
                # empty entry or invalid planet
                continue

            # calculate acceleration
            if planet == 'Earth':
                # get degree of SH coefficients
                lmax = self.EarthSHn
                # get SH-gravity if lmax provided
                if lmax:
                    objectAcc = self.get_EarthGravity(JDnow, y[:3], lmax, preloadedObjects=preloadedObjects)
                # if no lmax, use point gravity
                else:
                    r_norm = np.linalg.norm(SCpos)
                    objectAcc = - GM * SCpos/r_norm**3

            elif planet == 'Moon':
                # get degree of SH coefficients
                lmax = self.MoonSHn
                # get SH-gravity if lmax provided
                if lmax:
                    objectAcc = self.get_MoonGravity(JDnow, y, lmax, preloadedObjects=preloadedObjects)
                # if no lmax, use point gravity
                else:
                    if not preloadedObjects:
                        objectToEarth = objectPos - EarthPos
                    objectToSat =  objectToEarth - SCpos
                    objectAcc = self.get_ThirdBodyAcc(objectToEarth, objectToSat, SCpos, GM)

            else:
                if not preloadedObjects:
                    objectToEarth = objectPos - EarthPos
                objectToSat =  objectToEarth - SCpos
                objectAcc = self.get_ThirdBodyAcc(objectToEarth, objectToSat, SCpos, GM)

            # add up accelerations
            accPlanets += objectAcc

        return accPlanets

    def get_ThirdBodyAcc(self, bodyToEarth, bodyToSC, SCtoEarth, GM):
        """
        Internal function.
        Returns the acceleratin of a celestial body in GCRF.

        Parameters
        ----------
        bodyToEarth : array of floats
            Position of celestial body in GCRF
        bodyToSC : TYPE
            Position of celestrial body rel. to spacecraft in ICRF/GCRF
        SCtoEarth : TYPE
            Spacecraft position (generally from satte vector)
        GM : TYPE
            Gravity constant of celestial body

        Returns
        -------
        objectAcc : array of floats
            Acceleration caused by body expressed in ICRF/GCRF

        """

        # equations according to "Orbital Mechanics for Engineering Students", 12.11, page 713
        q = SCtoEarth.dot(2*bodyToEarth-SCtoEarth) / np.linalg.norm(bodyToEarth)**2 # 12.131
        Fq = q * (q**2 - 3*q + 3)/(1 + (1-q)**(3/2)) # F.3
        # calc object acceleration/perturbation on satellite in ECI frame
        objectAcc = GM / np.linalg.norm(bodyToSC)**3 * (Fq*bodyToEarth - SCtoEarth) # 12.130

        return objectAcc


    def get_sunVisibility(self, JDnow, GCRFpos, preloadedObjects=0):
        """
        Internal functions.
        Return the degree of Sun visibility (1=visible, 0=not visible), taking
        also into account the Penumbra (partial occultation).
        Required for precise simulation of SRP

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        GCRFpos : array of floats
            Position of spacecraft in GCRF (i.e. first half of state vector)
        preloadedObjects : int, optional
            Whether to compute the Moon pos/vel while execution (0) or use
            preloaded values (faster.)

        Returns
        -------
        v : float
            Visibility of the Sun (0-1, not visible to fully visible)
        SCtoSun : array of float
            Vector from SC to the Sun (in GCRF), needed for SRP-direction.
        occultingBodies : string
            Objects that lead to the occultation. For logfile.

        """

        # Sun pos in GCRF
        if preloadedObjects:
            SunPosGCRF = self.SunPosGCRF
        else:
            SunPos   = self.sun.at(ts.tdb_jd(JDnow)).position.m
            EarthPos = self.earth.at(ts.tdb_jd(JDnow)).position.m
            SunPosGCRF = SunPos - EarthPos

        # SC to Sun
        SCtoSun = SunPosGCRF - GCRFpos

        # get visibility for Earth
        vEarth = self.get_sunOccultation(SCtoSun, np.zeros(3), EARTH_RADIUS, GCRFpos)

        # default: no occultation by the Moon
        vMoon = 1
        # in case a moon is present, calculate also its occultation
        if 'Moon' in self.planets:
            # Moon Pos in GCRF
            if preloadedObjects:
                MoonPosGCRF = self.MoonPosGCRF
            else:
                MoonPos = self.moon.at(ts.tdb_jd(JDnow)).position.m
                MoonPosGCRF = MoonPos - EarthPos

            # get visibility for Moon
            vMoon = self.get_sunOccultation(SCtoSun, MoonPosGCRF, MOON_RADIUS, GCRFpos)

        # combine occultation from Earth and Moon
        v = np.min([vEarth, vMoon])


        # string for logfile
        if vEarth==1 and vMoon == 1:
            occultingBodies = ''
        elif vEarth<1 and vMoon ==1:
            occultingBodies = 'Earth'
        elif vEarth==1 and vMoon <1:
            occultingBodies = 'Moon'
        elif vEarth<1 and vMoon<1:
            occultingBodies = 'Earth, Moon'
        else:
            # shouldn't happen!
            occultingBodies = 'ERROR'


        return v, SCtoSun, occultingBodies

    def get_sunOccultation(self, SCtoSun, BodyPosGCRF, BodyRadius, SCpos):
        """
        Internal function.
        This function returns the amount of Sun visibility in
        regard of the provided celestial body.


        Parameters
        ----------
        SCtoSun : array of floats
            Vector from spacecraft to the Sun.
        BodyPosGCRF : array of floats
            Position of the body in GCRF.
        BodyRadius : float
            Radius of the occultation body in front of the Sun.
        SCpos : array of floats
            Position of the Spacecraft in GCRF.

        Returns
        -------
        v: float
            Sun visibility (0 <= v <= 1)

        """
        # Method by Montenbruck & Gill, Satellite Orbits, 2012

        # SC to Body
        SCtoBody = BodyPosGCRF - SCpos

        # terminate if spacecraft is below body radius (method doesnt work then)
        if np.linalg.norm(SCtoBody)<BodyRadius:
            return 0, SCtoSun

        # Method by Montenbruck & Gill, Satellite Orbits, 2012

        # angle at which Sun is seen from spacecraft
        a = np.arcsin(SUN_RADIUS/np.linalg.norm(SCtoSun))
        # angle at which Earth is seen from spacecraft
        b = np.arcsin(BodyRadius/np.linalg.norm(SCtoBody))
        # angle between the center of body and the center of the Sun
        c = np.arccos( (SCtoBody).dot(SCtoSun) / (np.linalg.norm(SCtoBody) * np.linalg.norm(SCtoSun)))


        # different cases
        if (a+b)<= c:
            v = 1 # no occultation
        elif c < (b-a) and a < b:
            v = 0 # total occulatation
        elif c < (a-b) and a > b:
            v = 1 - (np.pi * b**2)/(np.pi * a**2)
        else:
            x = (c**2 + a**2 - b**2)/(2*c)
            y = np.sqrt(a**2-x**2)
            A = a**2 * np.arccos(x/a) + b**2 * np.arccos((c-x)/b) - c * y
            v = 1 - A/(np.pi * a**2)

        return v


    def transform_GCRF_LLAgeodetic(self, JDnow, GCRFpos):
        """
        Internal function. Vectorized.
        Transform the GCRF position of a spacecraft at given Julian Date
        to its geodetic coordinates and altitude values.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        GCRFpos : array of floats
            GCRF position of spacecraft

        Returns
        -------
        Array of floats
            Latitude [°], longitude [°], geodetic altitude [m]

        """

        if GCRFpos.ndim==1:
            statepos = GCRFpos[:3]

        else:
            statepos = GCRFpos[:,:3].T

        t = ts.tdb_jd(JDnow)
        d = Distance(m=statepos)
        p = Geocentric(d.au, t=t) # GCRF
        g = wgs84.geographic_position_of(p)

        return np.squeeze([g.latitude.degrees, g.longitude.degrees, g.elevation.m]).T

    def transform_GCRF_LLAgeocentric(self, JDnow, GCRFpos):
        """
        Internal function.
        Transform the GCRF position of a spacecraft at given Julian Date
        to its geocentric coordinates and altitude values.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        GCRFpos : array of floats
            GCRF position of spacecraft

        Returns
        -------
        Array of floats
            Latitude [°], longitude [°], geocentric altitude [m]

        """

        t = ts.tdb_jd(JDnow)
        d = Distance(m=GCRFpos)
        p = Geocentric(d.au, t=t)  # GCRF
        lla = p.frame_latlon(itrs)

        return np.squeeze([lla[0].degrees, lla[1].degrees, lla[2].m])

    def transform_LLAgeodetic_GCRF(self, JDnow, lla):
        """
        Internal function.
        Transform geodetic LLA coordinates at given Julian Date to a state
        vector in the GCRF.
        See: https://rhodesmill.org/skyfield/api-topos.html

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        lla : array of floats
            Latitude [°], longitude [°], geocentric altitude [m]

        Returns
        -------
        y : array of floats
            State vector in GCRF

        """

        t = ts.tdb_jd(JDnow)
        itrsXYZ = wgs84.latlon(lla[0], lla[1], elevation_m=lla[2]).itrs_xyz
        p = ITRSPosition(itrsXYZ) # ITRS
        y = np.hstack([p.at(t).position.m, p.at(t).velocity.m_per_s]) # GCRF

        return y

    def get_EarthGravity(self, JDnow, GCRFpos, lmax, preloadedObjects=0):
        """
        Internal function.
        Returns the acceleration by the Earth including spherical harmonics.
        In- and output positions are in GCRF.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        GCRFpos : array of floats
            GCRF position of spacecraft
        lmax : int
            Max. amount of degree of spherical harmonics

        Returns
        -------
        aGCRF : array of floats
            Accleration by Earth in GCRF frame

        """

        if self.EarthTides:
            # get deltas for Cnm/Snm
            dCnm, dSnm = self.get_TidesCoeffs(JDnow, 'Earth', preloadedObjects=preloadedObjects)
            # load original coeffs
            self.clmEarth.coeffs = self.clmEarth.coeffsOrig * 1.0
            # apply deltas
            self.clmEarth.coeffs[0,:4,:4] += dCnm
            self.clmEarth.coeffs[1,:4,:4] += dSnm

        # get geocentric position of SC rel. to Earth
        LLAgeocentric = self.transform_GCRF_LLAgeocentric(JDnow, GCRFpos)

        # get spherical acceleration vector
        dUvalues = pysh.gravmag.MakeGravGridPoint(self.clmEarth.coeffs, self.clmEarth.gm, self.clmEarth.r0, LLAgeocentric[2], LLAgeocentric[0], LLAgeocentric[1], lmax)

        # calc acceleration in Earth ITRS frame
        aXYZ = self.transform_SPH_CART(LLAgeocentric[0]*DEG2RAD, LLAgeocentric[1]*DEG2RAD, dUvalues[0], -dUvalues[1],dUvalues[2])

        # rotation matrix from GCRF to ITRS
        R = itrs.rotation_at(ts.tdb_jd(JDnow))

        # apply rotation (ITRS to GCRF)
        aGCRF = np.linalg.inv(R).dot(aXYZ) # more accurate
        #aGCRF = R.T.dot(aXYZ) # faster

        return aGCRF

    def transform_GCRF_TrueEq(self, JDnow, y):
        """
        Internal function. Vectorized.
        Transforms the state vector from GCRF to "true equator and equinox
        of date", i.e. "ECI-Frame", only rotated by the GAST-angle rel. to
        the ECEF-frame.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        y : array of floats
            State vector in GCRF

        Returns
        -------
        yTrueEq : array of floats
            State vector in "true equator and equinox of date"

        """

        # get rotation matrix from GCRF to true equator and equinox of date
        # https://rhodesmill.org/skyfield/api-framelib.html#skyfield.framelib.true_equator_and_equinox_of_date
        R = true_equator_and_equinox_of_date.rotation_at(ts.tdb_jd(JDnow))

        # for single statevecs
        if y.ndim==1:
            # rotate position
            posTrueEq = R.dot(y[:3])
            # rotate velocity
            velTrueEq = R.dot(y[3:])
            # put together
            yTrueEq = np.hstack((posTrueEq, velTrueEq))
        else:
            yTrueEq = np.zeros((R.shape[2],6))
            for i in range(R.shape[2]):
                yTrueEq[i,:3] = R[:,:,i].dot(y[i,:3])
                yTrueEq[i,3:] = R[:,:,i].dot(y[i,3:])


        return yTrueEq

    def get_MoonGravity(self, JDnow, y, lmax, preloadedObjects=0):
        """
        Internal function.
        Returns the acceleration by the Moon including spherical harmonics.
        In- and output positions are in GCRF.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        y : array of floats
            State vector in GCRF
        lmax : int
            Max. amount of degree of spherical harmonics
        preloadedObjects : int, optional
            Whether to compute the Moon pos/vel while execution (0) or use
            preloaded values (faster.)

        Returns
        -------
        aXYZ_GCRF : array of floats
            Accleration by Moon in GCRF frame

        """

        if self.MoonTides:
            # get deltas for Cnm/Snm
            dCnm, dSnm = self.get_TidesCoeffs(JDnow, 'Moon', preloadedObjects=preloadedObjects)
            # load original coeffs
            self.clmMoon.coeffs = self.clmMoon.coeffsOrig * 1.0
            # apply deltas
            self.clmMoon.coeffs[0,:3,:3] += dCnm
            self.clmMoon.coeffs[1,:3,:3] += dSnm

        # get geocentric position of SC rel. to Moon (LLA and xyz, plus Moon position and rotation matrix)
        LLAgeocentric, MoonFixedPos, MoonPos, MoonVel, SatMoonXYZVel, R, _, _ = \
            self.transform_GCRF_MoonLLAplanetocentric(JDnow, y, 'PA', preloadedObjects=preloadedObjects)

        # get spherical acceleration vector
        dUvalues = pysh.gravmag.MakeGravGridPoint(self.clmMoon.coeffs, self.clmMoon.gm, self.clmMoon.r0, LLAgeocentric[2], LLAgeocentric[0], LLAgeocentric[1], lmax )

        # calc acceleration in Moon PA frame
        aXYZ = self.transform_SPH_CART(LLAgeocentric[0]*DEG2RAD, LLAgeocentric[1]*DEG2RAD, dUvalues[0], -dUvalues[1], dUvalues[2])

        # apply rotation back to GCRF
        aXYZ_GCRF = np.linalg.inv(R).dot(aXYZ)

        # calc third body acceleration
        accThirdBody = -self.MOON_GM*(MoonPos/np.linalg.norm(MoonPos)**3)

        # add third body acceleration
        aXYZ_GCRF += accThirdBody

        return aXYZ_GCRF

    
    def get_TidesCoeffs(self, JDnow, body, preloadedObjects=0):
        """
        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        body: string
            Name of celestial object for which the coefficients are calculated
            ('Earth' or 'Moon').
        preloadedObjects : int, optional
            Whether to compute the Moon pos/vel while execution (0) or use
            preloaded values (faster).

        Returns
        -------
        dCnm: array of floats
            delta values for spherical harmonics (cosine)
        dSnm: array of floats
            delta values for spherical harmonics (sine)
        """

        # GCRF position of Sun & Moon
        if preloadedObjects:
            SunPos_GCRF = self.SunPosGCRF
            MoonPos_GCRF = self.MoonPosGCRF
        else:
            # ICRF
            EarthPos_ICRF = self.earth.at(ts.tdb_jd(JDnow)).position.m
            SunPos_ICRF   = self.sun.at(ts.tdb_jd(JDnow)).position.m
            MoonPos_ICRF  = self.moon.at(ts.tdb_jd(JDnow)).position.m
            # GCRF
            SunPos_GCRF = SunPos_ICRF - EarthPos_ICRF
            MoonPos_GCRF = MoonPos_ICRF - EarthPos_ICRF
    
        # add zero speed to complete state vector 
        # (needed for transform_GCRF_MoonLLAplanetocentric()))
        SunPosVel_GCRF = np.append(SunPos_GCRF, np.zeros(3))
    
    
        if body == 'Moon':
        
            #
            # For Moon 
            #
            
            # LLA of bodies in Moon PA frame (latitude, longitude, altitude)
            object1LLA, _, _, _, _, _, _, _ = \
                self.transform_GCRF_MoonLLAplanetocentric(JDnow, np.zeros(6), preloadedObjects=preloadedObjects)
            object2LLA, _, _, _, _, _, _, _ = \
                self.transform_GCRF_MoonLLAplanetocentric(JDnow, SunPosVel_GCRF, preloadedObjects=preloadedObjects)
            
            # GM ratios
            object1GMratio = self.EARTH_GM / self.MOON_GM
            object2GMratio = self.SUN_GM / self.MOON_GM
            
            # Distance ratios
            object1Distratio = MOON_RADIUS / np.linalg.norm(MoonPos_GCRF)
            object2Distratio = MOON_RADIUS / np.linalg.norm(SunPos_GCRF-MoonPos_GCRF)
            
            # Solid Lunar Tide external Love numbers taken from GRAIL Primary Mission Data (2013)
            # https://ai-solutions.com/_help_Files/solid_tides_model.htm
            knm = np.array([[0,0,0],
                            [0,0,0],
                            [0.02408, 0.02414, 0.02394]])
            
        else:
            
            #
            # For Earth
            #
            
            # LLA 
            object1LLA = self.transform_GCRF_LLAgeocentric(JDnow, MoonPos_GCRF)
            object2LLA = self.transform_GCRF_LLAgeocentric(JDnow, SunPos_GCRF)
            
            # GM ratios
            object1GMratio = self.MOON_GM / self.EARTH_GM
            object2GMratio = self.SUN_GM / self.EARTH_GM
            
            # Distance ratios
            object1Distratio = EARTH_RADIUS / np.linalg.norm(MoonPos_GCRF)
            object2Distratio = EARTH_RADIUS / np.linalg.norm(SunPos_GCRF)
            
            # Solid Lunar Tide from IERS 2010, page 83
            knm = np.array([[0,0,0,0],
                            [0,0,0,0],
                            [0.30190,0.29830,0.30102,0],
                            [0.093,0.093,0.093,0.094]
                            ])
 
        #
        # Common processing 
        #
                
        # max degree
        nmax = knm.shape[0]-1
        
        # memory for new coeffs
        dCnm = np.zeros((nmax+1, nmax+1))
        dSnm = np.zeros((nmax+1, nmax+1))
        
        # normalized Legendre polynomials Pnm
        object1Pnm = pysh.legendre.legendre(nmax,np.sin(object1LLA[0] * DEG2RAD))
        object2Pnm = pysh.legendre.legendre(nmax,np.sin(object2LLA[0] * DEG2RAD))
        
        # start only at degree 2
        for n in range(2,nmax+1):
            for m in range(n+1):
                
                knmfactor = knm[n,m] / (2*n+1)
                
                object1 = object1GMratio * object1Distratio**(n+1) * object1Pnm[n,m]
                object2 = object2GMratio * object2Distratio**(n+1) * object2Pnm[n,m]
                           
                dCnm[n,m] = object1 * np.cos(m * object1LLA[1] * DEG2RAD) + \
                            object2 * np.cos(m * object2LLA[1] * DEG2RAD)
                            
                dSnm[n,m] = object1 * np.sin(m * object1LLA[1] * DEG2RAD) + \
                            object2 * np.sin(m * object2LLA[1] * DEG2RAD)            
                            
                dCnm[n,m] *= knmfactor            
                dSnm[n,m] *= knmfactor

        return dCnm, dSnm

    def transform_GCRF_MoonLLAplanetocentric(self, JDnow, y, refSystem='PA', preloadedObjects=0):
        """
        Internal function.
        Transforms the given GCRF position of a specraft in various
        Moon-centered positions.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        y : array of floats
            State vector in GCRF
        refSystem : str, optional
            Lunar reference system. The default is 'PA'. Possible are:
                PA: principal axis (for SH gravity calculations)
                ME: Mean Earth/Polar Axis (for LLA calculations)
        preloadedObjects : int, optional
            Whether to compute the Moon pos/vel while execution (0) or use
            preloaded values (faster.)

        Returns
        -------
        SatMoonLLA : array of floats
            Latitutde [°], longitude [°], radius to Moon center [m]
        SatMoonXYZ : array of floats
            xyz position of spacecraft in the chosen lunar reference system
        MoonPos : array of floats
            position of the Moon in GCRF
        MoonVel : array of floats
            velocity of the Moon in GCRF
        SatMoonXYZVel : array of floats
            xyz velocity of the spacecraft in the chosen lunar reference system
        R : 3x3 matrix
            rotation matrix from GCRF-frame to chosen lunar reference system
        Sat2Moon : array of floats
            Position of spacecraft rel. to Moon in GCRF
        Sat2MoonVel : array of floats
            Velocity of spacecraft rel. to Moon in GCRF

        """

        # get Moon position & velocity in ICRF/GCRF
        if preloadedObjects:
            MoonPos = self.MoonPosGCRF
            MoonVel = self.MoonVelGCRF
        else:
            MoonPos = self.moon.at(ts.tdb_jd(JDnow)).position.m - self.earth.at(ts.tdb_jd(JDnow)).position.m
            MoonVel = self.moon.at(ts.tdb_jd(JDnow)).velocity.m_per_s - self.earth.at(ts.tdb_jd(JDnow)).velocity.m_per_s


        # Satellite rel. to Moon in ICRF/GCRF
        Sat2Moon = y[:3] - MoonPos
        Sat2MoonVel = y[3:] - MoonVel

        # rotation matrix to get moon coordinates for velocity vector
        # PA (Principal Axis), as used for Gravity SH for the Moon
        if refSystem == 'PA':
            R = MoonFramePA.rotation_at(ts.tdb_jd(JDnow))
        # else, it's ME (Mean Earth/Polar Axis)
        else:
            R = MoonFrameME.rotation_at(ts.tdb_jd(JDnow))

        # Satpos in Moon XYZ
        SatMoonXYZ = R.dot(Sat2Moon)
        SatMoonXYZVel = R.dot(Sat2MoonVel)

        # Satpos in Moon LLA (geocentric, as needed for SH gravity)
        latitude  = np.arctan2(SatMoonXYZ[2],np.sqrt(SatMoonXYZ[0]**2+SatMoonXYZ[1]**2)) * RAD2DEG
        longitude = np.arctan2(SatMoonXYZ[1],SatMoonXYZ[0]) * RAD2DEG
        altitude  = np.linalg.norm(SatMoonXYZ)

        # put data together
        SatMoonLLA = np.array([latitude, longitude, altitude])

        return SatMoonLLA, SatMoonXYZ, MoonPos, MoonVel, SatMoonXYZVel, R, Sat2Moon, Sat2MoonVel


    def transform_SPH_CART(self, lat, lon, dUr, dUlat, dUlon):
        """
        Internal function.
        Transforms a 3d vector provided in spherical coordinates at given
        coordinates into a 3d cartesian vector.
        For geocentric coordinates only!

        Parameter
        ----------
        lat : float
            latiitude [°]
        lon : float
            longitude [°]
        dUr : float
            vector component in radial direction
        dUlat : float
            vector component along the meridians ("north/south")
        dUlon : float
            vector component along the parallels ("east/west")

        Returns
        -------
        aXYZ : array of floats
            3d vector in cartesian coordinates

        """

        clat = np.cos(lat)
        clon = np.cos(lon)
        slat = np.sin(lat)
        slon = np.sin(lon)

        rot_matrix = np.array([
            [clat*clon, -slon, -slat*clon],
            [clat*slon,  clon, -slat*slon],
            [slat,       0.   , clat]
            ])

        aXYZ = rot_matrix.dot(np.array([dUr, dUlon, dUlat]))

        return aXYZ

    def get_dynPress(self, MET, y, JDnow):
        """
        Internal function.
        Returns the dynamic pressure in Pascals [kg/(m*s^2)]

        Parameters
        ----------
        MET : float
            Mission Elapsed Time
        y : array of floats
            State vector
        JDnow : float
            Current Julian Date (TBD)

        Returns
        -------
        dynPress : array of floats
            dynamic pressure as a vector (in direction of velocity vector)

        """

        # get temperature/pressure/density/M1
        atmosvalues = self.get_atmos(JDnow, y[:3])

        # get relative velocity
        ECEFvel = self.get_relVelocityEarth(JDnow, y)

        # calc dynamic pressure (as a vector)
        # the linear transformation can be used to get density values similar to MSISE-90
        if self.transToMSISE90:
            dynPress = 0.5 * (self.RhoGain * atmosvalues[2] + self.RhoOffset) * np.linalg.norm(ECEFvel) * ECEFvel
        else:
            dynPress = 0.5 * atmosvalues[2] * np.linalg.norm(ECEFvel) * ECEFvel

        return dynPress

    def get_relVelocityEarth(self, JDnow, y):
        """
        Internal function.
        Returns the velocity vector relative to the (rotating) Earth. Needed
        for atmospheric influence, e.g. for drag calculations.
        Not vectorized for fast computation of single values

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        y : array of floats
            State vector in GCRF

        Returns
        -------
        ECEFvel : array of floats
            Relative velocity vector, described in GCRF

        """

        # Earth rotation axis in ITRS
        earthRotAxisITRS = np.array([0,0,1])
        #  GCRF to ITRS transformation
        R = itrs.rotation_at(ts.tdb_jd(JDnow))
        # Earth rotation axis in GCRF
        earthRotAxisGCRF = R.T.dot(earthRotAxisITRS)
        # relative velocity
        ECEFvel = y[3:] - np.cross(EARTH_ROT_SPEED * earthRotAxisGCRF, y[:3])

        return ECEFvel
    
    def get_liftVector(self, gVec, vrel):
        """
        Returns the vector of the lift direction.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF.
        gVec : array of floats
            Guidance vector / thrust direction / rocket pointing direction.
        vrel : array of floats
            Relative velocity of spacecraft w.r.t Earth surface/atmosphere.

        Returns
        -------
        liftVec: array of floats
            Direction of lift.

        """
        
        # normal to plane described by vrel and gVec
        gVecXvrel = np.cross(gVec, vrel)
        
        # lift direction
        liftVec = np.cross(vrel, gVecXvrel)
        liftVec = liftVec / np.linalg.norm(liftVec)
        
        return liftVec
    
    def get_AoAangle(self, vrel, gVec):
        """
        Returns Angle of Attack (AoA). Vectorized.

        Parameters
        ----------
        vrel : array of floats
            Velocity vector w.r.t. to Earth surface (relative velocity).
        gVec : array of floats
            Guidance vector / thrust direction / rocket pointing direction.

        Returns
        -------
        AoA : float / array of floats
            Angle between velocity and rocket pointing direction [rad].

        """
        
        # single values
        if vrel.shape == (3,):
            cosAoA = vrel.dot(gVec)/(np.linalg.norm(vrel)*np.linalg.norm(gVec))
            
            # sometimes, cos values are out of -1/1; correction needed
            if cosAoA>1:
                cosAoA = 1
            elif cosAoA<1:
                cosAoA = -1
            
            AoA = np.arccos(cosAoA)
        
        # multiple values
        else:
            cosAoA = np.diag(vrel.dot(gVec.T))/\
                     (np.linalg.norm(vrel, axis=1)*np.linalg.norm(gVec, axis=1))
            
            # sometimes, cos values are out of -1/1; correction needed
            cosAoA[cosAoA > 1.] =  1.
            cosAoA[cosAoA <-1.] = -1.
            
            AoA = np.arccos(cosAoA)
    
        return AoA

    def get_atmos(self, JDnow, y):
        """
        Internal function.
        Returns relevant atmospheric properties.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        y : array of floats
            State vector in GCRF

        Returns
        -------
        T : float
            Temperature [K]
        Pa : float
            Pressure [Pa]
        Rho : float
            Density [kg/m^3]
        M1 : float
            Mach 1 number [1]

        """

        # get geodetic coordinates (lat [°], lon [°], alt [m]), as needed for NRL MSISE-00
        LLAgeodetic = self.transform_GCRF_LLAgeodetic(JDnow, y[:3])

        # fix altitude if negative, otherwise MSISE model returns bad values (neg. temp...)
        if LLAgeodetic[2] < 0:
            LLAgeodetic[2] = 0

        # chose in function of atmospheric function what to do
        if self.atmosModel == 'nrlmsise00':

            # call nrlmsise00
            atmosTRho = msise_model(ts.tdb_jd(JDnow).utc_datetime(),
                                    LLAgeodetic[2]/1000,
                                    LLAgeodetic[0],
                                    LLAgeodetic[1],
                                    self.f107s , self.f107, self.Ap)
            # extract values
            T = atmosTRho[1][1]
            Rho = atmosTRho[0][5] * 1000
            # calculate pressure (T * Rho * R)
            Pa = T * Rho * R
            # calc sound of speed, i.e. Mach 1
            # https://en.wikipedia.org/wiki/Speed_of_sound
            M1 = 20.05 * np.sqrt(T)

        # no valid modell
        else:
            # temperature of Cosmic Microwave Background Radiation
            return 2.73, 0, 0, 33


        return T, Pa, Rho, M1

    def get_OrbElements(self, JDnow, y, GM):
        """
        Internal function. Vectorized.
        Returns the orbital elements of a spacecraft in the orbit of a planet
        https://rhodesmill.org/skyfield/elements.html

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        y : array of floats
            State vector in the reference frame of planet
        GM : float
            gravity constant of planet

        Returns
        -------
        orbElements : array of floats
            Various elements, see above URL.

        """

        if y.ndim==1:
            statepos = y[:3]
            statevel = y[3:]
        else:
            statepos = y[:,:3].T
            statevel = y[:,3:].T

        # make Skyfield position and velocity objects
        # https://rhodesmill.org/skyfield/api-units.html#skyfield.units.Distance
        position = Distance(m=statepos)
        velocity = Velocity(km_per_s = statevel/1000)

        # make Skyfield time series
        JDvec = ts.tdb_jd(JDnow)

        # get orbital elements
        orbElements = OsculatingElements(position, velocity, JDvec, GM/1000**3)

        return orbElements



    def get_ENUvec(self, JDnow, y):
        """
        Internal function. Vectorized.
        Calculates local ENU vectors for provided statevec's positions in an
        object's frame.
        Following use of these ENU vectors require data processing in object's
        frame (e.g. for HA/FPA)
        Assumes round object.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        y : array of floats
            State vector in the reference frame of planet


        Returns
        -------
        ENU : 3x3 array of floats
            East, North, Up vectors / rotation matrices in local ref. system

        """

        # check if single statevec or if list of statevecs
        if y.shape == (6,) or y.shape == (3,):
            y = np.expand_dims(y, axis=0)

        # get UP vectors
        UP = y[:,:3] * 1.0
        UP /= np.linalg.norm(UP, axis=1)[:,None]

        # get EAST vectors
        EAST = np.vstack([-y[:,1],y[:,0],np.zeros(len(y))]).T
        EAST /= np.linalg.norm(EAST, axis=1)[:,None]

        # get NORTH vectors
        NORTH = np.cross(UP, EAST)
        #NORTH /= np.linalg.norm(NORTH, axis=1)[:,None] # finetune length to 1; skipping to save time

        ENU = np.stack((EAST, NORTH, UP), axis=2) # E, N, U vecs are vertical, horizontal stacking

        return ENU


    def get_ENUvec_Earth(self, JDnow, y, frame='TOD'):
        """
        Internal function. Vectorized.
        Calculates ENU vectors of provided statevec positions for Earth, in
        GCRF (and not ITRS!)
        Following use of these ENU vectors require data processing to be done
        in GCRF (e.g. for HA/FPA)

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        y : array of floats
            State vector in the reference frame of planet
        frame: string
            Frame (TOD or WGS84) in which ENU is determined. WGS84 considers
            the shape of Earth, TOD assumes circular object. 
            Default is TOD. WGS84 only usueful for initial launch direction.


        Returns
        -------
        ENU : 3x3 array of floats
            East, North, Up vectors / rotation matrices in GCRF

        """
        # check if single statevec or if list of statevecs
        if y.shape == (6,):
            y = np.expand_dims(y, axis=0)
            
        if frame=='WGS84':

            t = ts.tdb_jd(JDnow)
            d = Distance(m=y[:,:3].T) # input is nxm, n = 3, m = number of provided statevecs
            p = Geocentric(d.au, t=t) # p is in GCRF
            g = wgs84.geographic_position_of(p) # WGS84 position
    
            NEU = g.rotation_at(t)
            # NEU returns a 3D matrix:
            # axis 0: the three components of the originating system (North, East, Up)
            # axis 1: the m different times (=amount of statevecs provided)
            # axis 2: the n (three) dimensions of the target systemt
    
            # New matrix to bring into right order
            ENU = np.zeros((NEU.shape))
            ENU[0,:,:] = NEU[1,:,:] # first row becomes East
            ENU[1,:,:] = NEU[0,:,:] # seconds row becomes North
            ENU[2,:,:] = NEU[2,:,:] # last row remains Up
    
            # for return value, the matrix is transposed
            # axis 0: the m different times (=amount of statevecs provided)
            # axis 1: thee three dimensions of the target system (GCRF)
            # axis 2: the three dimensions of the originating system (East, North, Up)
            
            ENU = ENU.T

        # true equator and equinox of date (TETE, TOD)
        else:
            R = true_equator_and_equinox_of_date.rotation_at(ts.tdb_jd(JDnow))

            # position in TOD
            pos_TOD = np.einsum('ij...,j...->i...', R, y[:,:3].T).T 
             
            # get local ENU in TOD
            ENU_TOD = self.get_ENUvec(JDnow, pos_TOD)
              
            # transform back to GCRF - not vectorized; #TODO
            if R.shape == (3,3):
                ENU = R.T.dot(np.squeeze(ENU_TOD))
            else:
                ENU = np.zeros((len(JDnow),3,3))
                for i in range(len(JDnow)):
                    ENU[i,:,:] = R[:,:,i].T.dot(ENU_TOD[i,:,:])

        return np.squeeze(ENU)

    def get_ENUvec_EarthGravity(self, JDnow, y, lmax):
        """
        Internal function. NOT vectorized
        Calculates ENU vectors at provided statevec position for Earth in GCRF,
        with UP = direction of gravity.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD)
        y : array of floats
            State vector in the reference frame of planet
        lmax : int
            Max. amount of degree of spherical harmonics

        Returns
        -------
        ENU_GCRF : 3x3 array of floats
            Earth, North, Up vectors / rotation matrices in GCRF

        """
        # get geocentric position of SC rel. to Earth
        LLAgeocentric = self.transform_GCRF_LLAgeocentric(JDnow, y[:3])

        # get spherical acceleration vector
        dUvalues = pysh.gravmag.MakeGravGridPoint(self.clmEarth.coeffs, self.clmEarth.gm, self.clmEarth.r0, LLAgeocentric[2], LLAgeocentric[0], LLAgeocentric[1], lmax)

        # calc acceleration in Earth ITRS frame
        aXYZ = self.transform_SPH_CART(LLAgeocentric[0]*DEG2RAD, LLAgeocentric[1]*DEG2RAD, dUvalues[0], -dUvalues[1],dUvalues[2])

        # the UP vector is in the opposite direction of the acceleration vector
        UP = -aXYZ * 1.0
        UP /= np.linalg.norm(UP)

        # get EAST vector
        EAST = np.array([-UP[1], UP[0], 0])
        EAST /= np.linalg.norm(EAST)

        # get NORTH vector
        NORTH = np.cross(UP, EAST)
        NORTH /= np.linalg.norm(NORTH)

        # compbine vectors
        ENU = np.vstack((EAST, NORTH, UP)).T

        # rotation matrix from GCRF to ITRS
        R = itrs.rotation_at(ts.tdb_jd(JDnow))

        # apply rotation (ITRS to GCRF)
        ENU_GCRF = np.linalg.inv(R).dot(ENU)

        return ENU_GCRF

    def get_LVLHframe(self, JDnow, y, planet='Earth'):
        """
        Internal function.
        Get local LVLH frame. Not vectorized, single statevecs only.
        https://ai-solutions.com/_freeflyeruniversityguide/attitude_reference_frames.htm


        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF
        planet : string, optional
            Name of central body. The default is: 'Earth'.

        Returns
        -------
        LVLH : 3x3 array of floats
            LVLH rotation matrix (LVLH to GCRF)


        """

        # get pos/vel in function of central object
        if planet=='Earth':
            pos = y[:3]
            vel = y[3:]
        elif planet=='Moon':
            _, _, _, _, _, _, pos, vel = self.transform_GCRF_MoonLLAplanetocentric(JDnow, y)

        # z value (towards planet)
        LVLH_z = -pos/np.linalg.norm(pos)

        # y value (normal to orbital plane)
        LVLH_y = np.cross(LVLH_z, vel)
        LVLH_y /= np.linalg.norm(LVLH_y)

        # x value (y x z)
        LVLH_x = np.cross(LVLH_y, LVLH_z)
        LVLH_x /= np.linalg.norm(LVLH_x)

        # assemble rotation matrix
        LVLH = np.vstack((LVLH_x,LVLH_y,LVLH_z)).T

        return LVLH


    def apply_deltaV_LVLH(self, JDnow, y, deltaV, planet='Earth'):
        """
        Internal function.
        Adds a delta-v (expressed in LVLH) to the velocity vector, itself
        returned in GCRF system.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF.
        deltaV : array of floads
            delta V components provided in LVLH.
        planet : string, optional
            Name of central body. The default is 'Earth'.

        Returns
        -------
        Array of flaots
            Updated state vector in GCRF

        """

        # get local LVLH frame
        LVLH = self.get_LVLHframe(JDnow, y, planet)

        # get velocity in LVLH frame and apply deltaV (described in LVLH)
        velLVLH = LVLH.T.dot(y[3:]) + deltaV

        # transform back into original frame and return statevec
        return np.hstack([y[:3], LVLH.dot(velLVLH)])


    def get_VNBframe(self, JDnow, y, planet='Earth'):
        """
        Internal function.
        Get local VNB frame. Not vectorized, single statevecs only.
        https://ai-solutions.com/_freeflyeruniversityguide/attitude_reference_frames.htm
        Works for:
            - GCRF (Earth)
            - Earth-fixed
            - Moon-centered

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF
        planet : string, optional
            Name of central body. The default is 'Earth'. Options are:
                'Earth' --> equals to GCRF
                'EarthFixed' --> w.r.t. to Earth surface
                'Moon' --> w.r.t. to the Moon

        Returns
        -------
        VNB : 3x3 array of floats
            VNB rotation matrix (VNB to GCRF).

        """

        # get pos/vel in function of central object
        if planet=='Earth':
            pos = y[:3]
            vel = y[3:]
        elif planet=='EarthFixed':
            # get GCRF position
            pos = y[:3]

            # relative velocity
            vel = self.get_relVelocityEarth(JDnow, y)

        elif planet=='Moon':
            _, _, _, _, _, _, pos, vel = self.transform_GCRF_MoonLLAplanetocentric(JDnow, y)


        # x value (equals velocity vector direction)
        VNB_x = vel/np.linalg.norm(vel)

        # y value (normal to orbital plane)
        VNB_y = np.cross(pos, vel)
        VNB_y /= np.linalg.norm(VNB_y)

        # z value (x x y)
        VNB_z = np.cross(VNB_x, VNB_y)
        VNB_z /= np.linalg.norm(VNB_z)

        # assemble rotation matrix
        VNB = np.vstack((VNB_x,VNB_y,VNB_z)).T

        return VNB

    def apply_deltaV_VNB(self, JDnow, y, deltaV, planet='Earth'):
        """
        Internal function.
        Adds a delta-v (expressed in VNB) to the velocity vector, itself
        returned in GCRF system.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF.
        deltaV : array of floads
            delta V components provided in VNB.
        planet : string, optional
            Name of central body. The default is 'Earth'. Options are:
                'Earth' --> equals to GCRF
                'EarthFixed' --> w.r.t. to Earth surface
                'Moon' --> w.r.t. to the Moon


        Returns
        -------
        Array of flaots
            Updated state vector in GCRF

        """

        # get local VNB frame
        VNB = self.get_VNBframe(JDnow, y, planet)

        # get velocity in VNB frame and apply deltaV (described in LVLH)
        velVNB= VNB.T.dot(y[3:]) + deltaV

        # transform back into original frame and return statevec
        return np.hstack([y[:3], VNB.dot(velVNB)])

    def get_FPAHAVEL_EF(self, JDnow, y):
        """
        Internal function. Vectorized.
        Returns velocity vector described w.r.t. to fixed Earth frame (EF).

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF.

        Returns
        -------
        FPA : float
            Flight Path Angle [rad] w.r.t. to Earth surface.
        HA : float
            Heading Angle [rad] w.r.t. to Earth surface.
        EFVEL : float
            Earth-fixed velocity [m/s].

        """
        
        # check if single statevec or if list of statevecs
        if y.shape == (6,):
            y = np.expand_dims(y, axis=0)

        # change velocity vectors from inertial to Earth-fixed

        # Earth rotation axis in ITRS
        earthRotAxisITRS = np.array([0,0,1])
        #  GCRF to ITRS
        R = itrs.rotation_at(ts.tdb_jd(JDnow))
        
        # add dimension to R if only a float value (instead of array) of JD
        # values was provided
        if R.shape == (3,3):
            R = np.expand_dims(R, axis=2)
            
        R = np.transpose(R,(2,0,1)) # first dimension is time

        # Earth rotation axis in GCRF
        earthRotAxisGCRF = np.transpose(R,(0,2,1)).dot(earthRotAxisITRS)

        yEF = y * 1.0
        yEF[:,3:] -= np.cross(EARTH_ROT_SPEED * earthRotAxisGCRF, yEF[:,:3])

        # get FPA / HA
        FPA, HA = self.get_FPAHA(JDnow, yEF, frame='Earth')

        # EF vel
        EFVEL = np.linalg.norm(yEF[:,3:], axis=1)

        return np.squeeze(FPA), np.squeeze(HA), np.squeeze(EFVEL)

    def get_FPAHA(self, JDnow, y, frame='Earth'):
        """
        Internal function. Vectorized.
        Returns FPA/HA in frame of provided statevectors. In case planet Earth
        is used, FPA/HA are calculated in the true local ENU frame in ITRS.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TDB).
        y : array of floats
            State vector in local frame (GCRF for Earth)
        frame : string, optional
            Frame in which FPA/HA are calculated. The default is 'Earth' 

        Returns
        -------
        FPA : float
            Flight Path Angle [rad] w.r.t. to state vector system / Earth surface
        HA : float
            Heading Angle [rad] w.r.t. to state vector system / Earth surface

        """

        # check if single statevec or if list of statevecs
        if y.shape == (6,):
            y = np.expand_dims(y, axis=0)

        if frame == 'Earth':
            # get local ENU values for Earth
            ENU = self.get_ENUvec_Earth(JDnow, y)
        elif frame == 'launchsite':
            ENU = self.ENU_liftoff * 1.0
        else:
            ENU = self.get_ENUvec(JDnow, y)
            
        # add dimension to ENU if only a 2D matrix was returned
        if ENU.shape == (3,3):
            ENU = np.expand_dims(ENU, axis=0)

        # calc angle between velocity vec and local UP
        # - calc dot product between UP and velocity
        # - divise by product of their norms to get costheta
        # - adjust to 1/-1 if necessary
        # - calc angle

        # calc angle using the formula for vector dot multiplication (A.B = cos(theta)*norm(A)*norm(B))
        # ENU[:,:,2] = UP vectors (last entry in highest dimension, i.e. columns at the very right at all times)
        UPdotV = np.diag(ENU[:,:,2].dot(y[:,3:].T))
        costheta = UPdotV / (np.linalg.norm(ENU[:,:,2], axis=1) * np.linalg.norm(y[:,3:], axis=1))

        # sometimes, cos values are out of -1/1; correction needed
        costheta[costheta > 1.] =  1.
        costheta[costheta <-1.] = -1.

        # calc flight path angle [rad]
        FPA = np.arcsin(costheta)

        # projected vector
        vProjected = y[:,3:] - np.diag(y[:,3:].dot(ENU[:,:,2].T))[:,None] * ENU[:,:,2]

        # save norm of projected vector, needed afterwards
        vProjectedNorm = np.linalg.norm(vProjected,axis=1)

        # find where no projected vector exists
        ind0vProjected = vProjectedNorm==0
        # correct when there is no projected vector (set it to NORTH vector)
        vProjected[ind0vProjected] = ENU[ind0vProjected,:,1]

        # calc cos of vProjected w.r.t. NORTH
        costheta_north = np.diag(ENU[:,:,1].dot(vProjected.T))[:,None] / (np.linalg.norm(ENU[:,:,1],axis=1)*vProjectedNorm)[:,None]
        # calc cos of vProjected w.r.t. EAST
        costheta_east  = np.diag(ENU[:,:,0].dot(vProjected.T))[:,None] / (np.linalg.norm(ENU[:,:,0],axis=1)*vProjectedNorm)[:,None]

        # calculate azimuth angle; 0 angle is the north, 90° is east [rad]
        HA = np.arctan2(costheta_east, costheta_north)

        # find negative values
        indNegAzimuth = HA<0
        # correct to positive values
        HA[indNegAzimuth] += 2*np.pi

        return np.squeeze(FPA), np.squeeze(HA)


    def get_RPSpos(self, JDnow, y, updateMissionDF=0):
        """
        Internal function. Vectorized.
        Gets the rotating-pulsating-system coordinates of the spacecraft w.r.t.
        Earth and Moon. Inlcudes possibility to get Moon position and velocity
        by function or call or by using data provided in mission dataframe.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF
        updateMissionDF : int, optional
            If set to 1, Moon position/velocity are taken from mission dataframe.
            The default is 0.

        Returns
        -------
        RPScoord : array of floats
            RPS position of spacecraft

        """

        # get MoonFrame information if not using MissionDF; assuming JDnow is a single float
        if not updateMissionDF:
            SatMoonLLA, SatMoonXYZ, MoonPos, MoonVel, SatMoonXYZVel, R, _, _ = self.transform_GCRF_MoonLLAplanetocentric(JDnow, y)
            # expand dimensions, needed for vectorization of later calculations
            MoonPos = np.expand_dims(MoonPos, axis=0)
            MoonVel = np.expand_dims(MoonVel, axis=0)
            y = np.expand_dims(y, axis=0)
        # otherwise, extract from dataframe
        else:
            MoonPos = self.missionDF[['MoonPosX', 'MoonPosY', 'MoonPosZ']].to_numpy()
            MoonVel = self.missionDF[['MoonVelX', 'MoonVelY', 'MoonVelZ']].to_numpy()

        # distance Earth-Moon
        MoonDist = np.linalg.norm(MoonPos, axis=1)

        # length (norm) of velocity of Moon
        MoonVelValue = np.linalg.norm(MoonVel, axis=1)

        # normalized spacecraft position w.r.t. Earth-Moon distance
        r_SC = y[:,:3] / MoonDist[:,None]
        # normalized Moon position w.r.t. Earth-Moon distance (r vector)
        MoonPosNorm = MoonPos / MoonDist[:,None]
        # normalized Moon velocity w.r.t.
        MoonVelNorm = MoonVel / MoonVelValue[:,None]

        # normal to ecliptic plane
        h = np.cross(MoonPosNorm, MoonVelNorm)

        # velocity vector normal to r/h-plane
        v = np.cross(h, MoonPosNorm)
        # get its length (norm)
        vValue = np.linalg.norm(v, axis=1)
        # normalize velocity vector normal to r/h-plane
        vNorm = v / vValue[:,None]

        # rotation matrix
        R = np.stack((MoonPosNorm, vNorm, h), axis=2) # vertical stacking

        # memory for new coordinates
        RPScoord = np.zeros((R.shape[0],3))

        # rotate SC position (normalized w.r.t. to Earth-Moon distance) into RPS frame
        for i in range(R.shape[0]):
            RPScoord[i,:] = R[i,:,:].T.dot(r_SC[i,:])

        return RPScoord
    
    def get_EarthRangeToLaunchsite(self, JDnow, y):
        """
        Internal function; not vectorized.
        Returns the ground distance between spacecraft and actual position of
        launchsite.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF.

        Returns
        -------
        rangeToLaunchsite : float
            Ground distance between spacecraft and launchsite.

        """
        
        # current GCRF position of launchsite
        y_launchsite = self.transform_LLAgeodetic_GCRF(JDnow, self.MD.launchsite_LLA)
        
        # calc angle between launchsite position and spacecraft position
        angle_between_positions = np.arccos(y[:3].dot(y_launchsite[:3])/\
                                            (np.linalg.norm(y[:3])*np.linalg.norm(y_launchsite[:3])))
            
        # get distance on ground (assuming round Earth)
        rangeToLaunchsite = angle_between_positions * EARTH_RADIUS
        
        return rangeToLaunchsite

    def transform_J2000SV(self, JDnow, y, targetFrame='GCRF'):
        """
        Internal function.
        Transforms MRS state vector from GCRF to EME2000/FK5.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF or J2000/FK5
        targetFrame: string
            Name of frame to transform to: GCRF (J2000 SV) or J2000 (GCRF SV)

        Returns
        -------
        yTargetSystem : array of floats
            State vectors in required target frame

        """


        # The IAU Resolutions on Astronomical Reference Systems, Time Scales,
        # and Earth Rotation Models, 2005, page 28
        J2000toGCRF = np.array([
            [ 1.00000000000000, 0.00000007078280, -0.00000008056149],
            [-0.00000007078280, 1.00000000000000, -0.00000003306041],
            [ 0.00000008056149, 0.00000003306041,  1.00000000000000]
            ])

        # check if single statevec or if list of statevecs
        if y.shape == (6,):
            y = np.expand_dims(y, axis=0)

        # when going from GCRF to J2000
        if targetFrame=='GCRF':
            newPos = J2000toGCRF.dot(y[:,:3].T).T
            newVel = J2000toGCRF.dot(y[:,3:].T).T
        # going from J2000 to GCRF
        elif targetFrame=='J2000':
            newPos = J2000toGCRF.T.dot(y[:,:3].T).T
            newVel = J2000toGCRF.T.dot(y[:,3:].T).T
        # catch error
        else:
            print('MRS:\t\tERROR transform_J2000SV(): unknown target frame ', targetFrame)
            newPos = y[:,:3]
            newVel = y[:,3:]

        # put position and velocity into one matrix
        yTargetSystem = np.append(newPos, newVel, axis=1)

        return yTargetSystem


    def make_TempDF(self, nrows):
        """
        Internal function.
        Makes a temporary dataframe to be used for propagation results.

        Parameters
        ----------
        nrows : int
            Number of rows in temporary dataframe.

        Returns
        -------
        None.

        """

        # preset column names
        colNames = [
            'JD_TBD',
            'MET',
            'segmentID','segmentType','configID',
            'propaMode', 'forcesID','atmosModel',
            'activeSC',
            'x','y','z',
            'vx','vy','vz',
            ]

        self.TempDF = pd.DataFrame(index = range(nrows), columns=colNames)

        # set dtype, because otherwise it's object
        self.TempDF[['JD_TBD','x','y','z','vx','vy','vz']] = self.TempDF[['JD_TBD', 'x','y','z','vx','vy','vz']].astype('float')


        return None

    def add_event(self, JDnow, y, eventType):
        """
        Internal function.
        Appends event statevec and name to eventsDF.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF.
        eventType : str
            Event type name.

        Returns
        -------
        None.

        """

        # pointer to new row in dataframe
        DFpointer = len(self.eventsDF)

        # add values to events dataframe
        self.eventsDF.loc[DFpointer, 'JD_TBD'] = JDnow
        self.eventsDF.loc[DFpointer, [ 'x','y','z','vx','vy','vz']] = y.tolist()
        self.eventsDF.loc[DFpointer, 'eventType'] = eventType

        return None


    def expand_DFname(self, typelist, DFname='missionDF'):
        """
        Adds calculated values for the given list of data types to the
        designated dataframe. Takes name of DF als argument.

        Parameters
        ----------
        typelist : array of strings
            List of data types to be added to the an existing dataframe
        DF : str. The default is 'missionDF'.
            Dataframe to be expanded.

        Returns
        -------
        None.

        """

        # different proceedures for different dataframes
        if DFname == 'missionDF':
            self.missionDF = self.expand_DF(typelist, self.missionDF)
        elif DFname == 'eventsDF':
            self.eventsDF = self.expand_DF(typelist, self.eventsDF)
        else:
            print('MRS:\t\tERROR expand_DF_new(): unknown kind of dataframe ', DFname)

        return None

    def expand_DF(self, typelist, DF):
        """
        Internal function.
        Adds calculated values for the given list of data types to the
        designated dataframe. Takes dataframe itself as argument

        Parameters
        ----------
        typelist : array of strings
            List of data types to be added to the an existing dataframe
        DF : pandas dataframe.
            Dataframe to be expanded.

        Returns
        -------
        DF : pandas dataframe
            Dataframe with some new values

        """

        # setup
        lenDF = len(DF)
        JDs = np.squeeze(DF[['JD_TBD']].to_numpy())
        METs = np.squeeze(DF[['MET']].to_numpy())
        statevecs = DF[['x','y','z','vx','vy','vz']].to_numpy()

        # loop through provided list of data types
        for datatype in typelist:

            # adding Earth LLA
            if datatype == 'EarthLLA':
                print('MRS:\t\tAdding Earth LLA to dataframe.')
                DF[['EarthLat', 'EarthLon', 'EarthAlt']] = self.transform_GCRF_LLAgeodetic(JDs, statevecs[:,:3])

            # adding atmospheric properties
            elif datatype == 'EarthAtmos':
                print('MRS:\t\tAdding Earth atmospheric values to dataframe.')

                # initial segment pointer
                segmentID = -1

                # preload space weather
                self.swDF = sw.sw_daily()

                # get step wise atmospheric values
                for i in range(lenDF):

                    # detect new segment and update atmospheric values
                    if DF.segmentID[i]>segmentID:
                        segmentID = DF.segmentID[i]

                        # load parameters for given atmosphere
                        self.atmosModel= DF.atmosModel[i]

                        # check if spaceweather needs to be loaded
                        if self.atmosModel == 'nrlmsise00':
                            # use F107 values for time at mission start (dirty, different values than in simulation)
                            datestring = ts.tdb_jd(JDs[i]).utc_strftime(format='%Y-%m-%d')

                            # save values if space weather in use
                            if self.use_spaceweather:
                                indices = self.swDF.loc[datestring]
                                self.f107s = indices['f107_81lst_adj']
                                self.f107  = indices['f107_adj']
                                self.Ap    = indices['Apavg']

                    # get and add atmospheric values to dataframe
                    DF.loc[i,['atmosT', 'atmosPa', 'atmosRho', 'atmosM1']] = \
                        self.get_atmos(JDs[i], statevecs[i])


            # Earth acceleration
            elif datatype == 'EarthAcceleration':
                print('MRS:\t\tAdding Earth gravity (incl. SH/tides if selected) to dataframe.')
                
                # initial segment pointer
                self.forcesID = -1

                # get step wise atmospheric values
                for i in range(lenDF):

                    # detect new forces config 
                    if DF.forcesID[i]!=self.forcesID:
                        self.forcesID = DF.forcesID[i]
                        
                        # set up tides
                        self.EarthTides = self.MD.forcesSettings.EarthTides[self.forcesID]

                        # get lmax
                        lmax = self.MD.forcesSettings.EarthSHn[self.forcesID]

                    # call function if SH number is provided
                    if lmax:
                        objectAcc = self.get_EarthGravity(JDs[i], statevecs[i,:3], lmax)
                    # if no lmax, use point gravity
                    else:
                        r_norm = np.linalg.norm(statevecs[i,:3])
                        objectAcc = - self.EARTH_GM * statevecs[i,:3]/r_norm**3


                    DF.loc[i,['EarthAccX','EarthAccY','EarthAccZ']] = objectAcc



            # EF velocity/HA/FPA
            elif datatype == 'EarthFixedFPAHAvel':
                print('MRS:\t\tAdding FPA/HA/VEL (w.r.t. to Earth surface) to dataframe.')
                FPA, HA, EFVEL = self.get_FPAHAVEL_EF(JDs, statevecs)
                DF['EarthFixedFPA'] = FPA * RAD2DEG
                DF['EarthFixedHA']  = HA * RAD2DEG
                DF['EarthFixedVEL'] = EFVEL

            # dyn. pressure
            elif datatype == 'dynPress':

                if not 'EarthFixedVEL' in DF.columns:
                    DF = self.expand_DF(['EarthFixedFPAHAvel'], DF)

                print('MRS:\t\tAdding Earth atmospheric dynamic pressure to dataframe.')
                DF['dynPress'] = .5 * DF['EarthFixedVEL']**2 * DF['atmosRho']

            # mach
            elif datatype == 'Mach':
                print('MRS:\t\tAdding Mach number to dataframe.')
                if not 'EarthFixedVEL' in DF.columns:
                    print('MRS:\t\tERROR: EarthFixedFPAHAvel needs to be added first; skipping dyn. pressure/Mach.')
                    continue
                if not 'atmosM1' in DF.columns:
                    print('MRS:\t\tERROR: EarthAtmos needs to be added first; skipping Mach.')
                    continue
                DF['Mach'] = DF['EarthFixedVEL'].to_numpy() / DF['atmosM1'].to_numpy()

            # angle of attack
            elif datatype == 'AoA':
               
                if not 'gVecX' in DF.columns:
                    DF = self.expand_DF(['GuidanceVec'], DF)
                
                print('MRS:\t\tAdding Angle of Attack to dataframe.')
                
                vrel = np.zeros((lenDF,3))
                for i in range(lenDF):
                    vrel[i,:] = self.get_relVelocityEarth(JDs[i], statevecs[i])
                
                DF['AoA'] = self.get_AoAangle(vrel, DF[['gVecX', 'gVecY', 'gVecZ']].to_numpy()) * RAD2DEG
                


            # static spacecraft properties
            elif datatype == 'SpacecraftStatic':

                if self.missionDFtype == 2:
                    print('MRS:\t\tERROR: External mission has no spacecraft; skipping Spacecraft data.')
                    continue

                print('MRS:\t\tAdding static spacecraft values to dataframe.')

                staticMass, staticDragArea, staticCd, staticCr, staticSRParea = self.SC.get_staticvalues()
                DF['SC_static_AreaCd'] = staticDragArea * staticCd
                DF['SC_static_Mass'] = staticMass
                DF['SC_static_Cr'] = staticCr
                DF['SC_static_SRParea'] = staticSRParea



            # static spacecraft properties
            elif datatype == 'SpacecraftActive':

                if self.missionDFtype == 2:
                    print('MRS:\t\tERROR: External mission has no spacecraft; skipping Spacecraft data.')
                    continue

                if self.SC.mode == 'static':
                    print('MRS:\t\tERROR: mission has only static spacecarft; skipping active spacecraft data.')
                    #continue

                print('MRS:\t\tAdding active spacecraft values to dataframe.')

                for i in range(lenDF):
                    SCthrust, SCmass, returnVals = \
                        self.SC.get_ThrustMassOther(DF.MET[i], pa=DF.atmosPa[i], verbose='')
                    DF.loc[i,['SC_active_Thrust', 'SC_active_Mass']] = SCthrust, SCmass

            # atmospheric drag
            elif datatype == 'DragForce':

                if not 'atmosRho' in DF.columns:
                    DF = self.expand_DF(['EarthAtmos'], DF)

                if not 'dynPress' in DF.columns:
                    DF = self.expand_DF(['dynPress'], DF)

                if not 'SC_static_AreaCd' in DF.columns:
                    DF = self.expand_DF(['SpacecraftStatic'], DF)

                print('MRS:\t\tAdding drag force to dataframe.')
                for i in range(lenDF):

                    # check if active or static SC
                    if DF.activeSC[i]:
                        DF.loc[i,'DragF'] = self.SC.get_DragF(DF.MET[i],\
                                                              DF.EarthFixedVEL[i],\
                                                              DF.atmosRho[i],\
                                                              DF.atmosM1[i])
                    # static
                    else:
                        DF.loc[i,'DragF'] = DF.dynPress[i] * DF.SC_static_AreaCd[i]


            # Orbital parameters Earth
            elif datatype == 'EarthOrbElements':
                if self.OE_TEME:
                    print('MRS:\t\tAdding Orbital Elements (Earth, true Equator/Equninox of date) to dataframe.')
                    statevecsTrueEq = self.transform_GCRF_TrueEq(JDs, statevecs)
                    EarthOrbElements = self.get_OrbElements(JDs, statevecsTrueEq, self.EARTH_GM)
                else:
                    print('MRS:\t\tAdding Orbital Elements (ICRF) to dataframe.')
                    EarthOrbElements = self.get_OrbElements(JDs, statevecs, self.EARTH_GM)
                DF['EarthOESMA'] = EarthOrbElements.semi_major_axis.m
                DF['EarthOEeccentricity'] = EarthOrbElements.eccentricity
                DF['EarthOERAAN'] = EarthOrbElements.longitude_of_ascending_node.degrees
                DF['EarthOEargPeriapsis'] = EarthOrbElements.argument_of_periapsis.degrees
                DF['EarthOEinclination'] = EarthOrbElements.inclination.degrees
                DF['EarthOEtrueAnomaly'] = EarthOrbElements.true_anomaly.degrees
                DF['EarthOEmeanAnomaly'] = EarthOrbElements.mean_anomaly.degrees

            elif datatype == 'EarthFPAHAvel':
                print('MRS:\t\tAdding FPA/HA/VEL (w.r.t. to inertial Earth) to dataframe.')
                FPA, HA = self.get_FPAHA(JDs, statevecs, frame='Earth')
                DF['EarthFPA'] = FPA * RAD2DEG
                DF['EarthHA'] = HA * RAD2DEG
                DF['EarthVEL'] = np.linalg.norm(statevecs[:,3:], axis=1)


            # SC position in Moon frame
            elif datatype == 'MoonFrame':
                print('MRS:\t\tAdding Moon positions to dataframe.')
                for i in range(lenDF):
                    SatMoonLLA, SatMoonXYZ, MoonPos, MoonVel, SatMoonXYZVel, R, SatMoonICRFPos, SatMoonICRFVel \
                        = self.transform_GCRF_MoonLLAplanetocentric(JDs[i], statevecs[i,:], 'ME')
                    DF.loc[i,['MoonLat', 'MoonLon', 'MoonAlt']] = SatMoonLLA
                    DF.loc[i,['MoonX', 'MoonY', 'MoonZ']] = SatMoonXYZ
                    DF.loc[i,['MoonVx', 'MoonVy', 'MoonVz']] = SatMoonXYZVel
                    DF.loc[i,['MoonICRFx', 'MoonICRFy', 'MoonICRFz']] = SatMoonICRFPos
                    DF.loc[i,['MoonICRFvx', 'MoonICRFvy', 'MoonICRFvz']] = SatMoonICRFVel
                    DF.loc[i,['MoonPosX', 'MoonPosY', 'MoonPosZ']] = MoonPos
                    DF.loc[i,['MoonVelX', 'MoonVelY', 'MoonVelZ']] = MoonVel

            # rotating pulsating system with Moon
            elif datatype == 'MoonRPS':
                print('MRS:\t\tAdding Rotating-Pulsating System positions to dataframe.')
                # check if Moon positions are already there, if not, complain
                if not 'MoonPosX' in DF.columns:
                    print('MRS:\t\tERROR: Moon Frame data needs to be added first; skipping MoonRPS.')
                    continue
                # RPS columns in DF
                RPScolumns = ['MoonRPSx', 'MoonRPSy', 'MoonRPSz']
                # get + write data
                DF[RPScolumns] = self.get_RPSpos(JDs, statevecs , updateMissionDF=1)

            # Orbital parameters around Moon
            elif datatype == 'MoonOrbElements':
                print('MRS:\t\tAdding Orbital Elements (Moon) to dataframe.')
                # check if Moon positions are already there, if not, complain
                if not 'MoonPosX' in DF.columns:
                    print('MRS:\t\tERROR: Moon Frame data needs to be added first; skipping MoonRPS.')
                    continue
                # Statevecs for Orb Elements in inertial frame
                #MoonStatevecs = statevecs - DF[['MoonPosX', 'MoonPosY', 'MoonPosZ','MoonVelX', 'MoonVelY', 'MoonVelZ' ]].to_numpy()
                # Statevecs for Orb Elements in Moon fixed frame
                MoonStatevecs = DF[['MoonX', 'MoonY', 'MoonZ','MoonVx', 'MoonVy', 'MoonVz' ]].to_numpy()
                MoonOrbElements = self.get_OrbElements(JDs, MoonStatevecs, self.MOON_GM)
                DF['MoonOESMA'] = MoonOrbElements.semi_major_axis.m
                DF['MoonOEeccentricity'] = MoonOrbElements.eccentricity
                DF['MoonOERAAN'] = MoonOrbElements.longitude_of_ascending_node.degrees
                DF['MoonOEargPeriapsis'] = MoonOrbElements.argument_of_periapsis.degrees
                DF['MoonOEinclination'] = MoonOrbElements.inclination.degrees
                DF['MoonOEtrueAnomaly'] = MoonOrbElements.true_anomaly.degrees


            elif datatype == 'MoonFPAHAvel':
                print('MRS:\t\tAdding FPA/HA/VEL (w.r.t. to Moon-frame) to dataframe.')
                Moonstatevecs = DF[['MoonX', 'MoonY', 'MoonZ','MoonVx', 'MoonVy', 'MoonVz']].to_numpy()
                FPA, HA = self.get_FPAHA(JDs, Moonstatevecs, frame='Moon')
                DF['MoonFPA'] = FPA * RAD2DEG
                DF['MoonHA'] = HA * RAD2DEG
                DF['MoonVEL'] = np.linalg.norm(DF[['MoonVelX', 'MoonVelY', 'MoonVelZ' ]].to_numpy()[:,3:], axis=1)


            elif datatype == 'UTC':
                print('MRS:\t\tAdding UTC time to dataframe.')
                DF['UTC'] = ts.tdb_jd(self.missionDF.JD_TBD.to_numpy()).utc_iso()

            elif datatype == 'Eclipse':
                print('MRS:\t\tAdding eclipse information to dataframe.')

                for i in range(lenDF):
                    self.planets = self.MD.forcesSettings.planets[DF.forcesID[i]]
                    DF.loc[i,['ShadowFunction']], DF.loc[i,['SC2SUN_x','SC2SUN_y','SC2SUN_z']], DF.loc[i,['OccultingBodies']] = self.get_sunVisibility(JDs[i], statevecs[i,:3])

            elif datatype == 'SRP':

                if not 'SC_static_Cr' in DF.columns:
                    DF = self.expand_DF(['SpacecraftStatic'], DF)

                if not 'ShadowFunction' in DF.columns:
                    DF = self.expand_DF(['Eclipse'], DF)

                print('MRS:\t\tAdding SRP force to dataframe.')

                for i in range(lenDF):
                    # add SRP only if required in forces settings
                    if self.MD.forcesSettings.SRP[DF.forcesID[i]]:
                        SRPforce = self.get_SRPForce(DF.ShadowFunction[i],\
                                                     DF.SC_static_Cr[i],\
                                                     DF.SC_static_SRParea[i],\
                                                     DF.loc[i,['SC2SUN_x','SC2SUN_y','SC2SUN_z']].to_numpy())
                        DF.loc[i,['SRP']] = np.linalg.norm(SRPforce)
                    else:
                        DF.loc[i, ['SRP']] = 0



            elif datatype == 'PlanetPositions':
                # Make list of Plaents used in Forcesettings, and making sure
                # that Earth is always first in list (because Earth Data used later on)
                PlanetsForceSettings = self.MD.forcesSettings.planets.to_numpy()
                PlanetList = ['Earth'] + [item for sublist in PlanetsForceSettings for item in sublist]
                PlanetList = list(dict.fromkeys(PlanetList))

                # loop through list
                for planet in PlanetList:
                    if planet=='Sun':
                        objectPos = self.sun.at(ts.tdb_jd(JDs)).position.m
                        DF[['Sun_ICRF_x','Sun_ICRF_y','Sun_ICRF_z']] = objectPos.T
                    elif planet=='Venus':
                        objectPos = self.venus.at(ts.tdb_jd(JDs)).position.m
                        DF[['Venus_ICRF_x','Venus_ICRF_y','Venus_ICRF_z']] = objectPos.T
                    elif planet=='Earth':
                        objectPos = self.earth.at(ts.tdb_jd(JDs)).position.m
                        DF[['Earth_ICRF_x','Earth_ICRF_y','Earth_ICRF_z']] = objectPos.T
                    elif planet=='Moon':
                        objectPos = self.moon.at(ts.tdb_jd(JDs)).position.m
                        DF[['Moon_ICRF_x','Moon_ICRF_y','Moon_ICRF_z']] = objectPos.T
                    elif planet=='Mars':
                        objectPos = self.mars.at(ts.tdb_jd(JDs)).position.m
                        DF[['Mars_ICRF_x','Mars_ICRF_y','Mars_ICRF_z']] = objectPos.T
                    elif planet=='Jupiter':
                        objectPos = self.jupiter.at(ts.tdb_jd(JDs)).position.m
                        DF[['Jupiter_ICRF_x','Jupiter_ICRF_y','Jupiter_ICRF_z']] = objectPos.T
                    else:
                        print('MRS:\t\tERROR: invalid planet name provided in forcesSettings.')
                        continue

            elif datatype == 'CompVNB_Earth':

                if not 'compX' in DF.columns:
                    print('MRS:\t\tERROR: no comparison trajectory loaded; skipping CompVNB_Earth.')
                    continue

                print('MRS:\t\tAdding VNB pos. rel. to comparison trajectory to DF (Earth orbit.')

                for i in range(lenDF):

                    # get VNB rotation matrix for given time/comparison state vector
                    VNBrot = self.get_VNBframe(JDs[i], \
                                               DF.loc[i,['compX', 'compY', 'compZ',\
                                                         'compVx', 'compVy', 'compVz']].astype('float').to_numpy() )
                    DF.loc[i,['EO_VNBx', 'EO_VNBy', 'EO_VNBz']] = \
                        VNBrot.T.dot(DF.loc[i,['compPosDiffx','compPosDiffy','compPosDiffz']].astype('float').to_numpy() )


            elif datatype == 'CompVNB_Moon':

                if not 'compX' in DF.columns:
                    print('MRS:\t\tERROR: no comparison trajectory loaded; skipping CompVNB_Moon.')
                    continue

                print('MRS:\t\tAdding VNB pos. rel. to comparison trajectory to DF (Moon orbit.')

                for i in range(lenDF):

                    # get VNB rotation matrix for given time/comparison state vector
                    VNBrot = self.get_VNBframe(JDs[i], \
                                               DF.loc[i,['compX', 'compY', 'compZ',\
                                                         'compVx', 'compVy', 'compVz']].astype('float').to_numpy(),\
                                               planet='Moon')
                    DF.loc[i,['MO_VNBx', 'MO_VNBy', 'MO_VNBz']] = \
                        VNBrot.T.dot(DF.loc[i,['compPosDiffx','compPosDiffy','compPosDiffz']].astype('float').to_numpy() )

            elif datatype == 'GuidanceVec':
                
                print('MRS:\t\tAdding guidance vector.')
                
                for i in range(lenDF):
                    DF.loc[i,['gVecX', 'gVecY', 'gVecZ']] = self.GO.get_guidance(JDs[i], statevecs[i], DF.MET[i], mode='TrueMET')
             
            elif datatype == 'RangeToLaunchsite':
                
                print('MRS:\t\tAdding range distance to launchsite.')
                
                for i in range(lenDF):
                    DF.loc[i,['range']] = self.get_EarthRangeToLaunchsite(JDs[i], statevecs[i])
                    
            
            elif datatype == 'GuidanceVecAngles':
                
                if self.GO.mode == 'static':
                    print('MRS:\t\tWARNING: no gVec angles available in static guidance; skipping.')
                    continue
       
                if not 'gVecX' in DF.columns:
                    DF = self.expand_DF(['GuidanceVec'], DF)
                    
                print('MRS:\t\tAdding guidance vector angle values.')
                    
                for i in range(lenDF):
                    gVecValues = self.GO.get_delta_gvec_to_vel(JDs[i], statevecs[i], \
                                    DF.loc[i,['gVecX', 'gVecY', 'gVecZ']].astype('float').to_numpy() )
                    
                    DF.loc[i,['gVec_Earth_ENU_abs_elev', 'gVec_Earth_ENU_abs_head']] = gVecValues[0,:2]
                    DF.loc[i,['gVec_EFvel_Earth_ENU_delta_elev', 'gVec_EFvel_Earth_ENU_delta_head']] = gVecValues[1,:2]
                    DF.loc[i,['gVec_SFvel_Earth_ENU_delta_elev', 'gVec_SFvel_Earth_ENU_delta_head']] = gVecValues[2,:2]
                    DF.loc[i,['gVec_Launch_ENU_abs_pitch', 'gVec_Launch_ENU_abs_head']] = gVecValues[3,:2]
                    DF.loc[i,['gVec_GCRF_abs_elev', 'gVec_GCRF_abs_head']] = gVecValues[4,:2]
                    DF.loc[i,['gVec_VUW_abs_elev', 'gVec_VUW_abs_head']] = gVecValues[5,:2]
                    DF.loc[i,['gVec_VNB_abs_elev', 'gVec_VNB_abs_head']] = gVecValues[6,:2]
                    
            
            # unknown kind of data requested
            else:
                print('MRS:\t\tERROR: unknown kind of data type requested: ', datatype)

        # clear temporary variables
        self.del_tempVars()

        # return modified dataframe
        return DF

    def load_externalMission(self, timevec, statevecs, name='extMission', atmosModel='', DEephemeris='DE421'):
        """
        This function generates a mission dataframe based on provided
        state vectors of an external mission. Further processing, e.g. adding
        atmospheric values is possible.

        Parameters
        ----------
        timevec : array of floats
            Julian Dates (TBD) of external mission.
        statevecs : array of flaots
            State vectors of external mission (in [m]/[m/s])
        name : string, optional
            Name of the external mission. The default is 'extMission'.
        atmosModel : string, optional
            Name of the atmospheric model to be used. The default is ''.

        Returns
        -------
        None.

        """


        # check sizes
        if timevec.shape[0] != statevecs.shape[0]:
            print('MRS:\t\tError loading comparison DF: times and state vecs do not match.')
            return 0
        if statevecs.shape[1] != 6:
            print('MRS:\t\tError loading comparison DF: state vecs incomplete.')
            return 0

        # check if timevec are UTC strings
        if type(timevec[0]) == str:
            # make vector of JD_TBD times
            JDvec = np.zeros(timevec.shape[0])
            # go through all rows
            for i in range(timevec.shape[0] ):
                UTCtime = datetime.datetime.fromisoformat(timevec[i])
                UTCtime = UTCtime.replace(tzinfo=utc)
                JDvec[i] = ts.from_datetime(UTCtime).tdb
        # if not, assume it's JD TBD time
        else:
            JDvec = timevec

        # set mission data (MD)
        self.MD.name = name
        print('MRS:\t\tImporting mission '+self.MD.name+'.')

        # make empty missionDF
        self.make_TempDF(int(0))
        self.missionDF = self.TempDF

        # fill data into tempDF
        self.missionDF['JD_TBD'] = JDvec
        self.missionDF[['x','y','z','vx','vy','vz']] = statevecs
        self.missionDF['atmosModel'] = atmosModel

        # set MET values
        self.missionDF['MET'] = (JDvec-JDvec[0]) * 24 * 3600

        # set not used values
        self.missionDF[['segmentID', 'segmentType', 'configID', 'propaMode','forcesID']] = 0

        # save kind of missionDF (1=MRS sim, 2=external data)
        self.missionDFtype = 2

        # load ephemeris
        self.DEephemeris = DEephemeris
        self.load_ephemeris()

        # set t0_JD
        self.MD.t0_JD = JDvec[0]

        # make a default spacecraft
        self.SC = MRSstaticSpacecraft()

        return None


    def load_comparisonStateVec(self, timevec, statevecs):
        """
        This function interpolates the state vectors from an external mission
        at the Julian Dates of the existing mission dataframe and and adds
        the interpolated state vectors to the mission dataframe.


        Parameters
        ----------
        timevec : array of floats
            Julian Dates (TBD) of external mission. If 0, timevec is assumed
            to be identical to mission data frame. Equal amount of state vecs
            is required.
        statevecs : array of flaots
            State vectors of external mission (in [m]/[m/s])

        Returns
        -------
        None.

        """

        # if timevec is provided, use the timestamps
        if np.sum(timevec) != 0:

            # check sizes
            if timevec.shape[0] != statevecs.shape[0]:
                print('MRS:\t\tError loading comparison DF: times and state vecs do not match.')
                return 0
            if statevecs.shape[1] != 6:
                print('MRS:\t\tError loading comparison DF: state vecs incomplete.')
                return 0

            # check if timevec are UTC strings
            if type(timevec[0]) == str:
                # make vector of JD_TBD times
                JDvec = np.zeros(timevec.shape[0])
                # go through all rows
                for i in range(timevec.shape[0] ):
                    UTCtime = datetime.datetime.fromisoformat(timevec[i])
                    UTCtime = UTCtime.replace(tzinfo=utc)
                    JDvec[i] = ts.from_datetime(UTCtime).tdb
            # if not, assume it's JD TBD time
            else:
                JDvec = timevec


            if JDvec[0]>self.missionDF.JD_TBD.iloc[0]:
                print('MRS:\t\tWARNING: Setting comparison start time to mission start time.')
                JDvec[0] = self.missionDF.JD_TBD.iloc[0]

            if JDvec[-1]<self.missionDF.JD_TBD.iloc[-1]:
                print('MRS:\t\tWARNING: Setting comparison end time to mission end time.')
                JDvec[-1] = self.missionDF.JD_TBD.iloc[-1]


            # set up bspline interpolation
            comparisonInterpolation = make_interp_spline(JDvec, statevecs, k=7)

            # get interpolated state vecs and put into dataframe
            comparisonStateVecs = comparisonInterpolation(self.missionDF['JD_TBD'].to_numpy())
            self.missionDF[['compX', 'compY', 'compZ', \
                            'compVx', 'compVy', 'compVz',\
                            ]] = comparisonStateVecs

        # no timevec provided; assumed to be identical to MRS data frame
        else:
            # check size
            if statevecs.shape[0] != len(self.missionDF):
                print('MRS:\t\tError loading comparison DF: state vecs no same length as DF.')
                return 0

            # copy statevecs
            self.missionDF[['compX', 'compY', 'compZ', \
                            'compVx', 'compVy', 'compVz',\
                            ]] = statevecs

        # calc error between state vectors (position)
        diffXYZ = self.missionDF[['x','y','z']].to_numpy() - self.missionDF[['compX','compY','compZ']].to_numpy()
        self.missionDF[['compPosDiffx','compPosDiffy','compPosDiffz']] = diffXYZ
        self.missionDF['compPosDiff'] = np.linalg.norm(diffXYZ, axis=1)

        # calc error between state vectors (velocity)
        diffVelXYZ = self.missionDF[['vx','vy','vz']].to_numpy() - self.missionDF[['compVx','compVy','compVz']].to_numpy()
        self.missionDF[['compVelDiffx','compVelDiffy','compVelDiffz']] = diffVelXYZ
        self.missionDF['compVelDiff'] = np.linalg.norm(diffVelXYZ, axis=1)


        print('MRS:\t\tSuccessfully imported comparison data.')

        return None

    def loadDataframes(self, missionDFfile='', eventsDFfile=''):
        """
        Loads external csv files to missionDF and eventsDF.

        Parameters
        ----------
        missionDFfile : string
            Relative path to missionDF csv-file.
        eventsDFfile : string
            Relative path to eventsDF csv-file.

        Returns
        -------
        None.

        """

        if missionDFfile:
            print('MRS:\t\tLoading missionDF from ', missionDFfile)
            self.missionDF = pd.read_csv(missionDFfile, sep=',')

        if eventsDFfile:
            print('MRS:\t\tLoading eventsDF from ', eventsDFfile)
            self.eventsDF = pd.read_csv(missionDFfile, sep=',')

        return None

    def exportDataframes(self, folder='./', usetimestamp=0, missionDFonly=0):
        """
        Saves the dataframes (mission & events) as csv files.

        Parameters
        ----------
        folder : string
            Relative path to folder where to save the files
        usetimestamp: int (0/1); default 0
            Wether to use a timestamp in the filename (1) or not (0)
        missionDFonly : int (0/1); default 0
            Whether to save only the missionDF (1) file or both (0)

        Returns
        -------
        None.

        """

        # get current date/time as string
        timestamp = datetime.datetime.now().strftime("%Y.%m.%d-%H.%M.%S")

        # prepare filename
        if usetimestamp:
            filename = folder + timestamp +' '+ self.MD.name
        else:
            filename = folder + self.MD.name

        if hasattr(self, 'missionDF'):
            filenameMission = filename +'_missionDF.csv'
            print('MRS:\t\tSaving missionDF to ', filenameMission)

            # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_csv.html
            self.missionDF.to_csv(filenameMission,
                                  index=False,
                                  sep=','
                                  )
        else:
            print('MRS:\t\tWARNING: no missionDF available for csv-export.')


        if hasattr(self, 'eventsDF') and not missionDFonly:
            filenameEvents = filename +'_eventsDF.csv'
            print('MRS:\t\tSaving eventsDF to ', filenameEvents)

            # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_csv.html
            self.eventsDF.to_csv(filenameEvents,
                                  index=False,
                                  sep=','
                                  )
        elif not hasattr(self, 'eventsDF') and not missionDFonly:
            print('MRS:\t\tWARNING: no eventsDF available for csv-export.')
        else:
            pass


        return None


    def return_StateVec(self):
        """
        Returns state vectors saved in mission dataframe

        Returns
        -------
        JDs : array of floats
            Julian Date (TBD)
        statevecs : n x 6 array of floats
            State vectors in GCRF.

        """

        # get timevec values
        JDs = self.missionDF.JD_TBD.to_numpy()

        # get statevec values
        statevecs = self.missionDF[['x', 'y', 'z', 'vx', 'vy', 'vz']].to_numpy()

        return JDs, statevecs

    def return_comparisonStateVec(self):
        """
        Returns the interpolated values of loaded comparison state vectors.

        Returns
        -------
        JDs : array of floats
            Julian Date (TBD)
        statevecs : n x 6 array of floats
            State vectors in GCRF.

        """
        #

        # get timevec values
        JDs = self.missionDF.JD_TBD.to_numpy()

        # get statevec values
        statevecs = self.missionDF[['compX', 'compY', 'compZ', 'compVx', 'compVy', 'compVz']].to_numpy()

        return JDs, statevecs

    def get_EventsList(self, eventNames=[]):
        """
        Checks mission dataframe agains provided list of events and, if found,
        interpoalates the event state vector and adds it to the event dataframe.

        Parameters
        ----------
        eventNames : list of strings
            List of event types to be searched for.

        Returns
        -------
        None.

        """

        # Potential events
        # - Mach1
        # - MaxQ
        # - closest point to Moon
        # - enters EI
        # - SC crashes on Earth
        # - Rocket Staging

        # loop through provided array of event names
        for eventType in eventNames:

            # skip if eventType already present in eventsDF
            if eventType in self.eventsDF.eventType.to_numpy():
                continue

            # get JD time for given eventType
            if eventType == 'Closest2Earth':

                if not 'EarthAlt' in self.missionDF.columns:
                    DF = self.expand_DFname(['EarthLLA'])

                JDevent = self.find_EventTime(0,'min','EarthAlt', eventType)

            elif eventtype == 'FartestFromEarth':

                if not 'EarthAlt' in self.missionDF.columns:
                    DF = self.expand_DFname(['EarthLLA'])

                JDevent = find_EventTime(0,'max','EarthAlt', eventType)

            else:
                print('MRS:\t\tERROR in get_EventList(): unknown kind of event type requested: ', eventtype)

            # in case the event was found, add it to eventsDF
            if JDevent:
                self.add_event(JDevent, np.zeros(6), eventType)

        # add MET values; round due to resolution of JD values
        self.eventsDF['MET'] = round((self.eventsDF['JD_TBD'] - self.MD.t0_JD) / SEC2DAYS, 4)

        # add event list from mission file
        if hasattr(self.MD, 'missionEvents'):
            # add JD values for given MET values
            self.MD.missionEvents['JD_TBD'] = self.MD.t0_JD + (self.MD.missionEvents.MET - self.MD.t0_MET) * SEC2DAYS
            # append to eventsDF
            self.eventsDF = self.eventsDF.append(self.MD.missionEvents[['MET','JD_TBD','eventType']], ignore_index=True)

        # get spacecraft events
        SCevents = self.SC.get_EventsList()
        # add JD to spacecraft events
        SCevents['JD_TBD'] = self.MD.t0_JD + (SCevents.MET - self.MD.t0_MET) * SEC2DAYS
        # add to eventsDF
        self.eventsDF = self.eventsDF.append(SCevents, ignore_index=True)

        # sort event list by JD
        self.eventsDF = self.eventsDF.sort_values(by=['JD_TBD']).reset_index(drop=True)

        # prepare interpolation of state vectors
        SVInterp = make_interp_spline(self.missionDF.JD_TBD.to_numpy(), \
                                            self.missionDF.loc[:,['x','y','z','vx','vy','vz']].to_numpy(), k=7)

        # names of colums to copy their value from missionDF
        colNames = ['segmentID', 'segmentType', 'configID', 'propaMode','forcesID', 'atmosModel']

        # add additional values to eventsDF
        for i in range(len(self.eventsDF)):

            # check if statevector is NOT provided (by adding x,y,z-values and comparing to 0)
            if not np.sum(self.eventsDF.loc[i,['x','y','z']]):
                # interpolate state vector values based on provided JD_TBD value
                self.eventsDF.loc[i,['x','y','z','vx','vy','vz']] = SVInterp(self.eventsDF.JD_TBD[i]).tolist()

            # get pointer to the closest missionDF entry ahead of event
            missionDFpointer = ((self.eventsDF.JD_TBD[i] - self.missionDF.JD_TBD.to_numpy())>=0).nonzero()[0][-1]

            # add values to eventsDF (needed to expand the DF later on)
            self.eventsDF.loc[i, colNames] = self.missionDF.loc[missionDFpointer, colNames]




    def find_EventTime(self, value, eventMode, DFcolumn, eventType):
        """
        Internal function.
        Looks for a specific position in a column mission dataframe and returns
        found time (as Julian Date).

        Parameters
        ----------
        value : float
            A reference value to be found in specified column (if needed)
        eventMode : string
            Kind of searched event mode. Possible are:
                max: looks for the max. value in column
                min: looks for the min. value in column
                above: finds the first occurence when column value is above
                    specified reference value
                below: finds the last occurence when column value is below
                    specified referencec value
        DFcolumn : string
            Name of column in mission dataframe (self.missionDF).
        eventType : string
            Event type (input not processed in current version)

        Returns
        -------
        JDevent : float
            Julian Date (TBD) of event occurence

        """
        # gets DF for events

        # prepare variables
        JDevent = 0

        if eventMode == 'max' or eventMode == 'min' :
            # find index of max value
            if eventMode == 'max':
                DFpointer = self.missionDF.loc[:,DFcolumn].idxmax()
            # find index of min value
            else:
                DFpointer = self.missionDF.loc[:,DFcolumn].idxmin()

            # check if pointer is first or last value of missionDF, return that JD value
            if DFpointer == 0 or DFpointer == len(self.missionDF)-1:
                # return valid JD value
                return self.missionDF.loc[DFpointer,'JD_TBD']

            # find offset of interpolated peak rel. to max value
            offset = self.get_QuadraticPeakInterpolation(self.missionDF.loc[DFpointer-1:DFpointer+1,'EarthAlt'].to_numpy())

            # get JD-differences to entry before (pre) and after (post)
            deltaJDpre = self.missionDF.loc[DFpointer,'JD_TBD'] - self.missionDF.loc[DFpointer-1,'JD_TBD']
            deltaJDpost= self.missionDF.loc[DFpointer+1,'JD_TBD'] - self.missionDF.loc[DFpointer,'JD_TBD']

            # if difference is bigger than 1%, two segments with different timestamps --> offset value not valid
            if np.abs(deltaJDpost-deltaJDpre)/deltaJDpre > 0.01:
                print('MRS:\t\tWARNING in find_EventLTime(): non-continous JD values for inpterpolation of: ', eventType)
                # return 0 (=event not found)
                return 0

            # calc new JD value for interpolated peak
            JDevent = self.missionDF.loc[DFpointer,'JD_TBD'] + offset * deltaJDpre


        elif eventMode == 'above':
            # finds first occurence where a DFcolumn is above given value

            # get array of index when required value was passed
            arrayValuesAbove = ((self.missionDF.loc[:,DFcolumn].to_numpy())-value>=0).nonzero()[0]

            # check that at least one index value is available; return 0 if not
            if not arrayValuesAbove.shape[0]:
                return 0

            # get first occurence
            missionDFpointer = arrayValuesAbove[0]

            # check if pointer is first value of missionDF;
            # if this is the case, return 0 (because no gradient can be calculated)
            if missionDFpointer == 0:
                return 0

            # calc gradient between DFpointer-value and the one before
            valueGradient = (self.missionDF.loc[missionDFpointer,  DFcolumn]-   \
                             self.missionDF.loc[missionDFpointer-1,DFcolumn]) / \
                            (self.missionDF.loc[missionDFpointer,  'JD_TBD']-   \
                             self.missionDF.loc[missionDFpointer-1,'JD_TBD'])

            # calc JD time delta rel. to value before passing the given value
            JDdelta = (value - self.missionDF.loc[missionDFpointer-1,DFcolumn])/valueGradient

            # linearly interpolated JD value for eent
            JDevent = self.missionDF.loc[missionDFpointer-1,'JD_TBD'] + JDdelta


        elif eventMode == 'below':
            # finds last occurence where a DFcolumn is right below given value

            # get array of index when required value is above threshold-values
            arrayValuesAbove = ((self.missionDF.loc[:,DFcolumn].to_numpy()-value)>=0).nonzero()[0]

            # check that at least one index value is available; return 0 if not
            if not arrayValuesAbove.shape[0]:
                return 0

            # get last occurence
            missionDFpointer = arrayValuesAbove[-1]

            # check if pointer is last value of missionDF;
            # if this is the case, return 0 (because no gradient can be calculated,
            # as there is no following value)
            if missionDFpointer == len(self.missionDF)-1 :
                return 0

            # calc gradient between DFpointer-value and the one before
            valueGradient = (self.missionDF.loc[missionDFpointer+1,  DFcolumn]-     \
                             self.missionDF.loc[missionDFpointer,DFcolumn]) / \
                            (self.missionDF.loc[missionDFpointer+1,  'JD_TBD']-     \
                             self.missionDF.loc[missionDFpointer,'JD_TBD'])

            # calc JD time delta rel. to value before passing the given value
            JDdelta = (value - self.missionDF.loc[missionDFpointer,DFcolumn])/valueGradient

            # linearly interpolated JD value for eent
            JDevent = self.missionDF.loc[missionDFpointer,'JD_TBD'] + JDdelta


        else:
            print('MRS:\t\tERROR in find_EventLTimeLinear(): unknown kind of event mode requested:', eventMode)


        return JDevent


    def get_QuadraticPeakInterpolation(self, magnitude):
        """
        Returns the position of an interpolated peak for a curve described by
        three points.
        Reference:
        https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html


        Parameters
        ----------
        magnitude : array of floats (3x1)
            Amplitudes of three neighboring values where the peak is estimated.

        Returns
        -------
        peakPosition : float
            Peak-position as offset from central x-value (factor 0-0.5 w.r.t to
                the distance to x-1 and x+1

        """

        peakPosition = 1/2 * (magnitude[0] - magnitude[2]) / \
                       (magnitude[0] - 2*magnitude[1] + magnitude[2])

        return peakPosition




