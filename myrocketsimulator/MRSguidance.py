#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 19:58:52 2023

@author: thibault
"""




import numpy as np
from skyfield.api import load, Distance, wgs84
from skyfield.positionlib import Geocentric
from skyfield.framelib import itrs, true_equator_and_equinox_of_date
from skyfield.toposlib import ITRSPosition

# load Skyfield timescale
ts = load.timescale()


# constants
EARTH_ROT_SPEED = 7292115.1467e-11 # [rad/s]
DEG2RAD = np.pi/180.
RAD2DEG = 180./np.pi

# frame IDs
F_Earth_ENU_abs = 1
F_EFvel_Earth_ENU_delta = 2
F_SFvel_Earth_ENU_delta = 3
F_Launch_ENU_abs = 4
F_GCRF_abs = 5
F_VUW_abs = 6
F_VNB_abs = 7
F_none = 10

F_SM0_abs = 100
F_SM1_abs = 101
F_SM2_abs = 102
F_SM3_abs = 103
F_SM4_abs = 104
F_SM5_abs = 105
F_SM6_abs = 106
F_SM7_abs = 107
F_SM8_abs = 108
F_SM9_abs = 109
   
# colum IDs
I_MET = 0
I_frame = 1
I_end_value = 2
I_start_value = 3
I_gradient = 4   

class Guidance():
    """
    The Guidance() class has one primary goal: provide a direction for the thrust
    vector of the spacecraft. The direction is described by two tables (heading 
    and elevation) in different reference frames. The output frame is always
    GCRF. 
    
    The Guidance() class can also be used outside of MyRocketSimulator. This allows
    other simulation tools to take advantage of its features. 
    
    Contrary to MRS, the Guidance() class completely relies on numpy for arrays
    and tables in order to increase its speed. 
    
    """
    
    def __init__(self, gd):
        """
        Initializes the MRS guidance object.

        Parameters
        ----------
        gd : class
            Guidance data stored in form of a class. 
    

        Returns
        -------
        None.

        """
        
        self.guidancename = gd.name # save guidance name
        self.gd = gd # save spacraft data
        self.mode = 'active'
        self.gElevPointer = -1 # pointer to elevation guidance table
        self.gHeadPointer = -1 # pointer to heading guidance table
        self.ENU_liftoff = np.eye(3) # ENU at lift-off
        
        # empty stable matrix frames
        self.SMframes = np.zeros((10,3,3))
        self.SMframes[:, None, None] = np.eye(3)
        
        # initialize spacecraft elements
        self.init_guidance()
        
        return None
        
    def init_guidance(self):
        """
        Initializes the guidance tables for elevation and heading.

        Returns
        -------
        None.

        """
        
        # Update guidance tables
        self.gd.gElevTab = self.init_gTables(self.gd.gElevTab)
        self.gd.gHeadTab = self.init_gTables(self.gd.gHeadTab)
        
        # empty rotation matrix from GCRF to ITRF
        self.GCRF2ITRF = np.eye(3)
        
        return None
        
    def init_ENU_liftoff(self, JDnow, lla):
        """
        Sets up the ENU frame for the launchsite. Not using actual gravity vector.

        Parameters
        ----------
        JDnow : float
            Julian data (TDB) at launch of spacecraft.
        lla : array of floats
            Latitude [°], longitude [°], geocentric altitude [m].

        Returns
        -------
        None.

        """
     
        t = ts.tdb_jd(JDnow)
        itrsXYZ = wgs84.latlon(lla[0], lla[1], elevation_m=lla[2]).itrs_xyz
        p = ITRSPosition(itrsXYZ) # ITRS position
        self.ENU_liftoff = self.get_ENUvec_Earth(JDnow, p.at(t).position.m, frame='WGS84')
        
        return None  
            
    def init_gTables(self, gTab):
        """
        Initializes distinct guidance tables (either elevation or heading).
        Primary tasks are to store start values and the gradient per guidance
        segment.

        Parameters
        ----------
        gTab : 2D array of floats
            Guidance table (either elevation or heading).

        Returns
        -------
        gTab : 2D array of floats
            Guidance table (either elevation or heading) with set up start values
            and gradients.

        """
        
        # add two columns for start + gradient
        gTab = np.append(gTab,np.empty((len(gTab),2))*float('nan'), axis=1)
        
        # for first guidance segment, start value = end value
        gTab[0, I_start_value] = gTab[0, I_end_value] 
        gTab[0, I_gradient] = 0.
       
        # insert start values + gradient when possible
        for i in range(1,len(gTab)):
            if gTab[i, I_frame] == gTab[i-1, I_frame]:
                # set start value equal to end value of previous segment
                gTab[i,I_start_value] = gTab[i-1, I_end_value] 
                    
                # set gradient of value (but not in last, because its fixed at start value)
                if i != len(gTab)-1:
                    gTab[i,I_gradient] = \
                        (gTab[i,I_end_value]-gTab[i, I_start_value])/ \
                        (gTab[i+1, I_MET]-gTab[i, I_MET] )
                # only for last segment 
                else:
                    gTab[i, I_gradient] = 0
            
            # in case the previous guidance segment was not defined, use
            # end value as start value, gradient = 0
            elif gTab[i-1, I_frame] == F_none:
                gTab[i,I_start_value] = gTab[i, I_end_value] 
                gTab[i, I_gradient] = 0.
                
            # in case of stable matrix frames (not following F_none or the identical 
            # type of SMF), set start=end and gradient = 0.
            elif gTab[i, I_frame] >= 100 and gTab[i, I_frame] < 110:
                gTab[i,I_start_value] = gTab[i, I_end_value] 
                gTab[i, I_gradient] = 0.
            
        
        return gTab
    
    def get_guidance(self, JDnow, y, MET, mode='TrueMET'):
        """
        This method returns the guidance vector in GCRF. It can be called 
        in mode 'TrueMET' or in mode 'intermediateMET'. Ahead of calling the 
        method get_gVec(), it may try to calculate missing start values for
        the current guidance segments.
        
        Mode 'TrueMET' is called to identify the current guidance segments (for 
        both elevation and heading) and calculate start values of segments if 
        it couldn't be pre-calculated during the initialization. This mode can
        be used at all times. 
        
        The mode 'intermediateMET' relies on the stored segment pointers for
        the elevation and heading. This mode is only recommended to be called
        during integration steps. These steps should not be longer than the smallest
        time resulotion of the gudiance data, e.g., 100 ms in case guidance 
        segments are defined at a resolution of 1/10 seconds.

        Parameters
        ----------
        JDnow : float
           Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF.
        MET : float
            Mission Elapsed Time.
        mode : string, optional
            Mode of guidance determination. The default is 'TrueMET'.

        Returns
        -------
        gVec : array of floats
            Guidance vector in GCRF.

        """
        
        # clean up MET value (10 micro seconds resolution)
        MET = round(MET,5)
        
        # check valid MET value is provided
        if MET < self.gd.gElevTab[0, I_MET] or MET < self.gd.gHeadTab[0, I_MET]:
            print('MRS:\t\tWARNING get_guidance(): MET provided before guidance starts.')
            # set MET to latest guidance start
            MET = max(self.gd.gElevTab[0, I_MET], self.gd.gHeadTab[0, I_MET])
        
        elif MET >= self.gd.gElevTab[-1, I_MET] or MET >= self.gd.gHeadTab[-1, I_MET]:
            gVec = y[3:]/np.linalg.norm(y[3:])
            return gVec
        
        
        # if not using guidance for intermediate interpolator steps but 
        # for true MET values, additional tasks need to be done.
        if mode=='TrueMET':
            self.gElevPointer = ((MET-self.gd.gElevTab[:,I_MET])>=0).nonzero()[0][-1]
            self.gHeadPointer = ((MET-self.gd.gHeadTab[:,I_MET])>=0).nonzero()[0][-1]
        
            # set rotation from ITRF to GCRF 
            self.GCRF2ITRF = itrs.rotation_at(ts.tdb_jd(JDnow))
        
            # check if start value and/or gradient need to be calculated (elevation)
            if np.isnan(self.gd.gElevTab[self.gElevPointer,[I_start_value, I_gradient]]).any()\
                and self.gd.gElevTab[self.gElevPointer,I_frame] != F_none:
                
                # first call of segment should be at MET of segment start
                if MET == self.gd.gElevTab[self.gElevPointer, I_MET]:
                    
                    # go to previous guidance segment 
                    self.gElevPointer -= 1
                    # only use previous Heading segment if the current hasn't a start value either
                    if np.isnan(self.gd.gHeadTab[self.gHeadPointer, [I_start_value,I_gradient]]).any():
                        self.gHeadPointer -= 1
                        use_previous_HeadPointer = 1
                    else:
                        use_previous_HeadPointer = 0
                    
                    # get gVec using the previous segment
                    gVec_previousSegment = self.get_gVec(JDnow, y, MET)
                    
                    # get gVec elev/head in different frames
                    gVec_values = self.get_delta_gvec_to_vel(JDnow, y, gVec_previousSegment)
                    
                    # update back to current guidance segment
                    self.gElevPointer += 1
                    if use_previous_HeadPointer == 1:
                        self.gHeadPointer += 1
                        use_previous_HeadPointer = 0
                    
                    # store relevant value as start_value
                    if self.gd.gElevTab[self.gElevPointer, I_frame] == F_Earth_ENU_abs:
                        self.gd.gElevTab[self.gElevPointer, I_start_value] = gVec_values[0,0]
                    elif self.gd.gElevTab[self.gElevPointer, I_frame] == F_EFvel_Earth_ENU_delta:
                        self.gd.gElevTab[self.gElevPointer, I_start_value] = gVec_values[1,0]
                    elif self.gd.gElevTab[self.gElevPointer, I_frame] == F_SFvel_Earth_ENU_delta:
                        self.gd.gElevTab[self.gElevPointer, I_start_value] = gVec_values[2,0]
                    elif self.gd.gElevTab[self.gElevPointer, I_frame] == F_Launch_ENU_abs:
                        self.gd.gElevTab[self.gElevPointer, I_start_value] = gVec_values[3,0]
                    elif self.gd.gElevTab[self.gElevPointer, I_frame] == F_GCRF_abs:
                        self.gd.gElevTab[self.gElevPointer, I_start_value] = gVec_values[4,0]    
                    elif self.gd.gElevTab[self.gElevPointer, I_frame] == F_VUW_abs:
                        self.gd.gElevTab[self.gElevPointer, I_start_value] = gVec_values[5,0]    
                    elif self.gd.gElevTab[self.gElevPointer, I_frame] == F_VNB_abs:
                        self.gd.gElevTab[self.gElevPointer, I_start_value] = gVec_values[6,0]  
                    else:
                        self.gd.gElevTab[self.gElevPointer, I_start_value] = \
                            self.gd.gElevTab[self.gElevPointer, I_end_value] 
                        
                    # calc gradient
                    self.gd.gElevTab[self.gElevPointer, I_gradient] = \
                        (self.gd.gElevTab[self.gElevPointer, I_end_value] - self.gd.gElevTab[self.gElevPointer, I_start_value]) / \
                        (self.gd.gElevTab[self.gElevPointer+1, I_MET] - self.gd.gElevTab[self.gElevPointer, I_MET])
                
                    
                    
                else:
                    print('MRS:\t\tERROR get_guidance(): MET after segment start in mode TrueMET.')
                    print('MRS:\t\tWARNING get_guidance(): Using end_value as start_value, gradient=0.')
                    self.gd.gElevTab[self.gElevPointer, I_start_value] = \
                        self.gd.gElevTab[self.gElevPointer, I_end_value]
                    self.gd.gElevTab[self.gElevPointer, I_gradient] = 0.
    
        
            # check if start value and/or gradient need to be calculated (heading)
            if np.isnan(self.gd.gHeadTab[self.gHeadPointer,[I_start_value, I_gradient]]).any() \
                and self.gd.gHeadTab[self.gHeadPointer,I_frame] != F_none:
                # first call of segment should be at MET of segment start
                if MET == self.gd.gHeadTab[self.gHeadPointer, I_MET]:
                    
                    # go to previous guidance segment 
                    self.gHeadPointer -= 1
                    # only use previous Heading segment if the current hasn't a start value either
                    if np.isnan(self.gd.gElevTab[self.gHeadPointer, [I_start_value,I_gradient]]).any():
                        self.gHeadPointer -= 1
                        use_previous_ElevPointer = 1
                    else:
                        use_previous_ElevPointer = 0
                    
                    # get gVec using the previous segment
                    gVec_previousSegment = self.get_gVec(JDnow, y, MET)
                    
                    # get gVec elev/head in different frames
                    gVec_values = self.get_delta_gvec_to_vel(JDnow, y, gVec_previousSegment)
                    
                    # update back to current guidance segment
                    self.gHeadPointer += 1
                    if use_previous_ElevPointer == 1:
                        self.gElevPointer += 1
                        use_previous_ElevPointer = 0
                    
                    # store relevant value as start_value
                    if self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_Earth_ENU_abs:
                        self.gd.gHeadTab[self.gHeadPointer, I_start_value] = gVec_values[0,1]
                    elif self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_EFvel_Earth_ENU_delta:
                        self.gd.gHeadTab[self.gHeadPointer, I_start_value] = gVec_values[1,1]
                    elif self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_SFvel_Earth_ENU_delta:
                        self.gd.gHeadTab[self.gHeadPointer, I_start_value] = gVec_values[2,1]
                    elif self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_Launch_ENU_abs:
                        self.gd.gHeadTab[self.gHeadPointer, I_start_value] = gVec_values[3,1]
                    elif self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_GCRF_abs:
                        self.gd.gHeadTab[self.gHeadPointer, I_start_value] = gVec_values[4,1]    
                    elif self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_VUW_abs:
                        self.gd.gHeadTab[self.gHeadPointer, I_start_value] = gVec_values[5,1]  
                    elif self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_VNB_abs:
                         self.gd.gHeadTab[self.gHeadPointer, I_start_value] = gVec_values[6,1]  
                    else:
                        self.gd.gHeadTab[self.gHeadPointer, I_start_value] = \
                            self.gd.gHeadTab[self.gHeadPointer, I_end_value] 
                        
                    # calc gradient
                    self.gd.gHeadTab[self.gHeadPointer, I_gradient] = \
                        (self.gd.gHeadTab[self.gHeadPointer, I_end_value] - self.gd.gHeadTab[self.gHeadPointer, I_start_value]) / \
                        (self.gd.gHeadTab[self.gHeadPointer+1, I_MET] - self.gd.gHeadTab[self.gHeadPointer, I_MET])
                else:
                    print('MRS:\t\tERROR get_guidance(): MET after segment start in mode TrueMET.')
                    print('MRS:\t\tWARNING get_guidance(): Using end_value as start_value, gradient=0.')
                    self.gd.gHeadTab[self.gHeadPointer, I_start_value] = \
                        self.gd.gHeadTab[self.gHeadPointer, I_end_value]
                    self.gd.gHeadTab[self.gHeadPointer, I_gradient] = 0.
        
        
        
        # if called by interpolator
        elif mode=='intermediateMET':
            # check pointers are available; if not, guidance = velocity 
            if not self.gElevPointer>-1 or not self.gHeadPointer>-1:
                print('MRS:\t\tERROR get_guidance(): no guidance table pointers in mode intermediateMET.')
                # calc guidance vector 
                gVec = y[3:]/np.linalg.norm(y[3:])
                return gVec
                
        
        # calculate guidance vector for given table pointers and MET
        gVec = self.get_gVec(JDnow, y, MET)
        
        return gVec
        
        
    def get_gVec(self, JDnow, y, MET):
        """
        Calculates the get_gVec guidance vector using the internally stored
        guidance pointers. At first, the elevation and heading angles are 
        determined. Afterwards, and in function of the used frames, the actual
        guidance vector in GCRF is determined.

        Parameters
        ----------
        JDnow : float
           Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF.
        MET : float
            Mission Elapsed Time.

        Returns
        -------
        gVec : array of floats
            Guidance vector in GCRF.

        """
        
        # calculate only if guidance is in used (i.e. frame is not F_none)
        
        if not self.gd.gElevTab[self.gElevPointer, I_frame] == F_none \
            and not self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_none:
        
            # get FPA/HA in Earth-frame and (inertial) space-frame
            SF_FPA, SF_HA, ENU =  self.get_FPAHA(JDnow, y, frame='Earth')
            EF_FPA, EF_HA, _, EF_VEL  = self.get_FPAHAVEL_EF(JDnow, y)    
            
            #
            # elevation calculations
            #
            
            # calc current elevation angle value in current frame/segment
            gElevNow = self.gd.gElevTab[self.gElevPointer, I_start_value] + \
                       self.gd.gElevTab[self.gElevPointer, I_gradient] * \
                       (MET-self.gd.gElevTab[self.gElevPointer, I_MET])
            gElevNow *= DEG2RAD
    
            # adjust elevation value to used reference value
            if self.gd.gElevTab[self.gElevPointer, I_frame] == F_EFvel_Earth_ENU_delta:
                gElevFrame = gElevNow + EF_FPA 
            elif self.gd.gElevTab[self.gElevPointer, I_frame] == F_SFvel_Earth_ENU_delta:
                gElevFrame = gElevNow + SF_FPA 
            else: 
                gElevFrame = gElevNow
          
            #
            # heading calculations
            #
            
            # calc current elevation angle value 
            gHeadNow = self.gd.gHeadTab[self.gHeadPointer, I_start_value] + \
                       self.gd.gHeadTab[self.gHeadPointer, I_gradient] * \
                       (MET-self.gd.gHeadTab[self.gHeadPointer, I_MET])
            gHeadNow *= DEG2RAD
        
            # adjust heading value to used reference value
            if self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_EFvel_Earth_ENU_delta:
                gHeadFrame = gHeadNow + EF_HA 
            elif self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_SFvel_Earth_ENU_delta:
                gHeadFrame = gHeadNow + SF_HA 
            else: 
                gHeadFrame = gHeadNow
                
        # no guidance
        else:
            # gVec = velocity vector direction in GCRF
            gVec = y[3:]/np.linalg.norm(y[3:])
            return gVec
            
      
        #
        # get gVec in function of used reference frames for Elevation and Heading
        #
        
        # Elevation:
        # - F_Earth_ENU_abs 
        # - F_EFvel_Earth_ENU_delta 
        # - F_SFvel_Earth_ENU_delta 
        # Heading:
        # - F_Earth_ENU_abs 
        # - F_EFvel_Earth_ENU_delta 
        # - F_SFvel_Earth_ENU_delta 
        if self.gd.gElevTab[self.gElevPointer, I_frame] in \
            (F_Earth_ENU_abs, F_EFvel_Earth_ENU_delta, F_SFvel_Earth_ENU_delta)\
            and self.gd.gHeadTab[self.gHeadPointer, I_frame] in \
            (F_Earth_ENU_abs, F_EFvel_Earth_ENU_delta, F_SFvel_Earth_ENU_delta):
               
            gVec = self.get_gVec_ENU_abs(JDnow, y, gElevFrame, gHeadFrame, ENU)
        
        # Elevation:
        # - F_Launch_ENU_abs 
        # Heading:
        # - F_Earth_ENU_abs 
        # - F_EFvel_Earth_ENU_delta 
        # - F_SFvel_Earth_ENU_delta 
        
        elif self.gd.gElevTab[self.gElevPointer, I_frame] == F_Launch_ENU_abs \
            and self.gd.gHeadTab[self.gHeadPointer, I_frame] in \
            (F_Earth_ENU_abs, F_EFvel_Earth_ENU_delta, F_SFvel_Earth_ENU_delta):
                
            gVec = self.get_gVec_Earth_ENU_abs_LaunchPitch(JDnow, y, gElevFrame, gHeadFrame, ENU)
        
        
        # Elevation:
        # - F_Launch_ENU_abs 
        # Heading:
        # - F_Launch_ENU_abs 
        
        elif self.gd.gElevTab[self.gElevPointer, I_frame] == F_Launch_ENU_abs \
            and self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_Launch_ENU_abs:
                
            gVec = self.get_gVec_ENU_abs(JDnow, y, (90*DEG2RAD-gElevFrame), gHeadFrame, self.ENU_liftoff)
        
        
        # Elevation:
        # - F_GCRF_abs 
        # Heading:
        # - F_GCRF_abs  
        
        elif self.gd.gElevTab[self.gElevPointer, I_frame] == F_GCRF_abs \
            and self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_GCRF_abs:
                
            # GCRF ENU is given by np.eye(3)
            # Heading value is transformed from mathematical notation (Right ascension) to compass
            gVec = self.get_gVec_ENU_abs(JDnow, y, gElevFrame, (90*DEG2RAD-gHeadFrame), np.eye(3))
            
            
        # Elevation:
        # - F_VUW_abs 
        # Heading:
        # - F_VUW_abs  
        
        elif self.gd.gElevTab[self.gElevPointer, I_frame] == F_VUW_abs \
            and self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_VUW_abs:
                
            # Heading value is transformed from mathematical notation (Right ascension) to compass
            gVec = self.get_gVec_ENU_abs(JDnow, y, gElevFrame, (90*DEG2RAD-gHeadFrame), \
                                         self.get_VUW(JDnow, y))
                
        
        # Elevation:
        # - F_VNB_abs 
        # Heading:
        # - F_VNB_abs  
        
        elif self.gd.gElevTab[self.gElevPointer, I_frame] == F_VNB_abs \
            and self.gd.gHeadTab[self.gHeadPointer, I_frame] == F_VNB_abs:
                
            # Heading value is transformed from mathematical notation (Right ascension) to compass
            gVec = self.get_gVec_ENU_abs(JDnow, y, gElevFrame, (90*DEG2RAD-gHeadFrame), \
                                         self.get_VNB(JDnow, y))
            
    
        # Elevation:
        # - F_SMx_abs
        # Heading:
        # - F_SMx_abs
        
        elif self.gd.gElevTab[self.gElevPointer, I_frame] == self.gd.gHeadTab[self.gElevPointer, I_frame] \
            and self.gd.gElevTab[self.gElevPointer, I_frame] >= F_SM0_abs \
            and self.gd.gElevTab[self.gElevPointer, I_frame] <= F_SM9_abs:
            
            # get ENU frame from stable matrix frames array for given index
            # Heading value is transformed from mathematical notation (used in SMF; right ascension) to compass (used in ENU)
            gVec = self.get_gVec_ENU_abs(JDnow, y, gElevFrame, (90*DEG2RAD-gHeadFrame), \
                                         self.SMframes[int(self.gd.gElevTab[self.gElevPointer, I_frame]-100)])
            
        
        
        # no valid combination of frames for elevation and heading 
        else:
            print('MRS:\t\tERROR get_gVec(): no valid combination of heading/elevation frames.')

            # calc guidance vector 
            gVec = y[:3]/np.linalg.norm(y[:3])
            
        
        
        return gVec



    def get_gVec_ENU_abs(self, JDnow, y, gElevFrame, gHeadFrame, ENU=0):
        """
        Returns the guidance vector gVec defined by its elevation and heading
        angles for given frame.        

        Parameters
        ----------
        JDnow : float
           Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF.
        gElevFrame : float
            Elevation angle in provided ENU frame [rad].
        gHeadFrame : float
            Heading angle in provided ENU frame [rad].
        ENU : 3x3 array of floats, optional
            A rotation matrix describing a frame in GCRF. The default is 0.
            ENU stands for East/North/Up, but may be any kind of inertial or
            non-inertial frames defined by three axes.
            In case no frame is provided, the local Earth-centered ENU frame 
            is assumed (not recommended).

        Returns
        -------
        gVec : array of floats
            Guidance vector in GCRF.

        """
        
        # get ENU if not provided
        if not ENU.any():
            ENU = self.get_ENUvec_Earth(JDnow, y)
            
        # calc v-vector projection on north/east-plane based on azimuth
        gVec_projected = ENU[:,0] * np.sin(gHeadFrame) + ENU[:,1] * np.cos(gHeadFrame)

        # calc v vector using v_projected and elevation     
        gVec = gVec_projected * np.cos(gElevFrame) + ENU[:,2] * np.sin(gElevFrame)
        
        # make sure length is 1 
        gVec /= np.linalg.norm(gVec)
       
        return gVec
        
    
    def get_gVec_Earth_ENU_abs_LaunchPitch(self, JDnow, y, gElevFrame, gHeadFrame, ENU=0):
        """
        Returns the guidance vector gVec. The heading angle is defined in the 
        provided frame (ENU), whereas the angle is used as pitch angle measured
        from the Z-axis of the pre-determiend ENU frame at liftoff (ENU_liftoff).

        Parameters
        ----------
        JDnow : float
           Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF.
        gElevFrame : float
            Pitch angle measured from UP of ENU_liftoff [rad].
        gHeadFrame : float
            Heading angle in provided ENU frame [rad].
        ENU : 3x3 array of floats, optional
            A rotation matrix describing a frame in GCRF. The default is 0.
            ENU stands for East/North/Up, but may be any kind of inertial or
            non-inertial frames defined by three axes.
            In case no frame is provided, the local Earth-centered ENU frame 
            is assumed (not recommended).

        Returns
        -------
        gVec : array of floats
            Guidance vector in GCRF.

        """
        
        # get ENU if not provided
        if not ENU.any():
            ENU = self.get_ENUvec_Earth(JDnow, y)
        
        # calc v-vector projection on north/east-plane based on azimuth
        gVec_projected = ENU[:,0] * np.sin(gHeadFrame) + ENU[:,1] * np.cos(gHeadFrame)
       
        # calculate the required elevation angle in local Earth ENU in order
        # to comply to the pitch angle relative to the launch ENU 
        
        # intermediate calculatios
        a = gVec_projected.dot(self.ENU_liftoff[:,2])
        b = ENU[:,2].dot(self.ENU_liftoff[:,2])
        c = np.cos(gElevFrame)
        
        # round values; makes sure to maintain a**2+b**2 >= c**2
        a = round(a,10)
        b = round(b,10)
        
        # calculate elevation in local ENU frame
        ENU_elev = 2*np.arctan((b-np.sqrt(a**2+b**2-c**2))/ (a+c))
        
        # calc v vector using v_projected and elevation     
        gVec = gVec_projected * np.cos(ENU_elev) + ENU[:,2] * np.sin(ENU_elev)
        
        # make sure length is 1 
        gVec /= np.linalg.norm(gVec)
        
        return gVec
        

    def get_delta_gvec_to_vel(self, JDnow, y, gVec):
        """
        This method is called to determine the heading and elevation values
        for a given guidance vector in different frames. Generally, this 
        method is only called to determine the start value of a guidance segment 
        if a different reference frame was used in the previous segment. 
        

        Parameters
        ----------
        JDnow : float
           Current Julian Date (TBD).
        y : array of floats
            State vector in GCRF.
        gVec : array of floats
            Guidance vector in GCRF.

        Returns
        -------
        gVec_values : array of floats
            Numerous angles (in degrees) describing the elevation and heading
            values of the provided guidance vector in different reference frames. 

        """
        
        # get angle values for velocity vec and gVec in EF and SF
        
        # make a state vector with gVec
        y_gVec = np.append(y[:3], gVec)
           
        # values in Earth-fixed frame
        EF_FPA_vel, EF_HA_vel, _, EF_EFVEL = self.get_FPAHAVEL_EF(JDnow, y)
        EF_FPA_gVec, EF_HA_gVec, _ = self.get_FPAHA(JDnow, y_gVec)
            
        # values in space-fixed frame
        SF_FPA_vel, SF_HA_vel, _ = self.get_FPAHA(JDnow, y, frame='Earth')
        SF_FPA_gVec, SF_HA_gVec, _ = self.get_FPAHA(JDnow, y_gVec, frame='Earth')
        
        # values in launch ENU
        LaunchENU_FPA_gVec, LaunchENU_HA_gVec, _ = self.get_FPAHA(JDnow, y_gVec, frame='launchsite')
        
        # values in GCRF
        GCRF_FPA_gVec, GCRF_HA_gVec, _ = self.get_FPAHA(JDnow, y_gVec, frame='gcrf')
        
        # values in VUW
        VUW_FPA_gVec, VUW_HA_gVec, _ = self.get_FPAHA(JDnow, y_gVec, frame='vuw')
        
        # values in VNB
        VNB_FPA_gVec, VNB_HA_gVec, _ = self.get_FPAHA(JDnow, y_gVec, frame='vnb')
        
        
        # calc angles for different frames
        
        # F_Earth_ENU_abs (1)
        F_Earth_ENU_abs_elev = SF_FPA_gVec
        F_Earth_ENU_abs_head = SF_HA_gVec
        
        # F_EFvel_Earth_ENU_delta (2)
        F_EFvel_Earth_ENU_delta_elev = EF_FPA_gVec - EF_FPA_vel
        F_EFvel_Earth_ENU_delta_head = EF_HA_gVec - EF_HA_vel
        
        # F_SFvel_Earth_ENU_delta (3)
        F_SFvel_Earth_ENU_delta_elev = SF_FPA_gVec - SF_FPA_vel
        F_SFvel_Earth_ENU_delta_head = SF_HA_gVec - SF_HA_vel
        
        # F_Launch_ENU_abs (4)
        F_Launch_ENU_abs_pitch = (np.pi/2 -LaunchENU_FPA_gVec) # angle from launch-UP
        F_Launch_ENU_abs_head = LaunchENU_HA_gVec
        
        # F_GCRF_abs (5)
        F_GCRF_abs_elev = GCRF_FPA_gVec
        F_GCRF_abs_head = (np.pi/2 - GCRF_HA_gVec) # transform to right ascension
        
        # F_VUW_abs (6)
        F_VUW_abs_elev = VUW_FPA_gVec
        F_VUW_abs_head = (np.pi/2 - VUW_HA_gVec) # transform to right ascension
        
        # F_VNB_abs (7)
        F_VNB_abs_elev = VNB_FPA_gVec
        F_VNB_abs_head = (np.pi/2 - VNB_HA_gVec) # transform to right ascension
        
        # put all values in some array and transform to degress
        gVec_values = np.array([
            [F_Earth_ENU_abs_elev, F_Earth_ENU_abs_head],
            [F_EFvel_Earth_ENU_delta_elev, F_EFvel_Earth_ENU_delta_head],
            [F_SFvel_Earth_ENU_delta_elev, F_SFvel_Earth_ENU_delta_head],
            [F_Launch_ENU_abs_pitch, F_Launch_ENU_abs_head],
            [F_GCRF_abs_elev, F_GCRF_abs_head],
            [F_VUW_abs_elev,F_VUW_abs_head],
            [F_VNB_abs_elev,F_VNB_abs_head]
            ]) * RAD2DEG
        
        return gVec_values





    def get_ENUvec(self, JDnow, y):
        """
        Internal function. Not vectorized in MRSguidance for performance.
        Calculates local ENU vectors for provided statevec's positions in an
        object's frame.
        Following use of these ENU vectors require data processing in object's
        frame (e.g. for HA/FPA)
        Assumes round object.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD).
        y : array of floats
            State vector in the reference frame of planet.


        Returns
        -------
        ENU : 3x3 array of floats
            East, North, Up vectors / rotation matrices in local ref. system

        """

        # get UP vectors from position vetor
        UP = y[:3] * 1.0
        UP /= np.linalg.norm(UP)

        # get EAST vectors
        EAST = np.array([-y[1], y[0], 0])
        EAST /= np.linalg.norm(EAST)

        # get NORTH vectors
        NORTH = np.cross(UP, EAST)
        NORTH /= np.linalg.norm(NORTH)

        ENU = np.stack((EAST, NORTH, UP), axis=1) # E, N, U vecs are vertical, horizontal stacking

        return ENU
    
    def get_ENUvec_Earth(self, JDnow, y, frame='TOD'):
        """
        Internal function. Not vectorized in MRSguidance for performance.
        Calculates ENU vectors of provided statevec positions for Earth using WGS84, 
        in GCRF (and not ITRS!)
        Following use of these ENU vectors require data processing to be done
        in GCRF (e.g. for HA/FPA)

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TBD).
        y : array of floats
            State vector in the reference frame of planet.
        frame: string
            Frame (TOD or WGS84) in which ENU is determined. WGS84 considers
            the shape of Earth, TOD assumes circular object. 
            Default is TOD. WGS84 only usueful for initial launch direction.


        Returns
        -------
        ENU : 3x3 array of floats
            East, North, Up vectors / rotation matrices in GCRF.

        """
        
        if frame=='WGS84':

            t = ts.tdb_jd(JDnow)
            d = Distance(m=y[:3]) # input is nxm, n = 3, m = number of provided statevecs
            p = Geocentric(d.au, t=t) # p is in GCRF
            g = wgs84.geographic_position_of(p) # WGS84 position
    
            NEU = g.rotation_at(t)
            # NEU returns a 3D matrix:
            # axis 0: the three components of the originating system (North, East, Up)
            # axis 1: the n (three) dimensions of the target systemt
    
            # New matrix to bring into right order
            ENU = np.zeros((NEU.shape))
            ENU[0,:] = NEU[1,:] # first row becomes East
            ENU[1,:] = NEU[0,:] # seconds row becomes North
            ENU[2,:] = NEU[2,:] # last row remains Up
           
            # for return value, the rotation matrix is transposed     
            ENU = ENU.T
            
        # true equator and equinox of date (TETE, TOD)
        else:
            R = true_equator_and_equinox_of_date.rotation_at(ts.tdb_jd(JDnow))

            # position in TOD
            pos_TOD = R.dot(y[:3])
            
            # get local ENU
            ENU_TOD = self.get_ENUvec(JDnow, pos_TOD)
            
            # back to GCRF
            ENU = R.T.dot(ENU_TOD)
            
           
        return ENU

    def get_FPAHAVEL_EF(self, JDnow, y):
        """
        Internal function. Not vectorized in MRSguidance for performance.
        Returns velocity vector described w.r.t. to fixed Earth frame (EF), 
        as well as the FPA and heading w.r.t. to fixed Earth.

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

        # Earth rotation axis in ITRS
        earthRotAxisITRS = np.array([0,0,1])
        
        # Earth rotation axis in GCRF
        earthRotAxisGCRF = self.GCRF2ITRF.T.dot(earthRotAxisITRS)

        yEF = y * 1.0
        # remove velocity component of Earth
        yEF[3:] -= np.cross(EARTH_ROT_SPEED * earthRotAxisGCRF, yEF[:3])

        # get FPA / HA / ENU(GCRF)
        FPA, HA, ENU = self.get_FPAHA(JDnow, yEF, frame='Earth')

        # EF vel
        EFVEL = np.linalg.norm(yEF[3:])

        return FPA, HA, ENU, EFVEL

    def get_FPAHA(self, JDnow, y, frame='Earth'):
        """
        Internal function. Not vectorized in MRSguidance for performance.
        Returns FPA/HA in frame of provided statevectors. In case planet Earth
        is used, FPA/HA are calculated in the true local ENU frame in ITRS.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TDB).
        y : array of floats
            State vector in local frame (GCRF for Earth).
        frame : string, optional
            Frame in which FPA/HA are calculated. The default is 'Earth'.

        Returns
        -------
        FPA : float
            Flight Path Angle [rad] w.r.t. to provided frame.
        HA : float
            Heading Angle [rad] w.r.t. to provided frame.
            Measured from North to East (compas).

        """


        if frame == 'Earth':
            # get local ENU values for Earth
            ENU = self.get_ENUvec_Earth(JDnow, y)
        elif frame == 'launchsite':
            ENU = self.ENU_liftoff * 1.0
        elif frame == 'gcrf':
            ENU = np.eye(3)
        elif frame == 'vuw':
            ENU = self.get_VUW(JDnow, y)
        elif frame == 'vnb':
             ENU = self.get_VNB(JDnow, y)
        else:
            ENU = self.get_ENUvec(JDnow, y)
            

        # calc angle between velocity vec and local UP
        # - calc dot product between UP and velocity
        # - divise by product of their norms to get costheta
        # - adjust to 1/-1 if necessary
        # - calc angle

        # calc angle using the formula for vector dot multiplication (A.B = cos(theta)*norm(A)*norm(B))
        # ENU[:,:,2] = UP vectors (last entry in highest dimension, i.e. columns at the very right at all times)
        UPdotV = ENU[:,2].dot(y[3:])
        costheta = UPdotV / (np.linalg.norm(ENU[:,2]) * np.linalg.norm(y[3:]))

        # sometimes, cos values are out of -1/1; correction needed
        if costheta > 1.:
            costheta = 1.
        elif costheta < -1.:
            costheta = -1.

        # calc flight path angle [rad]
        FPA = np.arcsin(costheta)

        # projected vector
        vProjected = y[3:] - y[3:].dot(ENU[:,2]) * ENU[:,2]

        # save norm of projected vector, needed afterwards
        vProjectedNorm = np.linalg.norm(vProjected)

        # if no vProjectedNorm, then set vProjected to North vector
        if not vProjectedNorm:
            vProjected = ENU[:,1]

        # calc cos of vProjected w.r.t. NORTH
        costheta_north = ENU[:,1].dot(vProjected) / (np.linalg.norm(ENU[:,1])*vProjectedNorm)
        # calc cos of vProjected w.r.t. EAST
        costheta_east  = ENU[:,0].dot(vProjected) / (np.linalg.norm(ENU[:,0])*vProjectedNorm)

        # calculate azimuth angle; 0 angle is the north, 90° is east [rad]
        HA = np.arctan2(costheta_east, costheta_north)

        # correct HA to positive value
        if HA<0:
            HA += 2*np.pi

        return FPA, HA, ENU
    
    
    def get_VUW(self, JDnow, y):
        """
        Calculates and returns current VUW frame.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TDB).
        y : array of floats
            State vector in local frame (GCRF for Earth).

        Returns
        -------
        VUW : 3x3 array of flaots
            Rotation matrix for current VUW frame.

        """
        
        # e1: along velocity vector
        e1 = y[3:] * 1.0
        e1 /= np.linalg.norm(e1)
        
        # e3: along normal to orbital plane (ECI plane)
        e3 = np.cross(y[:3],e1)
        e3 /= np.linalg.norm(e3)
        
        # e2: e3 x e1
        e2 = np.cross(e3,e1)
        e2 /= np.linalg.norm(e2)
        
        VUW = np.stack((e1, e2, e3), axis=1)
        
        return VUW
    
    
    def get_VNB(self, JDnow, y):
        """
        Calculates and returns current VNB frame.

        Parameters
        ----------
        JDnow : float
            Current Julian Date (TDB).
        y : array of floats
            State vector in local frame (GCRF for Earth).

        Returns
        -------
        VNB : 3x3 array of flaots
            Rotation matrix for current VNB frame.

        """
        
        # e1 / velocity: along velocity vector
        e1 = y[3:] * 1.0
        e1 /= np.linalg.norm(e1)
        
        # e2 / normal: along normal to orbital plane
        e2 = np.cross(y[:3], e1)
        e2 /= np.linalg.norm(e2)
        
        # e3 / binormal: e1 x e2
        e3 = np.cross(e1,e2)
        e3 /= np.linalg.norm(e3)
    
        VNB = np.stack((e1, e2, e3), axis=1)
    
        return VNB
    
    
    
    
    
    