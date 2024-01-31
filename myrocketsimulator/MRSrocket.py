#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 12:11:24 2024

@author: thibault
"""

import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable 
from scipy import interpolate 

# parts table
I_partsName = 0
I_partsStaging = 1
I_partDryMass = 2
I_partsFuelMass = 3
I_partsDragArea = 4
I_partsRemFuel = 5



# engines table 
I_enginesName = 0
I_enginesDesc = 1
I_enginesThrustSL = 2
I_enginesThrustVAC = 3
I_enginesFueFlow = 4

# throttle table
I_throttleMET = 0
I_throttleStart = 1
I_throttleEnd = 2
I_throttleEngType = 3
I_throttleEngAmount = 4
I_throttleDesc = 5
I_throttleGrad = 6 # table filled in init_part()
I_validStart = 7 # table filled in fix_invalidThrottleMET()
I_invalidFlag = 8 # table filled in fix_invalidThrottleMET()
I_throttleActStage = 9 # table filled in add_stagingEventThrottle()
I_throttleDryMassOther = 10 # table filled in add_precalcValues()
I_throttleFuelMassStart = 11 # table filled in add_precalcValues()
I_throttleFuelMassEnd = 12 # table filled in add_precalcValues()
I_throttleFuelFlowConst = 13 # table filled in add_precalcValues()
I_throttleFuelFlowGrad = 14 # table filled in add_precalcValues()
I_throttleThrustSLallEngines = 15 # table filled in add_precalcValues()
I_throttleThrustVACallEngines = 16 # table filled in add_precalcValues()
I_throttleDragArea = 17 # table filled in add_precalcValues()
I_throttleDescPrinted = 18 # table filled in get_ThrustMassThrottleFlowrateRemFuel
    
class SpaceCraft():
    """
    The MRS SpaceCraft() class is used to generate an instance of spacecraft. A
    spacecraft can be built open different elements, called spacecraft elements (SCE)
    SCEs are defined by one or more parts,
    Each element is an instance of the SpaceCraftElement() class. It provides
    thrust, mass, drag area, drag coefficient and other values at any given 
    time of the flight (MET). 


    Attributes
    ----------
    spacecraftnanme: String
        Name of the spacecraft.
    mode: string
        Defines the spacecraft as 'active'.


    Methods
    -------
    get_ThrustMassOther()
        Returns the summed up thrust and mass values of all spacecraft elements.
    get_DragF()
        Returns the summed up drag force of all spacecraft elements.
    get_AreaCd()
        Returnst he summed up products of drag areas with drag coefficient.
    get_staticvalues()
        Returns a fixed set of static values; used for static spacraft propagation.
    set_fixedThrottlePointer()
        Sets a fixed throttle segment; used for step-wise propagation with integrator.
    reset_fixedThrottlePointer()
        Resets the fixed trottle segment pointer.
    get_EventsList()
        Returns a list of all segments (only availabel after simulation).
    plot_Cd()
        Plots the drag cofficient for a given spacraft element.
    plot_Thrust()
        Plots the VAC thrust for a given spacecraft element.

    """
    
    def __init__(self, scd):
        """
        Initializes the spacraft object.

        Parameters
        ----------
        scd : class
            Class that contains all relelvant spacraft data.

        Returns
        -------
        None.

        """
        
        self.name = scd.name # save spacecraft name
        self.SCD = scd # save spacraft data
        self.mode = 'active'
        
        # initialize spacecraft elements
        self.init_SCelements()
        
    def init_SCelements(self):
        """
        Initializes a list of objects for all spacecraft elements. 
        

        Returns
        -------
        None.

        """
        
        # make empty list for element objects
        self.SCelementsList = []
        
        # list with part names
        self.SCelementsNames = []

        # go through elements and initialize them as object and attach to list
        for i, SCelement in enumerate(self.SCD.SCelements):
            
            # check that amount is at least one
            if SCelement[1]>0:
                self.SCelementsList.append(SpaceCraftElement(eval('self.SCD.'+SCelement[0]), SCelement[1]))
                self.SCelementsNames.append(SCelement[0])
             
    
    def get_ThrustMassOther(self, MET, pa=0.0, verbose='v'):
        """
        Returns relevant summed up values for the spacraft. To be called by 
        simulator.

        Parameters
        ----------
        MET : float
            Mission Elapsed Time.
        pa : float, optional
            Ambient pressure. The default is 0.0.
        verbose : string, optional
            Set to v to print additional information while execution. The default is 'v'.

        Returns
        -------
        SCthrust : float
            Sum of all thrust values of all spacraft elements.
        SCmass : float
            Sum of all mass values of all spacraft elements.
        returnVals : array of floats
            All returned values from the spacraft element's object method 
            get_ThrustMassThrottleFlowrateRemFuel().

        """
       
        # initialize cummulated values
        SCthrust = 0
        SCmass = 0
        
        # memory
        returnVals = np.zeros((len(self.SCelementsList),6))
        
        for i, SCelement in enumerate(self.SCelementsList):
             # thrust, mass, throttle, flowrate, fuelPart, amount
             returnVals[i] = SCelement.get_ThrustMassThrottleFlowrateRemFuel(MET, pa, verbose)
             
             # add to cummulated values for amount of given part
             SCthrust += returnVals[i,0] * returnVals[i,5]
             SCmass += returnVals[i,1] * returnVals[i,5]
             
                 
        return SCthrust, SCmass, returnVals
    
    
    def get_DragF(self, MET, vrel=np.zeros((3)), rho=0.0, mach1=343.0):
        """
        Returns the cumulated drag vector of all spacecraft elements. To be used
        by simulator.

        Parameters
        ----------
        MET : float
            Mission Elapsed Time [s].
        vrel : array of floats OR single fluat
            Velocity vector (w.r.t. atmosphere). The default is np.zeros((3)).
        rho : float, optional
            Density [kg / m^3]. The default is 0.0.
        mach1 : float, optional
            Mach 1 [m/s]. The default is 343.0.

        Returns
        -------
        SCdragF : array of floats OR single float
            Cumulated drag force vector for all spacecraft elements.

        """
        
        # initialize cummulated values
        if isinstance(vrel, float):
            SCdragF = 0
        else:
            SCdragF = np.zeros((3))
        
        #for SCelement in self.SCelementsList:
        for i, SCelement in enumerate(self.SCelementsList):
             # drag
             SCdragF += SCelement.get_DragF(MET, vrel, rho, mach1) * SCelement.SCED.amount
         
        return SCdragF

    
    def get_AreaCd(self, MET, mach=343.0):
        """
        Returns the cumulated product of drag area and drag cofficients of all
        spacecraft elements. 

        Parameters
        ----------
        MET : float
            Mission Elapsed Time [s].
        mach : float, optional
            Actual Mach speed. The default is 343.0.

        Returns
        -------
        AreaCd : float
            Sum of all area-cofficient products.

        """
        
        # returns sum of area * Cd values (needed for external drag computation)
    
        # initialize cummulated values
        AreaCd = 0
        
        #for SCelement in self.SCelementsList:
        for i, SCelement in enumerate(self.SCelementsList):
             # drag
             AreaCd += SCelement.get_Cd(MET, mach) * SCelement.get_DragArea(MET) * SCelement.SCED.amount
         
        return AreaCd

    def get_staticvalues(self):
        """
        Returns separate static values

        Returns
        -------
        Float
            Static mass.
        Float
            Static drag area.
        Float
            Static drag cofficient Cd.
        Float
            Static reflectivity coefficient Cr.
        Float
            Static SRP area.

        """
        
        # returns static values
        return self.SCD.staticValues.mass, \
               self.SCD.staticValues.dragarea, \
               self.SCD.staticValues.Cd, \
               self.SCD.staticValues.Cr, \
               self.SCD.staticValues.SRParea
        
        
    
    def set_fixedThrottlePointer(self, MET):
        """
        Sets fixed throttle pointers for all spacecraft elements. 

        Parameters
        ----------
        MET : float
            Mission Elapsed Time [s].

        Returns
        -------
        int
            Always returns 1.

        """
        
        # set fixed throttle pointer in all elements
        for SCelement in self.SCelementsList:
            SCelement.set_fixedThrottlePointer(MET)
            
        return 1
     
    def reset_fixedThrottlePointer(self):
        """
        Deactivates the fixed throttle pointer for all spacecraft elements.

        Returns
        -------
        int
            Always returns 1.

        """
        
        # set fixed throttle pointer in all elements
        for SCelement in self.SCelementsList:
            SCelement.fixedThrottlePointer = -1
        
        return 1
    
    def get_EventsList(self):
        """
        Gets the event lists from all spacecraft elements and sorts them by MET.

        Returns
        -------
        Array
            MET and description of events.

        """
        
        # make empty event dataframe; two columns for MET and string
        self.eventsDF = np.empty((0,2), dtype=object)
        
        # returns event list
        for SCelement in self.SCelementsList:
            self.eventsDF = np.vstack((self.eventsDF, SCelement.eventsDF))
            
        # sorts by MET
        self.eventsDF = self.eventsDF[self.eventsDF[:, 0].argsort()]
        
        # round MET values
        self.eventsDF[:,0] = np.round(self.eventsDF[:,0].astype(float),3)
        
        return self.eventsDF

    def get_METvalues(self):
        """
        Returns an array of all MET values of all element's throttleTable.

        Returns
        -------
        METvalues : array of floats
            MET values of all spacecraft elements.

        """
        
        # returns event list
        METvalues = np.array([])
        for SCelement in self.SCelementsList:
            METvalues = np.append(METvalues, SCelement.throttleTable[:,0])
            
        # remove duplicates and sort
        METvalues = np.unique(METvalues)
        
        return METvalues

    def plot_Cd(self, SCEname): 
        """
        Plots the drag coefficient for the given spacecraft element/part.
        

        Parameters
        ----------
        SCEname : string
            Name of the spacecraft element.

        Returns
        -------
        figCd : figure handle
            Fig handle returned for further plot manipulations.
        axCd : axis handle
            Axis handle returned for further plot manipulations.

        """
        
        # get position of part in SCelementsList
        partIndex = self.SCelementsNames.index(SCEname)
        
        # plot Cd
        figCd, axCd = self.SCelementsList[partIndex].plot_Cd()
        
        return figCd, axCd
    
    def plot_Thrust(self, SCEname): 
        """
        Plots the VAC thrust profile for the given spacecraft element.

        Parameters
        ----------
        SCEname : string
            Name of the spacecraft element.

        Returns
        -------
        figCd : figure handle
            Fig handle returned for further plot manipulations.
        axCd : axis handle
            Axis handle returned for further plot manipulations.

        """
        
        # get position of part in SCelementsList
        partIndex = self.SCelementsNames.index(SCEname)
        
        # plot Cd
        figThrust, axThrust = self.SCelementsList[partIndex].plot_Thrust()
        
        return figThrust, axThrust



class SpaceCraftElement():
    """
    The SpaceCraftElement() class is used to generate a element of spacecraft. 
    Each element is an instance of the SpaceCraftElement() class. It provides
    thrust, mass, drag area, drag coefficient and other values at any given 
    time of the flight (MET). Different elements of a spacecraft can be:
        - the main stages
        - boosters
        - payload
        - fairings, escape towers, ullage engines, ...


    Attributes
    ----------
    amount: int
        How often this element is part of the spaceacraft (e.g. 2 for SLS boosters)
    fixedThrottlePointer: int
        Currently set pointer for the throttle table
        
    Methods
    -------
    set_fixedThrottlePointer()
        Sets the pointer to the current throttle profile.
    get_ThrustMassThrottleFlowrateRemFuel()
        Returns thrust and mass, as well as other values at given time.
    get_DragF()
        Returns the drag force at given time for given vrel vector.
    get_Cd()
        Returns current drag coefficient (also used by get_DragF()).
    get_DragArea()
        Returns current drag area (also used by get_DragF()).
    plot_Cd()
        Plots drag coefficient.
    plot_Thrust()
        Plots thrust in VAC.  

    """
    
    def __init__(self, sced, amount=1):
        """
        Initialization.

        Parameters
        ----------
        sced : class
            Spacecraft element data.
        amount : int, optional
            Amount of spacecraft elements of this kind. The default is 1.

        Returns
        -------
        None.

        """
        
        # save spacraft element data
        self.SCED = sced  

        # save amount of this kind of element (e.g. for 2 boosters, 20 payloads, ...)         
        self.SCED.amount = amount
        
        # initialize fixedThrottlePointer
        self.fixedThrottlePointer = -1
        
        # initialize part
        self.init_SCE()
        
        # initialize Cd interpolation
        self.init_Cdinterp()
        
        # make empty event dataframe; two columns for MET and string
        self.eventsDF = np.empty((0,2), dtype=object)
        
    def init_SCE(self):
        """
        Performs an initialization of the spacecraft element. A major part is to 
        generate a throttle table that also includes the mass and other properties
        for any time during the flight. Also,  invalid MET values are corrected.

        Returns
        -------
        None.

        """
        
        # allowed time step
        self.valid_stepsize = 0.1

        # add empty table if no profile table was provided
        if not hasattr(self.SCED, 'throttleInit'): 
            self.SCED.throttleInit = np.array([
                 # MET    Start   End     EngineType  EngineAmount,   Description
                 [ -86400,   0.0,  0.0,     -1,          0,            ''],
                 ], dtype=object) 

        # make a copy of throttleInit table to work with
        self.throttleTable = self.SCED.throttleInit * 1
        
        # when throttle segment End-value == -1, replace by following rule start value
        for i in range(len(self.throttleTable)-1):
            if self.throttleTable[i,I_throttleEnd] == -1:
                self.throttleTable[i,I_throttleEnd] = self.throttleTable[i+1,I_throttleStart]
                
        # copy back to init table (only modified end values)
        self.SCED.throttleInit = self.throttleTable * 1
                
        # round MET times to milliseconds
        self.throttleTable[:,I_throttleMET] = np.round(self.throttleTable[:,I_throttleMET].astype(float),3)

        # last throttle segment has always 0 thrust
        self.throttleTable[-1,I_throttleStart] = 0.0
        self.throttleTable[-1,I_throttleEnd] = 0.0

        # insert gradients 
        gradients = (self.throttleTable[:-1,I_throttleEnd] - self.throttleTable[:-1,I_throttleStart]) /\
                      (self.throttleTable[1:,I_throttleMET] - self.throttleTable[:-1,I_throttleMET])
        self.throttleTable = np.hstack((self.throttleTable, np.append(gradients,0)[:,None]))

        # fix throttle table for invalid MET times 
        self.fix_invalidThrottleMET()
        
        # add staging events to throttle table
        self.add_stagingEventThrottle()
        
        # add precalculated values for fuel consumption and overall rocket specs
        self.add_precalcValues()
        
        # add remaining fuel after staging to parts table
        self.add_remFuel()
        
        # add flag if throttle segment desription was alreay shown
        descPrinted = np.zeros((len(self.throttleTable),1), dtype=int) 
        self.throttleTable = np.append(self.throttleTable, descPrinted, axis=1)


    def fix_invalidThrottleMET(self):
        """
        Invalid MET values (i.e. smaller than a 1/10 of a second) get compensated
        by using intermediate throttle profiles resulting in the same impulse [Ns]
        the with invalid MET values.

        Returns
        -------
        None.

        """
        
        # insert allowed start time and mark if actual start time is not allowed
        starttime = np.zeros((len(self.throttleTable),1))     
        invalidstarttime = np.zeros((len(self.throttleTable),1)) 
        
        for i in range(len(self.throttleTable)):
            starttime[i] = np.round(np.floor(np.round(self.throttleTable[i,I_throttleMET]\
                                                      /self.valid_stepsize,2))*self.valid_stepsize,2)
            if np.round(self.throttleTable[i,I_throttleMET]-starttime[i],3) != 0:
                invalidstarttime[i] = 1
        self.throttleTable = np.hstack((self.throttleTable, starttime, invalidstarttime))

        # correct invalid MET values by using new throttle segments with specific 
        # fixed throttle values at allowed MET and with the same impulse [Ns]
        i = 0
        while True:
            
            # if in last profile, break
            if i==len(self.throttleTable)-1:
                 break
            
            # skip first throttle segment 
            i += 1
                
            # skip is start time is ok
            if not self.throttleTable[i,I_invalidFlag]:
                continue

            # calc throttle secs for valid time step
            METrange = [self.throttleTable[i,I_validStart], self.throttleTable[i,I_validStart] + self.valid_stepsize ]
            valid_throttle_sec = self.get_throttleSecondsRange(METrange)

            
            # previous segment always gets updated
            self.throttleTable[i-1,I_throttleEnd] = self.get_throttle(METrange[0])
            
            # index values for throttle settings with same allowed start time
            index_nonvalid = np.logical_and((self.throttleTable[i,I_validStart]-self.throttleTable[:,I_validStart])==0,\
                                            self.throttleTable[:,I_invalidFlag]==1).nonzero()[0]
                
            print('Fixing non-valid MET times for throttle segment(s) of %s:' % (self.SCED.name))
            for j in range(len(index_nonvalid)):
                print('\tSegment %i, MET=%.5f ' % (index_nonvalid[j], self.throttleTable[index_nonvalid[j],I_throttleMET]  ))
            
            # new start value for last illegal segment
            self.throttleTable[index_nonvalid[-1],I_throttleStart] = self.get_throttle(METrange[1])
            # adjust start time for last illegal segment
            self.throttleTable[index_nonvalid[-1],I_throttleMET] = METrange[1]
            # reset invalid flag for last illegal segment
            self.throttleTable[index_nonvalid[-1],I_invalidFlag] = 0   
            
            # make new throttle segment
            throttle_new = np.empty((1,self.throttleTable.shape[1]),dtype=object)
            throttle_new[0,I_throttleMET] = METrange[0]
            throttle_new[0,I_throttleStart] = valid_throttle_sec/self.valid_stepsize
            throttle_new[0,I_throttleEnd] =   valid_throttle_sec/self.valid_stepsize
            throttle_new[0,I_throttleEngType] = self.throttleTable[i,I_throttleEngType] 
            throttle_new[0,I_throttleEngAmount] = self.throttleTable[i,I_throttleEngAmount] 
            throttle_new[0,I_throttleDesc] = 'MET-adjusted throttle segment.'
            throttle_new[0,I_throttleGrad] = 0.
            throttle_new[0,I_throttleGrad:] = 0.
            
            # insert new throttle segment
            self.throttleTable = np.insert(self.throttleTable,index_nonvalid[-1],throttle_new, axis=0)
         
            # delete all but the last illegal segments 
            self.throttleTable = np.delete(self.throttleTable, index_nonvalid[:-1], axis=0)
            
        # round MET times to milliseconds
        self.throttleTable[:,I_throttleMET] = np.round(self.throttleTable[:,I_throttleMET].astype(float),3)
          
        # delete throttle segments with identical MET start times (but keep the last)
        identicalMET = (self.throttleTable[1:,I_throttleMET] - self.throttleTable[:-1,I_throttleMET]) == 0
        self.throttleTable = np.delete(self.throttleTable, np.append(identicalMET, False), axis=0)

        # delete contents of now useless columns
        self.throttleTable[:,I_validStart] = -1
        self.throttleTable[:,I_invalidFlag] = -1
        
        return None
       
    def add_stagingEventThrottle(self):
        """
        Adds staging events to the throttle table, needed to accurately calculate
        the mass of the spacecraft and provide an updated drag area value.

        Returns
        -------
        None.

        """
                
        # add column for active stage
        actStage = np.zeros((len(self.throttleTable),1), dtype=int)   
        self.throttleTable = np.append(self.throttleTable, actStage, axis=1)
        
        
        # add staging
        for i in range(len(self.SCED.partsInit)):
            
            # get staging time to 1/10 of second 
            stagingTime = round(self.SCED.partsInit[i,I_partsStaging],1)
            
            # warning if provided staging time was not valid
            if stagingTime != self.SCED.partsInit[i,I_partsStaging]:
                print('Fixing non-valid MET time for staging of %s to %.2f seconds.' % \
                      (self.SCED.partsInit[i,I_partsName], stagingTime ))
            
            # find throttle segment starting before staging
            throttleIDpriorStaging = ((stagingTime - self.throttleTable[:,I_throttleMET])>=0).nonzero()[0][-1]
         
            # check if intermediate throttle segment is needed
            # if staging time is after MET of previous throttle segment --> new
            if stagingTime > self.throttleTable[throttleIDpriorStaging,I_throttleMET]:
                
                # update end value of previous throttle segment
                self.throttleTable[throttleIDpriorStaging,I_throttleEnd] = self.get_throttle(stagingTime)
               
                # make new segment for staging
                throttle_new = np.empty((1,self.throttleTable.shape[1]),dtype=object)
                throttle_new[0,I_throttleMET] = stagingTime
                throttle_new[0,I_throttleStart] = 0
                throttle_new[0,I_throttleEnd] =   0
                throttle_new[0,I_throttleEngType] = -1
                throttle_new[0,I_throttleEngAmount] = 0
                throttle_new[0,I_throttleDesc] = 'Staging of %s.' % (self.SCED.partsInit[i,I_partsName])
                throttle_new[0,I_throttleGrad] = 0.
                throttle_new[0,I_throttleGrad:] = 0.
                
                # insert new throttle segment
                self.throttleTable = np.insert(self.throttleTable,throttleIDpriorStaging+1,throttle_new, axis=0)
                
                # set new active stage
                self.throttleTable[throttleIDpriorStaging+1:,I_throttleActStage] = i+1
                
            elif stagingTime == self.throttleTable[throttleIDpriorStaging,I_throttleMET]:
                
                # add staging info to event 
                self.throttleTable[throttleIDpriorStaging,I_throttleDesc] = \
                    'Staging of %s. ' % (self.SCED.partsInit[i,I_partsName]) + \
                        self.throttleTable[throttleIDpriorStaging,I_throttleDesc]
                
                # set active stage
                self.throttleTable[throttleIDpriorStaging:,I_throttleActStage] = i+1
                
        return None
                
    def add_precalcValues(self):
        """
        Adds all kind of values (masses, drag area, fuel flow, ...) to the 
        throttle table.

        Returns
        -------
        None.

        """
        
        # add columns for new values
        self.throttleTable = np.append(self.throttleTable, np.zeros((len(self.throttleTable),8)), axis=1) 

        previousStage = -1
        for i in range(len(self.throttleTable)-1):
            
            # get current stage for current segment
            actStage = self.throttleTable[i, I_throttleActStage]
            
            # calc dry masses + fuel of non active stages for current throttle segment
            mass_dryNonActiveFuel = np.sum(self.SCED.partsInit[actStage:,I_partDryMass]) + \
                                    np.sum(self.SCED.partsInit[actStage+1:,I_partsFuelMass]) 
                                                            
            # calc fuel at beginning of throttle segment  
            # if new stage
            if actStage != previousStage:
                mass_fuelStart = self.SCED.partsInit[actStage, I_partsFuelMass]
                # set previousStage to current stage so that stating is detected 
                previousStage = actStage
            # same stage
            else:
                mass_fuelStart = self.throttleTable[i-1, I_throttleFuelMassEnd]
                
            # if no fuel left
            if mass_fuelStart <= 0:
                # set throttle to 0 
                self.throttleTable[i, I_throttleStart] = 0.
                self.throttleTable[i, I_throttleEnd] = 0.
                self.throttleTable[i, I_throttleGrad] = 0.
            
            # calc only if engine is provided and fuel available
            if self.throttleTable[i, I_throttleEngType] >=0 and mass_fuelStart > 0:
                    
                # fuel consumption (constant)
                fuelFlow_const = - self.SCED.engines[self.throttleTable[i, I_throttleEngType],I_enginesFueFlow] * \
                                        self.throttleTable[i, I_throttleEngAmount] * \
                                        self.throttleTable[i, I_throttleStart]
                                        
                # fuel consumption (gradient)
                fuelFlow_grad = - self.SCED.engines[self.throttleTable[i, I_throttleEngType],I_enginesFueFlow] * \
                                        self.throttleTable[i, I_throttleEngAmount] * \
                                        self.throttleTable[i, I_throttleGrad] * .5
                                        
                # pre-calculate SL and VAC thrust at 100% 
                thrustSLallEngines  = self.throttleTable[i, I_throttleEngAmount] * \
                                      self.SCED.engines[self.throttleTable[i, I_throttleEngType], I_enginesThrustSL]
                thrustVACallEngines = self.throttleTable[i, I_throttleEngAmount] * \
                                      self.SCED.engines[self.throttleTable[i, I_throttleEngType], I_enginesThrustVAC]
                                        
            else:
                fuelFlow_const = 0.
                fuelFlow_grad = 0.
                thrustSLallEngines = 0.
                thrustVACallEngines = 0.
            
            # calc remaining fuel at end of segment
            segDuration = self.throttleTable[i+1, I_throttleMET] - self.throttleTable[i, I_throttleMET]
            mass_fuelEnd = mass_fuelStart + fuelFlow_const * segDuration + \
                                            fuelFlow_grad * segDuration**2
                                            
            # check if fuel is left
            if mass_fuelEnd <= 0 and mass_fuelStart > 0 :
                mass_fuelEnd = 0.
                print('WARNING: stage %s runs out of full in throttle segment %i.' % \
                      (self.SCED.partsInit[actStage,I_partsName], i ))
                
            
            # update throttleTable
            self.throttleTable[i, I_throttleDryMassOther] = mass_dryNonActiveFuel
            self.throttleTable[i, I_throttleFuelMassStart] = mass_fuelStart
            self.throttleTable[i, I_throttleFuelMassEnd] = mass_fuelEnd
            self.throttleTable[i, I_throttleFuelFlowConst] = fuelFlow_const
            self.throttleTable[i, I_throttleFuelFlowGrad] = fuelFlow_grad
            self.throttleTable[i, I_throttleThrustSLallEngines] = thrustSLallEngines
            self.throttleTable[i, I_throttleThrustVACallEngines] = thrustVACallEngines
            self.throttleTable[i, I_throttleDragArea] = self.SCED.partsInit[actStage, I_partsDragArea]
            
            
        return None   
            

    def add_remFuel(self):
        """
        Adds remaining fuel to list of parts after they detached.

        Returns
        -------
        None.

        """
        
        # add column for remaining fuel in stage-table; column index: I_partsRemFuel
        remFuel = np.zeros((len(self.SCED.partsInit),1))
        self.SCED.partsInit = np.append(self.SCED.partsInit, remFuel, axis=1)
        
        previousPart= 0
        for i in range(len(self.throttleTable)):
            # detect stating 
            if self.throttleTable[i, I_throttleActStage] != previousPart: 
                # store remaining fuel in parts table
                self.SCED.partsInit[previousPart, I_partsRemFuel] = self.throttleTable[i-1, I_throttleFuelMassEnd]
                # update part id
                previousPart = self.throttleTable[i, I_throttleActStage] 
                
    
    def get_throttle(self, MET):
        """
        Returns current throttle value.

        Parameters
        ----------
        MET : float
            Mission Elapsed Time [s].

        Returns
        -------
        throttleCurrent : float
            Throttle value (0.0-1.).

        """
        # pointer to current row
        tP = ((MET - self.throttleTable[:,I_throttleMET])>=0).nonzero()[0][-1]
        
        # current throttle value
        throttleCurrent = (MET-self.throttleTable[tP,I_throttleMET]) * self.throttleTable[tP,I_throttleGrad] +\
                self.throttleTable[tP,I_throttleStart]
        
        return throttleCurrent
    
    def get_throttleSecondsSegment(self, segmentID):
        """
        Returns the throttle-value intergrated over time for a given throttle segment.

        Parameters
        ----------
        segmentID : int
            Throttle segment ID.

        Returns
        -------
        throttle_sec : float
            Throttle-value intergrated over time.

        """
        
        throttle_sec = (self.throttleTable[segmentID, I_throttleStart] + self.throttleTable[segmentID, I_throttleEnd])/2 * \
                       (self.throttleTable[segmentID+1, I_throttleMET] - self.throttleTable[segmentID, I_throttleMET])
        
        return throttle_sec
    
    
    def get_throttleSecondsRange(self, METrange):
        """
        Returns the throttle-value intergrated over time for a given MET range.

        Parameters
        ----------
        METrange : array of floats
            MET start and end values.

        Returns
        -------
        throttle_sec : float
            Throttle-value intergrated over time..

        """
        
        # find start and end settings rows
        tPstart = ((METrange[0] - self.throttleTable[:,I_throttleMET])>=0).nonzero()[0][-1]
        tPend = ((METrange[1] - self.throttleTable[:,I_throttleMET])>=0).nonzero()[0][-1]
    
        
        # number of relevant rows
        tPcount = tPend-tPstart+1
        
        # all happens in one row
        if tPstart == tPend:
            throttle_start = self.get_throttle(METrange[0])
            throttle_end = self.get_throttle(METrange[1])
            throttle_sec = (METrange[1]-METrange[0]) * (throttle_start+throttle_end)/2
            
        # two rows involved
        elif tPend > tPstart:
            # first/left settings
            throttle_start_first = (METrange[0]-self.throttleTable[tPstart,I_throttleMET]) \
                * self.throttleTable[tPstart,I_throttleGrad] + self.throttleTable[tPstart,I_throttleStart]
            throttle_end_first =  self.throttleTable[tPstart, I_throttleEnd]
            duration_first =  self.throttleTable[tPstart+1, I_throttleMET] - METrange[0]
            throttle_sec_first = duration_first * (throttle_start_first+throttle_end_first)/2
            
            # last/right settings
            throttle_start_last = self.throttleTable[tPend, I_throttleStart]
            throttle_end_last =  (METrange[1]-self.throttleTable[tPend,I_throttleMET]) \
                * self.throttleTable[tPend,I_throttleGrad] + self.throttleTable[tPend,I_throttleStart]
            duration_last =  METrange[1] - self.throttleTable[tPend, I_throttleMET] 
            throttle_sec_last = duration_last * (throttle_start_last+throttle_end_last)/2
               
            # add the two throttle segments
            throttle_sec = throttle_sec_first + throttle_sec_last
            
            # add the throttle segments inbetween
            for i in range(tPstart+1,tPend):
               throttle_sec += self.get_throttleSecondsSegment(i)
        
        else:
            print('ERROR.')
            throttle_sec = -1
       
        return throttle_sec     
 
    def set_fixedThrottlePointer(self, MET):
        """
        Sets the trhottle segment pointer.

        Parameters
        ----------
        MET : float
            Mission Elapsed Time [s].

        Returns
        -------
        None.

        """
        
        
        self.fixedThrottlePointer = ((MET - self.throttleTable[:,I_throttleMET])>=0).nonzero()[0][-1]
        
        return None
        
    def get_ThrustMassThrottleFlowrateRemFuel(self, MET, pa=0.0, verbose=''):
        """
        Returns all relevant spacecraft values at given time. Pressure is provided
        for thrust calculations (SL/VAC).        
        
        Parameters
        ----------
        MET : float
            Mission Elapsed Time [s].
        pa : float, optional
            Pressure in Pascal. The default is 0.0.
        verbose : string, optional
            Set to 'v' for additional output while execution. The default is ''.

        Returns
        -------
        thrust : float
            Current thrust [N].
        mass : float
            Current mass of whole SCE [kg].
        throttle : float
            Current throttle [%].
        flowrate : float
            Fuel flow rate [kg/s].
        fuelPart : float
            How much fuel the current part has left.
        int
            Amount of instaleld spacecraft elements of this kind.

        """
        
        # get throttle pointer (tP) in throttle profile if not preset 
        if self.fixedThrottlePointer == -1:
            tP =  ((MET - self.throttleTable[:,I_throttleMET])>=0).nonzero()[0][-1]
        # otherwhise tP should be already set
        else:
            tP = self.fixedThrottlePointer 
            
        # print + save throttle segment description 
        if not self.throttleTable[tP, I_throttleDescPrinted] and \
            self.throttleTable[tP, I_throttleDesc] != '':
            
            # get current stage name for current segment
            actStage = self.throttleTable[tP, I_throttleActStage]    
            
            # make sure its a valid stage
            if actStage < len(self.SCED.partsInit):
                partName = self.SCED.name
            else:
                partName = self.SCED.name + ' terminated'
            
            # define message
            LogMessage = partName + ': ' + self.throttleTable[tP, I_throttleDesc]
            
            # show only if required
            if verbose=='v':
                print('%.2f:\t\t%s' % (MET, LogMessage))
            
            # set description as printed
            self.throttleTable[tP, I_throttleDescPrinted] = 1 
            
            # add event to events dataframe of spacecraft element
            evententry = np.array([MET, LogMessage], dtype=object)
            self.eventsDF = np.vstack((self.eventsDF, evententry))
            
                
        # calc time in current profile
        timeInProfile = MET - self.throttleTable[tP,I_throttleMET]
        
        # calc current throttle value
        throttle = self.throttleTable[tP, I_throttleStart] +\
                   self.throttleTable[tP, I_throttleGrad] * timeInProfile
        # calc thrust (for SL and VAC)
        thrustSL = throttle * self.throttleTable[tP, I_throttleThrustSLallEngines] 
        thrustVAC = throttle * self.throttleTable[tP, I_throttleThrustVACallEngines] 
        
        # calc thrust of spacecraft element
        # 1. Calc gradient for Thrust w.r.t to atmospheric pressure)
        thrustVSpressure = (thrustSL-thrustVAC) / 101325.0
        # 2. Calc thrust for given pressure
        thrust = thrustVAC + thrustVSpressure * pa
        
        # calc mass of spacecraft element
        # 1. Calc mass of fuel of active part
        fuelPart = self.throttleTable[tP, I_throttleFuelMassStart] +\
                   self.throttleTable[tP, I_throttleFuelFlowConst] * timeInProfile + \
                   self.throttleTable[tP, I_throttleFuelFlowGrad] * timeInProfile**2
        # 2. Add mass of other masses (dry mass of current part + complete other parts)
        mass = fuelPart + self.throttleTable[tP, I_throttleDryMassOther]
        
        # calc flowrate
        flowrate = - (self.throttleTable[tP, I_throttleFuelFlowConst] +\
                    2*self.throttleTable[tP, I_throttleFuelFlowGrad] * timeInProfile)
        
        # reset values if no fuel is left
        if fuelPart <=0:
            fuelPart = 0 
            thrust = 0
            flowrate = 0
       
        return thrust, mass, throttle, flowrate, fuelPart, self.SCED.amount
        
        
    
    def get_DragF(self, MET, vrel=np.zeros((3)), rho=0.0, mach1=343.0, AoA=0.0):
        """
        Calculates current drag force vector.

        Parameters
        ----------
        MET : float
            Mission Elapsed Time [s].
        vrel : array of floats 
            Velocity vector (w.r.t. atmosphere). The default is np.zeros((3)).
        rho : float, optional
            Density of local atmosphere [kg/m^3]. The default is 0.0.
        mach1 : float, optional
            Mach 1 speed [m/s]. The default is 343.0.
        AoA : float, optional
            Angle of attack [°]. Not implemented yet The default is 0.0.

        Returns
        -------
        DragF : float
            Drag force vector.

        """
        
        # norm of vrel vector
        vrel_norm = np.linalg.norm(vrel)
        
        # calc drag force for given velocity, density and mach number
        DragF = - 0.5 * vrel * vrel_norm * self.get_Cd(MET, vrel_norm/mach1, AoA)\
                * self.get_DragArea(MET) * rho
        
        return DragF
       
    
    def init_Cdinterp(self):
        """
        Initializes the scipy interp1d-interpolation object.

        Returns
        -------
        int
            Always returns 1.

        """
        
        # if no Cd table provided
        if not hasattr(self.SCED, 'C_D'): 
            self.CDinterpActive = 0
        # if Cd table provided
        else:
            self.CDinterpActive = 1
            self.CDinterp = interpolate.interp1d(self.SCED.C_D[:,0], self.SCED.C_D[:,1], kind='cubic')
            
        return 1
        
    
    def get_Cd(self, MET, mach, AoA=0.0):
        """
        Returns drag coefficient for single values or arrays of Mach values.

        Parameters
        ----------
        MET : float
            Mission Elapsed Time [s].
        mach : float or array of floats
            Mach alue.
        AoA : float, optional
            Angle of attack. Not implemented yet. The default is 0.0.

        Returns
        -------
        Float or array of floats
            Drag coefficient(s) for given Mach value(s).

        """
        
        # return 0 if no Cd interpolation not available
        if self.CDinterpActive == 0:
            return 0
        
        # if faster than highest mach value, provide last C_D value (only single values)
        if np.size(mach) == 1:
            if mach>self.SCED.C_D[-1,0]:
                return self.SCED.C_D[-1,1]
        
        return self.CDinterp(mach)
    
    def plot_Cd(self):
        """
        Plots drag cofficient.

        Returns
        -------
        figCd : figure handle
            Fig handle returned for further plot manipulations.
        axCd : axis handle
            Axis handle returned for further plot manipulations.

        """
        
        # return 0 if no Cd interpolation not available
        if self.CDinterpActive == 0:
            return 0
        
        # get min/max values for Mach number
        MachMinMax = np.array([self.SCED.C_D[0,0], self.SCED.C_D[-1,0]])
        
        # make vector of Mach values
        MachRange = np.linspace(MachMinMax[0], MachMinMax[1], 1000)
        
        # get Cd Values (hard coded for AoA=0 at the moment)
        Cd = self.get_Cd(0, MachRange)
        
        # make plot
        figCd, axCd = plt.subplots(1,1)
        axCd.set_xlabel('Mach Number, M')
        axCd.set_ylabel('Drag Coefficient, Cd')
        axCd.grid()
        axCd.set_title('Drag coefficient Cd for %s' % (self.SCED.name))
        axCd.set_xticks(np.linspace(0,np.ceil(MachMinMax[1]),int(np.ceil(MachMinMax[1]+1))))
        
        axCd.plot(MachRange, Cd, 'b')
        axCd.plot(self.SCED.C_D[:,0], self.SCED.C_D[:,1], 'rx')
        
        return figCd, axCd
        
    
    def get_DragArea(self, MET=0):
        """
        Returns current drag area.

        Parameters
        ----------
        MET : float
            Mission Elapsed Time [s]. Default is 0.

        Returns
        -------
        Float
            Drag area.

        """
        
        
        # get throttle pointer (tP) in throttle profile if not preset 
        if self.fixedThrottlePointer == -1:
            tP =  ((MET - self.throttleTable[:,I_throttleMET])>=0).nonzero()[0][-1]
        # otherwhise tP should be already set
        else:
            tP = self.fixedThrottlePointer 
            
        # return drag area
        return self.throttleTable[tP,I_throttleDragArea]
        
    
    def get_totalImpulse(self):
        """
        Returns total impulse of SCE in VAC in [Ns]. Returns faulty values if
        stage runs out of fuel within a throttle segment.

        Returns
        -------
        impulse : Float
            Total impulse of spacecraft element.

        """
        
        # return the total VAC impulse (Ns) for spacecraft element
        # Faulty values if stage runs out of fuel within a segment
        
        # calculate total impulse (by summing up integrals of trapezoids)
        impulse = ((self.throttleTable[:-1, I_throttleStart] + self.throttleTable[:-1, I_throttleEnd])/2 * \
                       (self.throttleTable[1:, I_throttleMET] - self.throttleTable[:-1, I_throttleMET]) * \
                       self.throttleTable[:-1, I_throttleThrustVACallEngines]).sum()

        return impulse
        
        
    def plot_Thrust(self):
        """
        Plots VAC thrust of SCE.

        Returns
        -------
        figCd : figure handle
            Fig handle returned for further plot manipulations.
        axCd : axis handle
            Axis handle returned for further plot manipulations.

        """
        
        # prepare plot window
        figThrust, axThrust = plt.subplots(1,1)
        axThrust.set_xlabel('MET [s]')
        axThrust.set_ylabel('Thrust [kN]')
        axThrust.grid()
        axThrust.set_title('VAC thrust profile for %s' % (self.SCED.name))
       
        # loop through throttle segments
        for i in range(len(self.throttleTable)-1):
            # calc thrust values in segment [kN]
            thrustStart = self.throttleTable[i,I_throttleThrustVACallEngines] * \
                          self.throttleTable[i,I_throttleStart]/1000
            thrustEnd =  self.throttleTable[i,I_throttleThrustVACallEngines] * \
                         self.throttleTable[i,I_throttleEnd]/1000
            thrustStartFollowing = self.throttleTable[i+1,I_throttleThrustVACallEngines] * \
                          self.throttleTable[i+1,I_throttleStart]/1000
                                         
            
            # plot throttle segment
            axThrust.plot([self.throttleTable[i,I_throttleMET], self.throttleTable[i+1,I_throttleMET]], \
                          [thrustStart, thrustEnd],'b')
            
            # if end value of previous and start value of following differ:
            if thrustEnd != thrustStartFollowing:
                axThrust.plot([self.throttleTable[i+1,I_throttleMET], self.throttleTable[i+1,I_throttleMET]], \
                         [thrustEnd, thrustStartFollowing], 'b')
        
    
        return figThrust, axThrust