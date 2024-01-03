#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:08:23 2023

@author: thibault
"""

import numpy as np
import pandas as pd
from decimal import Decimal # used for precise modulo calculations
import matplotlib.pyplot as plt
from scipy import interpolate 



class SpaceCraft():
    
    def __init__(self, scd):
        
        self.spacecraftname = scd.name # save spacecraft name
        self.SCD = scd # save spacraft data
        self.mode = 'active'
        
        # initialize spacecraft elements
        self.init_SCelements()
        
    def init_SCelements(self):
        
        # make empty list for element objects
        self.SCelementsList = []

        # go through elements and initialize them as object and attach to list
        for i, SCelement in enumerate(self.SCD.SCelements.values):
            
            # check that amount is at least one
            if SCelement[1]>0:
                self.SCelementsList.append(SpaceCraftElement(eval('self.SCD.'+SCelement[0]), SCelement[1]))
                

    def get_ThrustMassOther(self, MET, pa=0.0, verbose='v'):
       
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
    
    def get_DragF(self, MET, vrel=0.0, rho=0.0, mach1=0.0):
        # returns drag force in [N], as a direction vector
        
        # initialize cummulated values
        SCdragF = np.zeros((3))
        
        #for SCelement in self.SCelementsList:
        for i, SCelement in enumerate(self.SCelementsList):
             # drag
             SCdragF += SCelement.get_DragF(MET, vrel, rho, mach1) * SCelement.SCED.amount
         
        return SCdragF

    def get_AreaCd(self, MET, mach=0.0):
        # returns sum of area * Cd values (needed for external drag computation)
    
        # initialize cummulated values
        AreaCd = 0
        
        #for SCelement in self.SCelementsList:
        for i, SCelement in enumerate(self.SCelementsList):
             # drag
             AreaCd += SCelement.get_Cd(MET, mach) * SCelement.get_DragArea(MET) * SCelement.SCED.amount
         
        return AreaCd

    def get_staticvalues(self):
        
        # returns static values
        return self.SCD.staticValues.mass, \
               self.SCD.staticValues.dragarea, \
               self.SCD.staticValues.Cd, \
               self.SCD.staticValues.Cr, \
               self.SCD.staticValues.SRParea
        
        
    
    def set_fixedThrottlePointer(self, MET):
        
        # set fixed throttle pointer in all elements
        for SCelement in self.SCelementsList:
            SCelement.set_fixedThrottlePointer(MET)
            
        return 1
     
    def reset_fixedThrottlePointer(self):
        
        # set fixed throttle pointer in all elements
        for SCelement in self.SCelementsList:
            SCelement.fixedThrottlePointer = -1
        
        return 1
    
    def get_EventsList(self):
        
        # empty list
        self.eventsDF = pd.DataFrame(index = range(0), columns=['MET','eventType'])
        
        # returns event list
        for SCelement in self.SCelementsList:
            self.eventsDF = pd.concat([self.eventsDF, SCelement.eventsDF], ignore_index=True)
            
        # sorts by MET
        self.eventsDF = self.eventsDF.sort_values(by=['MET']).reset_index(drop=True)
        
        return self.eventsDF
        
        




class SpaceCraftElement():
    
    def __init__(self, sced, amount=1):
        
        # save spacraft element data
        self.SCED = sced  

        # save amount of this kind of element (e.g. for 2 boosters, 20 payloads, ...)         
        self.SCED.amount = amount
        
        # initialize fixedThrottlePointer
        self.fixedThrottlePointer = -1
        
        # initialize thrust/mass/area table
        self.init_throttleProfile()
        
        # initialize Cd interpolation
        self.init_Cdinterp()
        
        # make empty event dataframe
        self.eventsDF = pd.DataFrame(index = range(0), columns=['MET','eventType'])
        # set dtype, because otherwise it's object
        self.eventsDF[['MET','eventType']] = self.eventsDF[['MET','eventType']].astype('float')
        
        
    def set_fixedThrottlePointer(self, MET):
        
        self.fixedThrottlePointer = ((MET - self.SCED.throttleProfile.MET.values)>=0).nonzero()[0][-1]
       
        
    def init_throttleProfile(self):
        
        print('MRS:\t\tInitializing spacecraft throttle profile: ', self.SCED.name)

        # add empty table if no profile table was provided
        if not hasattr(self.SCED, 'throttleProfileInit'): 
            self.SCED.throttleProfileInit = pd.DataFrame([
                            [-20,       0,     -1,         0,          0,                0,   ''],     # Sim Start
                            ], 
                   columns= ['MET', 'partID', 'engineID', 'actEngines', 'throttle', 'gradMode', 'eventDesc']) 
        

        #
        # Adjust for non conform time steps
        #
                
        # add interpolation steps for sub-timestamp profiles (hardcoded to fixed step size)
        # valid step size used in the simulator
        validstep = 0.1
        
        # make a copy of throttleProfile to work with
        self.SCED.throttleProfile = self.SCED.throttleProfileInit.copy()
        
        # pointer to profile 
        i = 1 # starting in second profile, first one needs to be ok!
        # for loop not possible, because profiles are added during exeuction 
        while True: 
            
            # break when through 
            if i == len(self.SCED.throttleProfile):
                break
            
            # difference of MET to valid step size; going through Decimal numbers to be more accurate 
            diff2stepsize = float(Decimal(str(self.SCED.throttleProfile.MET.values[i])) % Decimal(str(validstep)))
            # check if MET is NOT rounded to 1/10 of a second, i.e. diff2tenth != 0
            if diff2stepsize:
                
                #check if stepwise change
                if self.SCED.throttleProfile.gradMode[i-1] == 0:
                    
                    print('MRS:\t\tWARNING: invalid MET value in throttle profile of ', self.SCED.name)
                    print('MRS:\t\tWARNING: - trying to recover static throttle profile.')
                    
                    # throttle before change
                    throttleStart = self.SCED.throttleProfile.throttle.values[i-1]
                    # throttle after change 
                    throttleEnd = self.SCED.throttleProfile.throttle.values[i]
                    
                    # ratio of difference rel. to valid step size
                    ratioDiff2Step = diff2stepsize/validstep
                    
                    # get present throttle profile settings
                    profileSettings = self.SCED.throttleProfile.loc[i].copy()
                    
                    # get MET with valid step size before the change
                    validMET = profileSettings.MET - diff2stepsize
                    
                    # new profiles
                    profileStart = profileSettings.copy()
                    profileIntermediate = profileSettings.copy()
                    profileEnd = profileSettings.copy()
                    
                    # if change occurs in first half of valid step size
                    if ratioDiff2Step <= .5:
                        # calculated intermediate value for new profile 
                        throttleIntermediate = throttleStart * (.5 + ratioDiff2Step) + \
                                               throttleEnd * (.5 - ratioDiff2Step)
                                               
                        # set up new profiles: 
                        # start profile
                        profileStart.MET = round(validMET - validstep, 3) # MET
                        profileStart.throttle = throttleStart # throttle start value
                        profileStart.gradMode = 1 # set to gradient
                        profileStart.eventDesc = profileStart.eventDesc + ' (adjusted)' 
                        
                        # intermediate profile
                        profileIntermediate.MET = round(validMET, 3) # MET
                        profileIntermediate.throttle  = throttleIntermediate # throttle intermediate value
                        profileIntermediate.gradMode = 1 # set to gradient
                        profileIntermediate.eventDesc = '' #profileIntermediate.eventDesc + ' (adjusted)' 
                        
                        # end profile
                        profileEnd.MET =  round(validMET + validstep, 3) # MET
                        profileEnd.throttle  = throttleEnd # throttle end value
                        profileEnd.gradMode *= 1 # keep unchanged
                        profileEnd.eventDesc = '' #profileEnd.eventDesc + ' (adjusted)' 
                                               
                    # change of power in second half of valid step size
                    else:
                        # calculated intermediate value for new profile 
                        throttleIntermediate = throttleStart * (-.5 + ratioDiff2Step) + \
                                               throttleEnd * (1.5 - ratioDiff2Step)
                    
                        # set up new profiles: 
                        # start profile
                        profileStart.MET = round(validMET, 3) # MET
                        profileStart.throttle = throttleStart # throttle start value
                        profileStart.gradMode = 1 # set to gradient
                        profileStart.eventDesc = profileStart.eventDesc + ' (adjusted)' 
                        
                        # intermediate profile
                        profileIntermediate.MET = round(validMET + validstep, 3)  # MET
                        profileIntermediate.throttle  = throttleIntermediate # throttle intermediate value
                        profileIntermediate.gradMode = 1 # set to gradient
                        profileIntermediate.eventDesc = '' #profileIntermediate.eventDesc + ' (adjusted)' 
                        
                        # end profile
                        profileEnd.MET = round(validMET + 2 * validstep, 3) # MET
                        profileEnd.throttle  = throttleEnd # throttle end value
                        profileEnd.gradMode *= 1 # keep unchanged
                        profileEnd.eventDesc = '' #profileEnd.eventDesc + ' (adjusted)' 
                    
                    # drop old throttle profile
                    self.SCED.throttleProfile.drop(i)
                    
                    # insert new profiles
                    self.SCED.throttleProfile.loc[i] = profileStart
                    self.SCED.throttleProfile.loc[i + .25] = profileIntermediate
                    self.SCED.throttleProfile.loc[i + .5] = profileEnd
                    self.SCED.throttleProfile = self.SCED.throttleProfile.sort_index().reset_index(drop=True)
                    
                    # jump by three values, because two where added
                    i += 3
            
                elif self.SCED.throttleProfile.gradMode[i-1] == 1:
                    print('MRS:\t\tWARNING: invalid MET value in throttle profile of ', self.SCED.name)
                    print('MRS:\t\tWARNING: - recovery not possibe, because gradient profiles not implemented yet.')
                    
                    # not implemented, therefore i += 1
                    i += 1
                    
                
            # good MET values, do nothing   
            else:
                
                # jump to next profile
                i += 1

   
        # add staging event (i.e. a throttle profile with updated masses, but w/o thrust)
        for i in range(len(self.SCED.parts)):
            # find throttle profile with new stages
            profilePointer = ((self.SCED.parts.stagingTime[i] - self.SCED.throttleProfile.MET.values)>=0).nonzero()[0][-1]
            # add a throttle profile for the new active stage (between last old stage and first new stage profile)
            self.SCED.throttleProfile.loc[profilePointer +.5,['MET','partID','engineID','actEngines','throttle','gradMode','eventDesc', 'kindOfProfile']] = \
                self.SCED.parts.stagingTime[i], i+1, -1, 0, 0, 0, '', 'staging'
            self.SCED.throttleProfile = self.SCED.throttleProfile.sort_index().reset_index(drop=True)

        # add throttle gradients
        for i in range(len(self.SCED.throttleProfile)-1):
            # no gradient, fixed throttle
            if self.SCED.throttleProfile.gradMode[i] == 0:
                self.SCED.throttleProfile.loc[i,'throttleGrad'] = 0
            # linear gradient
            elif self.SCED.throttleProfile.gradMode[i] == 1:
                throttleGrad = (self.SCED.throttleProfile.throttle[i+1]-self.SCED.throttleProfile.throttle[i])/ \
                               (self.SCED.throttleProfile.MET[i+1]-self.SCED.throttleProfile.MET[i])
                self.SCED.throttleProfile.loc[i,'throttleGrad'] = throttleGrad
        # last entry has a zero gradient
        self.SCED.throttleProfile.loc[i+1,'throttleGrad'] = 0
        
        # add thrust and fuel values, gradients and drag areas
        partID = -1 # init partID pointer
        for i in range(len(self.SCED.throttleProfile)-1):
            
            # check if new part is active
            if self.SCED.throttleProfile.partID[i]>partID:
                partID = int(self.SCED.throttleProfile.partID[i]) # update current pointer
                self.SCED.throttleProfile.loc[i,'fuelInit'] =  self.SCED.parts.fuelMass[partID] # add initial fuellMass
            # if part was already active
            else:
                # use previously calculated remaining fuel value
                self.SCED.throttleProfile.loc[i,'fuelInit'] = self.SCED.throttleProfile.fuelRemain[i-1]
                
            # if not in staging and engine is provided, add values for running engines 
            if self.SCED.throttleProfile.kindOfProfile[i] != 'staging':
                
                # add fuel + thrust only if engine is provided
                if self.SCED.throttleProfile.engineID[i]>=0:
                    
                    # write kind of element
                    self.SCED.throttleProfile.loc[i,'kindOfProfile'] = 'throttling'  
                    
                    # fuel gradient (static throttle)
                    fuelGrad = - self.SCED.throttleProfile.actEngines[i] * self.SCED.throttleProfile.throttle[i] \
                               * self.SCED.engines.fuelFlow[self.SCED.throttleProfile.engineID[i]]
                    self.SCED.throttleProfile.loc[i,'fuelGrad'] = fuelGrad        
                    
                    # fuel gradient (gradient throttle)
                    fuelGrad2 = - self.SCED.throttleProfile.actEngines[i] * self.SCED.throttleProfile.throttleGrad[i] \
                                * self.SCED.engines.fuelFlow[self.SCED.throttleProfile.engineID[i]] * .5
                    self.SCED.throttleProfile.loc[i,'fuelGrad2'] = fuelGrad2   
                    
                    # add thrust values (per engine)
                    self.SCED.throttleProfile.loc[i,'thrustSL'] =  self.SCED.engines.thrustSL[self.SCED.throttleProfile.engineID[i]]
                    self.SCED.throttleProfile.loc[i,'thrustVAC'] =  self.SCED.engines.thrustVAC[self.SCED.throttleProfile.engineID[i]]
        
                # otherwise, if no engine was provided
                else:
                    # write kind of element
                    self.SCED.throttleProfile.loc[i,'kindOfProfile'] = 'staticSC'  
                    
                    self.SCED.throttleProfile.loc[i,['fuelGrad','fuelGrad2', 'thrustSL','thrustVAC']] = 0
        
            # otherwise, if in staging
            elif self.SCED.throttleProfile.kindOfProfile[i] == 'staging':
                self.SCED.throttleProfile.loc[i,['fuelGrad','fuelGrad2', 'thrustSL','thrustVAC']] = 0
            
            # remaining fuel
            profileDur = self.SCED.throttleProfile.MET[i+1] - self.SCED.throttleProfile.MET[i]
            fuelRemain = self.SCED.throttleProfile.fuelInit[i] + \
                         self.SCED.throttleProfile.fuelGrad[i] * profileDur + \
                         self.SCED.throttleProfile.fuelGrad2[i] * profileDur**2
            
            # check if remaining fuel below zero
            if fuelRemain <= 0:
                print('MRS:\t\tWARNING: out of fuel in profile ', i, ' of ', self.SCED.name)
                fuelRemain = 0
            
            # store remaining fuel mass
            self.SCED.throttleProfile.loc[i,'fuelRemain'] = fuelRemain  
                
            # add dry masses + fuel of remaining stages     
            drymass    = self.SCED.parts.dryMass[partID:].values.sum()
            fuelRemain = self.SCED.parts.fuelMass[partID+1:].values.sum()
            self.SCED.throttleProfile.loc[i,'remainingMass'] = drymass + fuelRemain
            
            # add drag area
            self.SCED.throttleProfile.loc[i,'dragArea'] = self.SCED.parts.dragArea[partID] 
            
        # additional values for last row, when everything has staged     
        self.SCED.throttleProfile.loc[i+1,['fuelInit','fuelGrad','fuelGrad2','fuelRemain','remainingMass','dragArea', 'thrustSL','thrustVAC']] = 0
            
        # adding descriptions to the profiles
        for i in range(len(self.SCED.throttleProfile)):
            # when staging
            if self.SCED.throttleProfile.kindOfProfile[i] == 'staging':
                profileDesc = self.SCED.parts.name[self.SCED.throttleProfile.partID[i-1]] + ': Staging!' 
            elif self.SCED.throttleProfile.kindOfProfile[i] == 'staticSC':
                profileDesc = 'No use of engine(s).'   
            elif self.SCED.throttleProfile.kindOfProfile[i] == 'throttling':
                
                # when fixed throttle value
                if self.SCED.throttleProfile.throttleGrad[i] == 0:
                    profileDesc = str(int(self.SCED.throttleProfile.actEngines[i])) + ' engine(s) at ' + str(self.SCED.throttleProfile.throttle[i]*100) + '%.'
                # when ramping down
                elif self.SCED.throttleProfile.throttleGrad[i] < 0:
                    profileDesc = str(int(self.SCED.throttleProfile.actEngines[i])) + ' engine(s) at ' + str(self.SCED.throttleProfile.throttle[i]*100) + '%,' + \
                        ' throttling down to ' + str(self.SCED.throttleProfile.throttle[i+1]*100) +  '%.'
                # when ramping up
                elif self.SCED.throttleProfile.throttleGrad[i] > 0:
                    profileDesc = str(int(self.SCED.throttleProfile.actEngines[i])) + ' engine(s) at ' + str(self.SCED.throttleProfile.throttle[i]*100) + '%,' + \
                        ' throttling up to ' + str(self.SCED.throttleProfile.throttle[i+1]*100) +  '%.'
                else:
                    profileDesc = 'Error (unknown kind of throttle profile)'   
                    
                # add part name
                partname = self.SCED.parts.name[self.SCED.throttleProfile.partID[i]]
                profileDesc = partname + ': ' + profileDesc
                
            else:
                profileDesc = 'Error (unknown kind of throttle profile)'   
            
            # add MET
            profileDesc = str(self.SCED.throttleProfile.MET[i]) + '\t\t' + profileDesc
            
            # save to table, including bool variable if message was displayed
            self.SCED.throttleProfile.loc[i,['profileDesc','descPrinted']] = profileDesc, 0
            
        
        return 0 


    
    def get_ThrustMassThrottleFlowrateRemFuel(self, MET, pa=0.0, verbose=''):
        
        # get throttle pointer (tP) in throttle profile if not preset 
        if self.fixedThrottlePointer == -1:
            tP =  ((MET - self.SCED.throttleProfile.MET.values)>=0).nonzero()[0][-1]
        # otherwhise tP should be already set
        else:
            tP = self.fixedThrottlePointer 
            
        # show message for current profile
        if not self.SCED.throttleProfile.descPrinted[tP]:
            # show only if required
            if verbose=='v':
                print('MRS:\t\t', self.SCED.throttleProfile.profileDesc[tP])
            # save that description was printed
            self.SCED.throttleProfile.loc[tP,'descPrinted'] = 1
            
            # if event description provided
            if not self.SCED.throttleProfile.eventDesc[tP] == '':
                print('MRS: Event at MET=',MET,'seconds :',self.SCED.throttleProfile.eventDesc[tP])
                # add event to events dataframe of spacecraft element
                
                self.eventsDF.loc[len(self.eventsDF), ['MET','eventType']] = \
                                   MET, self.SCED.throttleProfile.eventDesc[tP]
                
                
                
        # calc time in current profile
        timeInProfile = MET - self.SCED.throttleProfile.MET[tP]
        
        # calc current throttle value
        throttle = self.SCED.throttleProfile.throttle[tP] +\
                   self.SCED.throttleProfile.throttleGrad[tP] * timeInProfile
        # calc thrust (for SL and VAC)
        thrustSL = throttle * self.SCED.throttleProfile.thrustSL[tP] * self.SCED.throttleProfile.actEngines[tP]
        thrustVAC = throttle * self.SCED.throttleProfile.thrustVAC[tP] * self.SCED.throttleProfile.actEngines[tP]
        
        # calc thrust of spacecraft element
        # 1. Calc gradient for Thrust w.r.t to atmospheric pressure)
        thrustVSpressure = (thrustSL-thrustVAC) / 101325.0
        # 2. Calc thrust for given pressure
        thrust = thrustVAC + thrustVSpressure * pa
        
        # calc mass of spacecraft element
        # 1. Calc mass of fuel of active part
        fuelPart = self.SCED.throttleProfile.fuelInit[tP] +\
                   self.SCED.throttleProfile.fuelGrad[tP] * timeInProfile + \
                   self.SCED.throttleProfile.fuelGrad2[tP] * timeInProfile**2
        # 2. Add mass of other masses (dry mass of current part + complete other parts)
        mass = fuelPart + self.SCED.throttleProfile.remainingMass[tP]
       
        # calc flowrate
        flowrate = - (self.SCED.throttleProfile.fuelGrad[tP] +\
                    2*self.SCED.throttleProfile.fuelGrad2[tP] * timeInProfile)
            
        # reset values if no fuel is left
        if fuelPart <=0:
            fuelPart = 0 
            thrust = 0
            flowrate = 0
       
        return thrust, mass, throttle, flowrate, fuelPart, self.SCED.amount
    
    def get_DragF(self, MET, vrel=np.zeros((3)), rho=0.0, mach1=0.0):
        
        # calc drag force for given velocity, density and mach number
        DragF = - 0.5 * vrel * np.linalg.norm(vrel) * self.get_Cd(MET, vrel/mach1)\
                * self.get_DragArea(MET) * rho
        
        return DragF
    
    def init_Cdinterp(self):
        
        # if no Cd table provided
        if not hasattr(self.SCED, 'C_D'): 
            self.CDinterpActive = 0
        # if Cd table provided
        else:
            self.CDinterpActive = 1
            self.CDinterp = interpolate.interp1d(self.SCED.C_D[:,0], self.SCED.C_D[:,1], kind='cubic')
            
        return 1
    
    def get_Cd(self, MET, mach):
        
        # return 0 if no Cd interpolation not available
        if self.CDinterpActive == 0:
            return 0
        
        # if faster than highest mach value, provide last C_D value
        if mach>self.SCED.C_D[-1,0]:
            return self.SCED.C_D[-1,1]
        
        return self.CDinterp(mach)
    
    def get_DragArea(self, MET):
        
        # get throttle pointer (tP) in throttle profile if not preset 
        if self.fixedThrottlePointer == -1:
            tP =  ((MET - self.SCED.throttleProfile.MET.values)>=0).nonzero()[0][-1]
        # otherwhise tP should be already set
        else:
            tP = self.fixedThrottlePointer 
            
        # return drag area
        return self.SCED.throttleProfile.dragArea[tP]
    
    def get_totalImpulse(self):
        # return the total impulse (Ns) for spacecraft element
        
        # get throttleprofile
        throttleProfile = self.get_throttleProfile()
        
        # calculate total impulse (by summing up integrals of trapezoids)
        impulse = ((throttleProfile.throttle[1:].values + throttleProfile.throttle[:-1].values)/2 * \
                   (throttleProfile.MET[1:].values - throttleProfile.MET[:-1].values) * \
                   throttleProfile.fullthrust[:-1].values).sum()

        return impulse
        
    def get_throttleProfile(self):
        # Adds end values for fixed power profiles and returns the throttle profile.
        # Used to calculate the total impulse or to plot the throttle profile.
        
        # make a copy to process profile
        throttleProfile = self.SCED.throttleProfile[['MET','throttle','actEngines', 'engineID', 'gradMode']].copy()
        
        # pointer to profile ID
        i = 0
        
        # normale looping through elemeents not possible, because adding elements while looping through table
        while True:
           # if in last profile, break
            if i==len(throttleProfile)-1:
                break
           
            # add full thrust value (actEngines * thrust for given engine)
            # get enginesID
            enginesID = throttleProfile.engineID[i]
            if enginesID >-1:
                fullthrust = throttleProfile.actEngines[i] * self.SCED.engines.thrustVAC[enginesID]
            else:
                fullthrust = 0
            
            # add fullthrust value to dataframe
            throttleProfile.loc[i,'fullthrust'] = fullthrust
            
            # detect non-gradient profiles
            if throttleProfile.gradMode[i] == 0:
                
                # make new profile
                profileSettings = throttleProfile.loc[i].copy()
                
                # new settings
                profileSettings.MET = throttleProfile.MET[i+1]
                
                # store new profile
                throttleProfile.loc[i + .5] = profileSettings
                throttleProfile = throttleProfile.sort_index().reset_index(drop=True) 
                 
                # increase by to to jump over the newly created profile
                i += 2
            else:
                i += 1
                
        # add fullthrust value to last profile (is zero per definition)
        throttleProfile.fullthrust[i] = 0
                
        # return cleaned throttle profile
        return throttleProfile
    
    
    
    
    
    