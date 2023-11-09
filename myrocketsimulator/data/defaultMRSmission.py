import datetime
import numpy as np
import pandas as pd

class MRSmissionData():
    """
    The MRSmissionData class contains all the mission specific settings. 
    Additionally, MRS default settings can also be overwritten.
    
    This file can be used as a template for your own mission design. It 
    therefore contains useful information how to enter the data.
    
    Future updates of MRS might require adjusted MRSmissioData structures.
    
    """
    
    # Mission name: this name will be used for display and for the filenames 
    # of mission and event data frames. Can be any kind of string.
    name = 'Default MRS Mission (ISS)'
    
    # The launchtype defines from where the trajectory starts:
    #   - launchtype 0: propagation starts from ICRF state vector (defined below)
    #   - launchtype 1: rocket sits on launchpad (position defined below)
    #
    # MRS 1.0/1.1 does only support launchtype 0 (due to missing active rocket)
    launchtype = 0 
    
    # t0 is the initial reference time. It can be provided as Julian Date TDB
    # (JD TBD) or as UTC time. 
    # 
    # The use of the reference time depends from the chosen launch type. In case of:
    #   - launchtype 0: t0 is the time of the initial state vector y0 (see below);
    #                   If t0_JD>0, it is used. Otherwise, t0_UTC is used.
    #   - launchtype 1: t0 is the time at which MET (mission ellapsed time) = 0;
    #                   t0_UTC is always used, t0_JD is ignored
    t0_JD = 0 # JD TDB 
    t0_UTC = datetime.datetime.fromisoformat('2023-01-09T12:00:00.000') # UTC
    #t0_UTC = datetime.datetime(2000,1,1,12,0,0,0) # [y,m,d,h,m,s,mus],UTC alternative
    
    # The state vector y0 is relevant for missions starting in an orbit or while
    # ascent (launchtype 0). It provides three position and three velocity values.
    # All values are defined in ICRF/J2000 and are [m]/[m/s].
    # The default mission provides a state vector of the ISS for UTC provided above.
    y0 = np.array([
        6091.160573649400, 2321.823227910050, -1926.816516037270 ,   
        -0.21429216175361, 5.21216680526399, 5.60755423362593        
        ]) * 1000 # multiply by 1000 to get meters / meters per second
    
    # t0_MET defines the MET-value at the moment of t0_JD/t0_UTC. Default is 0.
    # The use depends from the launchtype:
    #   - launchtype 0: value is used as provided
    #   - launchtype 1: value is set to 0; for launches MET=0 for given t0_UTC
    t0_MET = 0. # MET seconds
    
    # It might be helpful to stop the mission at at time earlier than the 
    # MET-value in the last row of the mission segments (see below), e.g. to
    # conduct just a short flight to test new parameters.
    # The propagation is terminated at the given tend_MET time (if tend_MET != 0)
    tend_MET = 0. # MET seconds
 
    # In case of launchtype 1, a launchsite needs to be provided.
    # The launchsite_name can be any kind of string and is used for display only.
    # The launchsite_LLA is the actual latitude, longitude, and altitude of the
    # launchsite.
    #
    # MRS 1.0/1.1 does not support launches, this data is therefore not used.
    launchsite_name = 'No launchsite in use'
    launchsite_LLA = np.array([0.,0.,0])  # [deg, geodetic], [deg], [m]
    
    # MRSlib has several default values that configure numerous parts of the 
    # program and the libraries it's using. These values may be altered to
    # adjust the behaviour of the simulator. 
    
    ### solve_ivp settings; default mission provides low values for fast processing
    integrator_atol = 1e-5 # solve_ivp atol value [m]; default: 1e-9
    integrator_rtol = 1e-5 # solve_ivp rtol value; default: 1e-9
    integrator_max_step = 60 # solve_ivp max_step value [s]; default: 10.
    
    ### gravity settings
    # EarthGravModel = 'EGM96'  # possible are: 'EGM96', 'EGM2008'
    # EarthGravZonalOnly = 0    # use only zonal spherical harmonics
    # DEephemeris = 'DE421'     # possible are: 'DE421', 'DE440'
    # MoonTides = 0             # experimental; not fully implemtend yet
    # OE_TEME = 0               # if 0: OE calculated in ICRF, if 1: OE calc. in TEME
    # fastEphemeris = 0         # if set to 1, planet positions will be fixed (fast)

    ### atmospheric settings
    # use_spaceweather = 1      # use actual spaceweather (through 'spaceweather')
    # f107s = 150               # fidxed f107s to be used when not using spaceweather
    # f107  = 150               # fidxed f107 to be used when not using spaceweather
    # Ap = 3                    # fidxed Ap/Kp to be used when not using spaceweather
    # transToMSISE90 = 0        # experimental; modify density in drag function
    # RhoGain = 1.07            # density gain
    # RhoOffset = 3.0495e-14    # density offset
    

    # The missionSegments table provides a list of mission segments with unique 
    # settings. At least, two segments are needed (start, stop). Four values are
    # provided per segment:
    #   - MET: Segment start time (including), expressed in MET; in case of the 
    #          last segment, processing stops at its given MET 
    #   - type: what kind of propagation is required in this segment:
    #       - type 0: propagation in space (launchtype 0) or on ground (launchtype 1)
    #       - type 1: delta V maneuver, executed at the given time; propagation
    #                 then continues with settings of proevious segment
    #   - configID: row number for settings of this segment, either in:
    #       - for type 0 segments: propaSettings table
    #       - for type 1 segments: maneuverSettings table
    #   - comment: a string, only used for display purposes
    # In case a t0_MET value is provided, the simulation starts at t0_MET with
    # the segment in which includes t0_MET. 
    # In case tend_MET value is provided, the simulation ends with the execution
    # of the current segment when tend_MET is reached.
    #
    # The default mission contains a few (commented) segments to demonstrate
    # the use of the missionSegments table.    
    #   - 1st segment in which the simulation starts
    #   - 2nd/3rd segment (commented): delta-v for orbit raise by 100 km 
    #   - 4th segment (commented): different propagation settings
    #   - 5th segment: only relevant for end time of previous segment
    missionSegments = pd.DataFrame([
                    [0.,       0,      1,        'Start of propagation'],
                    #[20000,   1,      0,        '1st deltaV maneuver'],
                    #[22819.5, 1,      1,        '2nd deltaV maneuver'],
                    #[70000,   0,      2,        'Hi-Res'],
                    [86400,    0,      1,        'End of simulation'], 
                    ], 
           columns= ['MET','type','configID', 'comment'])
    
    # The propaSettings table is used to store differen propagation settings.
    # This might be required for different phases of the flight.
    # Parameters are:
    #   - mode: how the propagation is performed:
    #           - mode 0: spacecraft is not propgated but standing on launchpad
    #           - mode 1: spacecraft is propagated by integration of acceleration;
    #                     This mode is prefered for static spacecrafts.
    #           - mode 2: spacecraft is propagated in discreete steps;
    #                     This mode is prefered for active spacecrafts, in order
    #                     to change guidance/thrust settings 
    #   - method: kind of integrator; possible are:
    #           - '-': no integrator; only possible for mode 0 
    #           - 'DOP853': solve_ivp DOP853; no other possible at the moment
    #   - stepsizePropa: integrator step size in mode 2; not relevant for 
    #                    the integrator in mode 0 or 1
    #   - forcesID: row number of forcesSettings table to be used 
    #   - downsampling: stepsizePropa * downsampling = step size for logging (all modes)
    #   - comment: string, only used for display
    
    propaSettings = pd.DataFrame([
                    [0,     '-',      0.1,            0,        10,     'Keeping LLA postion at launchsite'],
                    [1,     'DOP853',   1,            0,        60,     'LEO, auto-step size, no thrust, DOP853'],
                    [1,     'DOP853',   60,            0,        30,     'GMAT'],
                    [2,     'DOP853',   60,            1,       1,      '100 ms steps, DOP853']
                    ], 
           columns= ['mode','method','stepsizePropa','forcesID','downsampleLog','comment'])
  
    # The forcesSettings table contains different settings for forces acting
    # on the spacecraft. Depending on the flight phase, different settings might 
    # be helpful (only gravity and drag for launch, but all in orbit). 
    # Settings:
    #   - EarthSHn: degree/order for Earth spherical harmonics. If 0, point mass is used
    #   - MoonSHn: degree/order for Moon spherical harmonics. If 0, point mass is used
    #   - Planets: Planets to be used (for gravity and sun occultation); possible are (MRS 1.0/1.1):
    #       - Sun
    #       - Venus
    #       - Earth
    #       - Moon (even if it's not a planet :-)
    #       - Mars
    #       - Jupiter
    #   - atmosModel: only nrlmsise00 implemented in MRS 1.0/1.1. Can be '-' if no drag used.
    #   - drag: 0/1 wether to not use or use atmospheric drag force
    #   - SRP: 0/1 wether to not use or use solar radiation pressure force
    #   - activeSC: 0/1 wether to not use or to use active spacecraft; only 0 allowed for MRS 1.0/1.1
    #   - comment: string, only used for display 
    
    forcesSettings = pd.DataFrame([
                    [0,       0,   ['Earth'], '-',                0, 0, 0,  'No perturbating forces'],
                    [35,      0,   ['Earth', 'Sun', 'Moon'], 'nrlmsise00', 1, 1, 0, 'High fidelity simulation.'],
                    ], 
           columns= ['EarthSHn','MoonSHn','planets','atmosModel','drag','SRP', 'activeSC', 'comment'])
    
    # delta-v maneuvers are described in the maneuverSettings-table. Parameters:
    #   - frame: frame that described the xyz-axis of the maneuver. Implemented are:
    #           - VNB
    #           - LVLH
    #           - GCRF
    #   - dx/dy/dz: velocity deltas in x/y/z [m/s]
    #   - args: additional arguments within an array (not used in MRS 1.0/1.1)
    #   - planet: frame is defined w.r.t. to provided planet (Earth or Moon)
    #   - comment: string, only for display 
    
    maneuverSettings = pd.DataFrame([
                    ['VNB',    27.91,    0, 0,     [0,0], 'Earth',   '1st delta-v for Hohmann transfer'],
                    ['VNB',    27.81,    0, 0,     [0,0], 'Earth',   '2nd delta-v for circularization']
                    ],
            columns = ['frame', 'dx', 'dy', 'dz', 'args', 'planet', 'comment'])
    
    # It might be required to store trajectory data for givent MET values in the
    # event data frame. The missionEvent table defines these events.
    #   - MET: MET of event
    #   - eventType: string, name of event, only used for display/table entry 
    missionEvents = pd.DataFrame([
                    [10,   '10 seconds event'],
                    [999,  '999 seconds event']
                    ],
            columns = ['MET', 'eventType'])
    
    
    class spacecraft():
        """
        The spacecraft class contains data structures that describe the properties
        of the spacecraft. This information is mostly use to calculate acceleration
        values in other functions. 
        
        A spacecraft can be static (activeSC=0) or active (activeSC=1). The static
        values are always provided. Active spacecraft are not available in MRS 1.0/1.1.
        
        """
        
        # name of the spacecraft; string; only used for display
        name = 'Default MRS spacecraft (ISS)'
        
        # the staticValues-class provides fixed values for a static spacraft simulation.
        class staticValues():
            # the mass of the spacecraft in [kg]. Highly relevant for acceleration
            mass = 472355.00 
            # the area A used in the drag formula [m^2]
            dragarea = 1444
            # drag coefficient for the drag formula (no unit)
            Cd = 2.40 
            # reflectivity coefficient for the SRP formula (no unit)
            Cr = 1.8
            # the area A for the SRP formular in [m^2]
            SRParea = 2500  
            
            
    
            
            
            
            