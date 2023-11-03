import datetime
import numpy as np
import pandas as pd

class MRSmissionData():
    
    # Mission settings
    name = 'Default MRS Mission'
    launchtype = 0 # 0 from starts from state vector, 1 starts from launchsite at first mission segment
    
    # t0_UTC is needed for launchtype 0 (start from state vector) and launchtype 1 (start from launch pad)
    # t0_JD has priority if t0_JD>0, otherwise it gets calculated for given t0_UTC
    t0_JD = 0 # TDB 
    #t0_UTC = datetime.datetime(2000,1,1,12,0,0,0) # [y,m,d,h,m,s,mus]
    t0_UTC = datetime.datetime.fromisoformat('2023-01-09T12:00:00.000')
   
    # For launchtype 1, is ignored and set to 0; for launchtype 0, t0 might be offseted to MET=0.
    t0_MET = 0. # indicates the MET for given t0 time and y0 state vector
    tend_MET = 0 # MET preliminary end of simulation (!=0 to be used)
 
    # launch settings; mission starts at very first mission segment 
    launchsite_name = 'no launchsite available'
    launchsite_LLA = np.array([0.,0.,0])  # [deg, geodetic], [deg], [m]
    

    
    y0 = np.array([
        #0., 0., 0., 
        6091.160573649400, 2321.823227910050, -1926.816516037270 ,   # [m]
        -0.21429216175361, 5.21216680526399, 5.60755423362593    # [m/s]
        ]) * 1000

    
    class spacecraft():
        
        name = 'Default MRS spacecraft'
        
        """
        # SCelements in here for debug 
        SCelements = pd.DataFrame([
                        ['stages', 1],
                        ],
               columns= ['name','amount'])
        """
        
        class staticValues():
            mass = 472355.00 # [kg]
            dragarea = 1622.90 * 0.89 # [m^2]
            Cd = 2.40 
            Cr = 1.8
            SRParea = 2500  # [m^2] 
      
      
    # starts in in first line if launching, otherwise starts in row of t0_MET and uses state vector
    missionSegments = pd.DataFrame([
                    [0.,      0,      1,        'Start of propagation'],
                    #[20000,   1,      0,        'deltaV maneuver'],
                    #[22790,   1,      1,        'deltaV maneuver'],
                    #[39200,   0,      2,        'Hi-Res'],
                    [86400,   0,      3,        'End of simulation'], # 86400 for whole day
                    ], 
           columns= ['MET','type','configID', 'comment'])
    
    propaSettings = pd.DataFrame([
                    [0,     '-',      0.1,            0,        10,     'Keeping LLA postion at launchsite'],
                    [1,     'DOP853',   1,            1,        30,     'LEO, auto-step size, no thrust, DOP853'],
                    [1,     'DOP853',   60,            0,        30,     'GMAT'],
                    [2,     'DOP853',   60,            1,       1,      '100 ms steps, DOP853']
                    ], 
           columns= ['mode','method','stepsizePropa','forcesID','downsampleLog','comment'])
  
    forcesSettings = pd.DataFrame([
                    [0,       0,   ['Earth'], 'nrlmsise00',                1, 0, 0,  'GMAT tests'],
                    [50,       0,   ['Earth', 'Sun', 'Moon'], 'nrlmsise00', 1, 0, 0, 'Earth orbit forces'],
                    ], 
           columns= ['EarthSHn','MoonSHn','planets','atmosModel','drag','SRP', 'activeSC', 'comment'])
    
    maneuverSettings = pd.DataFrame([
                    ['VNB',    2000,    0, 0,     [0,0], 'Earth',   'Test burn'],
                    ['LVLH',    -160,    0, 0,     [0,0], 'Earth',   'Test burn']
                    ],
            columns = ['frame', 'dx', 'dy', 'dz', 'args', 'planet', 'comment'])
    
    
    missionEvents = pd.DataFrame([
                    [10,   '10 seconds event'],
                    [999,  '999 seconds event']
                    ],
            columns = ['MET', 'eventType'])
    
    