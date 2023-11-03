import datetime
import numpy as np
import pandas as pd

class MRSmissionData():
    
    # Mission settings
    name = 'CREME Default Mission'
    launchtype = 0 # 0 from starts from state vector, 1 starts from launchsite at first mission segment
    
    # t0_UTC is needed for launchtype 0 (start from state vector) and launchtype 1 (start from launch pad)
    # t0_JD has priority if t0_JD>0, otherwise it gets calculated for given t0_UTC
    t0_JD = 29945.50080073974 + 2430000.0  # TDB; here expresse in ModJulian (GMAT) + Offset
    #t0_UTC = datetime.datetime(2000,1,1,12,0,0,0) # [y,m,d,h,m,s,mus]
    t0_UTC = 0
   
    # For launchtype 1, is ignored and set to 0; for launchtype 0, t0 might be offseted to MET=0.
    t0_MET = 0. # indicates the MET for given t0 time and y0 state vector
    tend_MET = 0 # MET preliminary end of simulation (!=0 to be used)
 
    # launch settings; mission starts at very first mission segment 
    launchsite_name = ''
    launchsite_LLA = np.array([0.,0.,0])  # [deg, geodetic], [deg], [m]
    

    
    y0 = np.array([
        4160.654753643215, -5602.074997808763, 12.06646228442704,   # [m]
        -0.8311739944172626, -0.6011835366131142, 7.48792556529475  # [m/s]
        ]) * 1000

      
    
    class spacecraft():
        
        name = 'CREME'
    
        class staticValues():
            mass = 4.00 # [kg]
            dragarea = 0.033 # [m^2]
            Cd = 2.2 
            Cr = 1.5
            SRParea = 0.033  # [m^2] 
      
    # starts in in first line if launching, otherwise starts in row of t0_MET and uses state vector
    missionSegments = pd.DataFrame([
                    [0.,      0,      0,        'Start of propagation'],
                    [5000,   1,      0,        'deltaV maneuver #0'],
                    [5000+2938.69,   1,      1, 'deltaV maneuver #1'],
                    [86400,   0,      0,        'End of simulation'] 
                    ], 
           columns= ['MET','type','configID', 'comment'])
    
    propaSettings = pd.DataFrame([
                    [1,     'DOP853',   1,            0,        10,     'auto-step size, no thrust, DOP853'],
                   ], 
           columns= ['mode','method','stepsizePropa','forcesID','downsampleLog','comment'])
  
    forcesSettings = pd.DataFrame([
                    [0,       0,   ['Earth'], 'nrlmsise00',                0, 0, 0,  'Forces Settings']
                   ], 
           columns= ['EarthSHn','MoonSHn','planets','atmosModel','drag','SRP', 'activeSC', 'comment'])
    
    maneuverSettings = pd.DataFrame([
                    ['VNB', 32.640, 0.0, 0.0, 'args', 'Earth', '-'],
                    ['VNB', 32.500, 0.0, 0.0, 'args', 'Earth', '-']
                    ],
            columns = ['frame', 'dx', 'dy', 'dz', 'args', 'planet', 'comment'])
    
    
    missionEvents = pd.DataFrame([
                    ],
            columns = ['MET', 'eventType'])
    
    