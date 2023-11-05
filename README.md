# Welcome to MyRocketSimulator!

MyRocketSimulator is a Python library for spacecraft orbit propagation around the Earth and missions towards and around the Moon. It provides the simulation environment for
the author’s plan to replicate the complete trajectory of Apollo 11, from launch to splashdown. Necessary parts for such an undertaking a continuously added to MRS. 

The simulator relies on a high-fidelity propagator that includes all relevant perturbating forces, such as drag, SRP, gravity with spherical harmonics and third body gravity. A published in-depth comparison with GMAT proved its accuracy for Earth orbiting spacecrafts. Download the paper [Spacecraft Orbit Propagation in an Open-Source Python Environment]( https://www.researchgate.net/publication/375293398_Spacecraft_Orbit_Propagation_in_an_Open-Source_Python_Environment).

An early, non-published release of MRS was used to perform a preliminary flight simulation and analysis of the Moon-bound Artemis I mission. Download the paper [Preliminary Launch Trajectory Simulation for Artemis I with the Space Launch System]( https://www.researchgate.net/publication/362270344_Preliminary_Launch_Trajectory_Simulation_for_Artemis_I_with_the_Space_Launch_System).

Relevant features of MRS 1.0:
-	Multi-segment simulations providing distinct force configurations and delta-v maneuvers. 
-	Earth and Moon gravity with spherical harmonics using [pyshtools]( https://shtools.github.io/SHTOOLS/).
-	Third bodies (Sun, Venus, Mars, Jupiter, Moon).
-	Precise ephemeris and coordinate transformations through [Skyfield]( https://rhodesmill.org/skyfield/).
-	Atmospheric modelling with [PyNRLMSISE-00](https://github.com/st-bender/pynrlmsise00).
-	DOP853 integration.
-	Data frame export, including more than 60 run-time variables.
-	Import of external state vectors for 1:1 trajectory comparison.
-	Pre-configured visualizations. 

Examplary visualization of flight with two delta-v maneuvers in order to raise the spacecraft's altitude:
![GCRF view of satellite with Hohmann transfer to higher altitude](https://raw.githubusercontent.com/ThibaultBS/MyRocketSimulator/main/MRS_examples/MRSoutput/MRSexample3_GCRForbit.svg)

## Installation
### PyPI
MRS can be automatically installed trough its release on PyPI ([Link](https://pypi.org/project/myrocketsimulator/)). Simply execute the following command: 

`pip install myrocketsimulator`

### Github
The MRS library can be manually downloaded from Gibthub ([Link](https://github.com/ThibaultBS/MyRocketSimulator)). After downloading the repository, simply execute the following command in its top folder:

`pip install .`

## First simulation
After installation, you are ready to start with your first simulation. The following code performs a 24-hour propagation of the ISS, using the included MRS default mission. Please note that Skyfield and spaceweather will require an internet connection to download relevant data, such as ephemeris files.

```python
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
MRSdemo0viewer.plot_GroundtrackEarth()
```
The output will be:
```
MRS:		Using default mission.
MRS:		Loading mission object 'defaultMRSmission'.
MRS:		Mission 'Default MRS Mission' loaded.
MRS:		Checking mission data validity.
MRS:		Loading Default MRS spacecraft as static spacecraft.
MRS:		Mission data valid.
MRS:		Running mission Default MRS Mission Update.
MRS:		Processing mission segment 0.
MRS:		Mission ended. Processing time: 80.276 seconds.
MRS:		Adding Earth LLA to dataframe.
MRSviewer:	Loading dataframe of mission Default MRS Mission Update.
```






