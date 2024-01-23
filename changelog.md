# Changelog

Relevant changes to the MyRocketSimulator are summarized in this changelog. 


## [1.2.0] - 2024-02-01

### Added

- Numpy based MRSrocket.py
- Mission data can be provided as a class to __init__() or load_mission()
- get_FPAHAVEL_EF() can be called for round or elliptical Earth
- Several new expand_DF() types
- MRSvislib offers now a universal view for missionDF values: missionDFplotTS()
- MRSvislib get_eventProperty() returns a specific value for a given event

### Fixed

- Calculation of t0_JD_liftoff for pure launchtype=1 missions.
- Bugs in get_guidance() in MRSguidance
- MRSspacraft's SpacCraft() get_DragF() now returns single value when vrel is a float

### Changed

- Definitin of TempDF takes now place in propagate_Mode0()/propagate_Mode1()

### Removed

- Pandas based MRSrocket.py
- Mode 2 propagation (step wise integration loop)



## [1.1.1] - 2023-12-31

### Added

- MRSlib get_liftVector() (will be used later for lift force vector)
- MRSlib get_AoAangle() (will be used later for drag/lift calculation)
- MRSlib get_EarthRangeToLaunchsite()
- MRSguidance library in its first release.
- Guidance object for static guidance in MRSlib.
- Guidance vector can be added to mission dataframe.
- Guidance vector angles can be added to mission dataframe.
- Range (ground distance to launch site) added to mission dataframe.
- Angle of Attack added to mission dataframe.

### Fixed

- Nothing

### Changed

- DefaultMRSmission contains guidance data with explanations.
- MRSlib get_ENUvec_Earth() has frame as new option (WGS84 or TOD/TETE)
- get_acceleration() gets thrust direction from guidance object.
- Readme.md updates.

### Removed

- Nothing


## [1.1.0] - 2023-11-20

### Added

- First public release of MRS.
- Solid tides for Earth and Moon, to be activated in forcesSettings.
- Changelog added (changelog.md).

### Fixed

- Import of comparison state vectors with own timestamps.

### Changed

- Earth position is only preloaded if really needed (speed improvement).
- Readme.md update.

### Removed

- Hard-coded experimental tides function for Moon.

## [1.0.4] - 2023-11-11 [YANKED]

## [1.0.2] - 2023-11-03 [YANKED]
