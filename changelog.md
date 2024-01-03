# Changelog

Relevant changes to the MyRocketSimulator are summarized in this changelog. 

## [1.1.1] - 2023-12-31

### Added

- MRSlib get_liftVector()
- MRSlib get_AoAangle()
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
