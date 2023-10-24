# 1.1

1.1.0, 2023-10-23

* New functions and updates
    - Updated the angulardistance() function to use an updated form of the Vincenty equation, which is a special case for perfect spheres and now valid for arbitrarily small angles.  (The previous form suffered from occasional numerical instability and small angle errors)  Now users can specify returncomponents='cartesian' to get longitude component calculated along line of constant latitude (i.e. with cos(latitude) taken into account) or returncomponents='spherical' to get spherically orthogonal components (following orthogonal great circle arcs).
    - Added function calc_optimal_slew_loop() for determining an optimal order for slewing between a handful of targets. If slew rates (degrees/minute) are known for AZ/EL axes, the time per slew loop is calculated for each permutation.
    - Added simple cone search functions to find sources within a certain radius of a reference position/source.  sources_within_radius() uses an input list of ephem targets; srctable_within_radius() tests source separations from a table (e.g. pandas/astropy) of source coordinates, which is simpler for working with more traditional (and large!) source lists.
    - Added functionality to nearest_from_target_list() and nearest_standard_fringefinder() -- now they can accept a list of obstargets (for the nearest average/median source pos from a group of input targets).  Options to return the ephem target or the separations instead of simply the target name have been included.
    - Changed name of plot_night_observing_tracks() to the slightly more general plot_observing_tracks() -- without the 'night'. And added new convenience function for plotting the nighttime version, called  plot_night_observing_tracks() , which simply sets light_fill=False
    - Also added convenience function for plot_day_observing_tracks, which simply calls plot_observing_tracks with light=True
    - Upgraded obs.skysep_fixed_single() to have an option to return the longitude/latitude components in addition to the total separation, like with obs.angulardistance.  Also added simple function wrap_pmPI() to wrap radian values to range +/- PI.  
    - Added convenience functions to make a plot of sun(&moon) separations from a target source over the course of a specified year
    - Added convenience functions for printing daily sunseps and moonseps, with option to only print values for a candence of every N days
    - Added convenience function print_VLBA_observability_summary().  With minimal configuring, prints source observability every N days in the specified range, and also estimates start times dynamic scheduling.
    - Added convenience function eph2c() to take in an ephem source and output the RA,DEC coordinates. Output styles are degrees, radians, sexagesimal using obs.dec2sex and optional kwargs, or astropy.coordinates.SkyCoord


# 1.0

1.0.2, 2023-09-24

* Minor fixes and tweaks (1.0.2)
    - Merged pull request to remove numpy.warnings filter, due to it being removed from recent np versions and crashing obsplanning on import
    - package multicolorfits has been moved to 'extras' instead of a requirement.  Functions for plotting images (finder plots) will fail without it, but it requires QT - maintaining QT compatibility on various systems and software versions has become a headache, and the rest of the obsplanning package need not be held hostage by the image plotting dependencies.


1.0.0, 2021-10-03

* Initial release
    - pyephem-based time conversions, observers (telescope), and target source (sky) objects
    - Plotting/visualization tools for target elevation tracks and finder plots
    - Transit time and target visibility calculations for single telescopes or multiple stations (VLBI/VLBA)
    - Miscellaneous coordinate handling, radio astronomy info, and other helpful functions
    - Preliminary versions use functional approach
