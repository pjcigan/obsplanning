Determine observable targets from a source list
===============================================

If you have a list of target sources, and want to determine which of
them are observable within a particular observing time window,
obsplanning can help.

First, generate a source list (text file, or numpy array) in which each
row consists of the source name and coordinates. The coordinates can be
in either segagesimal or decimal. For example, here is a short list of
radio sources across the full RA range, which is saved in a file we will
call here 'some\_sources.txt':

::

    # [name, RA, DEC]
    J0001-1551, 00:01:05.32, -15:51:06.7
    J0137+0923, 01:37:51.10, +09:23:19.3
    J0319+5142, 03:19:56.06, +51:42:29.0
    J0518+2901, 05:18:06.50, +29:01:58.0
    J0733+5022, 07:33:52.53, +50:22:08.8
    J0909+0835, 09:09:12.15, +08:35:41.2
    J1100+0444, 11:00:11.46, +04:44:01.9
    J1230+1223, 12:30:49.42, +12:23:28.3
    J1404-0130, 14:04:45.88, -01:30:21.7
    J1550-0538, 15:50:29.85, -05:38:10.7
    J1728+0427, 17:28:24.96, +04:27:05.2
    J1910+2305, 19:10:45.10, +23:05:58.8
    J2108-2452, 21:08:12.32, -24:52:32.7
    J2257-3627, 22:57:10.61, -36:27:43.9

You could load this list into a numpy array with:

.. code:: python

    import numpy as np

    srclist=np.genfromtxt('some_sources.txt',delimiter=',',dtype=str)

If you have a list of source names but not their coordinates, you can
try automatically querying positions with
query\_object\_coords\_simbad(), like so:

.. code:: python

    import obsplanning as obs

    #obs.query_object_coords_simbad(stringname, return_fmt='dec')

    obs.query_object_coords_simbad('M1')
    # -->  [83.63308333333333, 22.0145]


    obs.query_object_coords_simbad('NGC1275',return_fmt='sex')
    # -->  ['03 19 48.1597', '+41 30 42.114']

Now, say you are scheduling observations with the Pietown antenna, on
2021 October 31, from 5:00 to 6:30 UTC, and you wish to observe one of
the sources from your list during this session. To determine the
observable targets at a single antenna, use
get\_visible\_targets\_from\_source\_list().

This function takes an ephem.Observer (telescope), beginning and end
times for the observations, and the source list. The input source list
can be given as either a path to the text file, or just input an Nx3
numpy array directly, e.g. how srclist was loaded above. If a text file
path is given, also specify the coordinate format, so that
numpy.genfromtxt can parse it correctly.

You can also specify the minimum telescope elevation angle at which to
consider a target 'visible', as well as the minimum number of minutes a
target must be up to be considered observable for that session. (If
specified as 'full', the default, then a target must be visible for the
full duration of the observation window.)

Here is an example, using the short target list from above, specifying
that a target must be up for at least 30 minutes:

.. code:: python

    Pietown = obs.create_ephem_observer('Pietown', '251:52:50.94', '34:18:03.61', 2365) 
    #obs.vlbaPT is also built in, but you can create an ephem observer for any station like above

    # The start/end times can be input as dt.datetime, ephem.Date, or 
    # strings formatted for ephem, like below:
    beginobs = '2021/10/31 05:00:00'
    endobs = '2021/10/31 06:30:00'

    obs.get_visible_targets_from_source_list(Pietown, beginobs, endobs, './some_sources.txt', 
        elevation_limit=15., minimum_observability_minutes=30., coord_format='sexagesimal')

Six of the 14 targets satisfy this condition observable for at least 30
minutes from Pietown in this window, due to being near the optimal RA
range, or at reasonably high declination.

.. code:: python

    #output = 
    [array(['J0001-1551', '00:01:05.32', '-15:51:06.7'], dtype='<U11'),
     array(['J0137+0923', '01:37:51.10', '+09:23:19.3'], dtype='<U11'),
     array(['J0319+5142', '03:19:56.06', '+51:42:29.0'], dtype='<U11'),
     array(['J0518+2901', '05:18:06.50', '+29:01:58.0'], dtype='<U11'),
     array(['J0733+5022', '07:33:52.53', '+50:22:08.8'], dtype='<U11'),
     array(['J1910+2305', '19:10:45.10', '+23:05:58.8'], dtype='<U11')]

But what if you want to select targets that are observable concurrently
among several different telescopes, for simultaneous observations? This
is a requirement for interferometry, and in particular VLBI. For this,
we can use get\_visible\_targets\_from\_\\source\_list\_multistation()

Consider the same observing setup as above, except that now we want to
observe with three antennas simultaneously: the Mauna Kea, Brewster, and
Saint Croix VLBA antennas. Now we can also specify the minimum number of
stations a target must be simultaneously visible from, in the uptime
duration calculation. Here we will use the default 'full', meaning it
must be visible from all three stations. But when using the full VLBA
array, for example, you might determine that 8 out of the 10 stations
are sufficient, etc.

.. code:: python

    station_array=[obs.vlbaMK, obs.vlbaBR, obs.vlbaSC]

    obs.get_visible_targets_from_source_list_multistation(station_array, beginobs, endobs,  
        './some_sources.txt', elevation_limit=15., decbin_limits_deg=[-90,90], 
        minimum_observability_minutes=30, minimum_mutual_vis_observers='full', coord_format='sex')

Now the number of targets that are observable with these constraints has
dropped to three - the requirement to be visible at both Mauna Kea and
St Croix is much more restrictive than for a single station on its own.

.. code:: python

    #output = 
    [array(['J0001-1551', '00:01:05.32', '-15:51:06.7'], dtype='<U11'),
     array(['J0137+0923', '01:37:51.10', '+09:23:19.3'], dtype='<U11'),
     array(['J0319+5142', '03:19:56.06', '+51:42:29.0'], dtype='<U11')]

Finally, a convenience function for determining sources observable from
the full VLBA, from a list of candidate targets,
get\_visible\_targets\_from\_source\_list\_VLBA()

.. code:: python

    obs.get_visible_targets_from_source_list_VLBA(beginobs, endobs, srclist, elevation_limit=15.,
        decbin_limits_deg=[-90,90], minimum_observability_minutes='full', 
        minimum_mutual_vis_observers='full', coord_format='sex', skip_header=None, 
        delimiter=',', nsteps=100)

    # --> [array(['J0137+0923', ' 01:37:51.10', ' +09:23:19.3'], dtype='<U12')]

Naturally, using the same input as for the general multistation version
above, the result is the same three target sources.

A simple declination cut can be enforced by specifying the minimum and
maximum range in degrees with the decbin\_limits\_deg keyword. It
defaults to [-90,90], which returns targets of any declination. If we wanted to restrict the
targets in our test list to being between -10 degrees and +60 degrees
declination, set decbin\_limits\_deg=[-10,60]

.. code:: python

    obs.get_visible_targets_from_source_list_multistation(station_array, beginobs, endobs,  
        './some_sources.txt', elevation_limit=15., decbin_limits_deg=[-10,60], 
        minimum_observability_minutes=30, minimum_mutual_vis_observers='full', coord_format='sex')

    #[array(['J0137+0923', ' 01:37:51.10', ' +09:23:19.3'], dtype='<U12'),
    # array(['J0319+5142', ' 03:19:56.06', ' +51:42:29.0'], dtype='<U12')]
    # --> J0001-1551, at DEC = -15:51:06.7, was cut.

Currently, no advanced optimization features are implemented, such as
clustering by nearest-neighbor searches, or weighting of sources for
preferential selection. These are wishlist todo items which may
hopefully come in the future.
