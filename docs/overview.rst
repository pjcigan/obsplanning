Introduction to working with pyephem and obsplanning objects
============================================================

A brief overview of how the various observer and target objects are
created and used in various computations. Most of obsplanning's inner
workings rely on the pyephem package - in particular, the
telescope/observatory locations, sky targets, and many of the dates &
times are instantiated as ephem Observer(), FixedBody() and Date()
objects.

Let's start by importing various packages used below.

.. code:: python

    import numpy as np
    import ephem
    import datetime as dt
    from astropy.time import Time
    import pytz

    import obsplanning as obs


Observers (telscope sites)
~~~~~~~~~~~~~~~~~~~~~~~~~~

To start things off, here's a look at how telescope or observatory
objects are created in obsplanning. As an example, we will plan a set of
observations of the Crab Nebula from the William Herschel Telescope on
lovely La Palma. The WHT is located at -17:52:53.8 East in longitude,
+28:45:37.7 in latitude, and at an elevation of 2344 meters:

.. code:: python

    # The basic format: 
    # observer_object = obs.create_ephem_observer( Name, longitude[+E], latitude, altitude[m] )

    wht = obs.create_ephem_observer('WHT', '-17 52 53.8', '28 45 37.7', 2344)

    #Equivalent ways to create this object using different coordinate formats:
    wht = obs.create_ephem_observer('WHT', '-17:52:53.8', '28:45:37.7', 2344)
    wht = obs.create_ephem_observer('WHT', -17.88161, 28.760472, 2344)
    wht = obs.create_ephem_observer('WHT', -0.31209297, 0.50196493, 2344, decimal_format='rad') 

String sexagesimal coordinates are automatically parsed by pyephem,
delimited by spaces or colons. Decimal coordinates can also be used.
Pyephem Observer objects internally store these coordinates in radians,
but as degrees are far more commonly used for decimal representation of
RA/DEC coordinates, this is the default for float input in obsplanning.
The keyword decimal\_format can be set to 'radians' for native pyephem
format if desired.

The numerical values in a pyephem object can be accessed in different
ways/formats, depending on how they are called. Simply calling the
object parameter (e.g., wht.lat for the latitude of our WHT observer
above) returns its internal numerical value. In the case of Observer()
objects, the .lat and .lon are stored as radians. If instead the user
calls it inside of print(), it prints the values in a more normal
human-readable format, such as sexagesimal coordinates for Observers and
sky targets.

.. code:: python

    wht.lat,wht.lon          
    #-->  (0.501964934706148, -0.3120929894500905) [in radians]

    print(wht.lat,wht.lon)   
    #-->  28:45:37.7 -17:52:53.8 

Dates
~~~~~

Internal calculations involving dates & times are mainly performed on
ephem.Date() objects. The common datetime.datetime format can also be
used for input, as well as strings formatted as 'YYYY/MM/DD HH:MM:SS.s'
which are automatically parsed by pyephem. Internally, the date info in
an ephem.Date object is stored as a floating point number, which is the
Dublin Julian Date (days since 1899 Dec. 31 at noon). The date/time can
then be returned from the object in a variety of formats.

.. code:: python

    obsdate_string='2025/01/01 23:59:00'
    obsdate_ephem=ephem.Date(obsdate_string)

    obsdate_ephem
    # --> 45657.49930555555  [In Dublin Julian days]

    obsdate_ephem.datetime()
    # --> datetime.datetime(2025, 1, 1, 23, 58, 59, 999999)

    print(obsdate_ephem)  #or, str(obsdate_ephem)
    # --> 2025/1/1 23:59:00

If you want to convert times to Modified Julian Days, it's easy to do so
from datetime, ephem.Date, or string format:

.. code:: python

    obs.MJD(obsdate_ephem)
    # --> 60676.9993055556

    obs.MJD(dt.datetime(2025,1,1,23,59,0))
    # --> 60676.9993055556

    obs.MJD('2025/01/01 23:59:00')
    # --> 60676.9993055556

    ### And similarly for regular Julian Days...
    obs.JD('2025/01/01 23:59:00')
    # --> 2460677.4993055556

Note that all input times are generally assumed to be in UTC. Local time
can be used in many cases with timezone-aware dt.datetime objects, which
is discussed below.

To perform calculations of things that are dependent on time at a
specific location (such as az/alt, position of Sun & Moon in sky from
the observer...), you would typically set the Observer.date before the
body.compute(Observer) calculations are performed.

.. code:: python

    wht.date='2025/01/01 23:59:59'

    ephem.Sun(wht).rise_time
    # --> 45656.83757346579

    print(ephem.Sun(wht).rise_time)
    # --> 2025/1/1 08:06:06

This process is handled automatically by the relevant calculation
functions in obsplanning.

Once an Observer is defined, and a date/time is applied to it, several
useful pieces of information can be computed, such as sunrise and
sunset, various stages of twilight, phases of the moon, etc.

.. code:: python

    ### Some calculations for the Sun:

    # Compute sunset/sunrise, civil twilights, nautical twilights, and astronomical twilights:
    # Twilight definitions: Civil = Sun is 6 degrees below horizon, 
    #                       Nautical = -12 deg, Astronomical = -18 deg

    
    sunset, twi_civil, twi_naut, twi_astro = obs.calculate_twilight_times(wht, 
        '2025/01/01 23:59:00', verbose=True)
    
    #  Sunset :   2025/1/1 18:24:27
    #  Sunrise :  2025/1/2 08:06:20
    #  Twilights
    #  Civil:        previous start at  2025/1/1 18:52:28 ,  next end at  2025/1/2 07:38:20
    #  Nautical:     previous start at  2025/1/1 19:21:09 ,  next end at  2025/1/2 07:09:41
    #  Astronomical: previous start at  2025/1/1 19:49:59 ,  next end at  2025/1/2 06:40:51


    print(sunset)
    # -->  [45657.26698041 45657.83773626]
    # each of the returned objects are list pairs


    ### Some calculations for the Moon:

    moonrise = obs.calculate_moon_times(wht,'2025/01/01 23:59:00',outtype='datetime', 
        verbose=True);
    #  Previous moonrise : 2025/1/1 09:36:32
    #  Next moonset :      2025/1/1 20:03:24

    obs.compute_moonphase(obsdate_string,return_fmt='perc')
    # --> 5.068473554934938  [percent]

    #return_fmt can be 'perc' for percent, 'frac' for fraction, or 'name' for colloquial name
    obs.compute_moonphase('2025/01/01 23:59:00',return_fmt='name')
    # --> 'Waxing Crescent'
    obs.compute_moonphase('2024/12/30 23:59:00',return_fmt='name')
    # --> 'New'


--------------


Local time & timezone-aware datetimes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The dates and times internally stored in ephem.Date objects are not
timezone-aware. That is, input times are assumed to be UTC. Local
timezone information *can* be incorporated for input by using pytz with
datetime.datetime. A standard datetime.datetime object defaults to UTC,
but timezones can be applied with either pytz builtin defined objects or
a standard Olson database timezone name (e.g., 'US/Mountain' or
'America/Chicago'), which will be resolved by pytz.

.. code:: python

    obsstart = '2025/01/01 18:00:00'  # 6PM local time
    obsend   = '2025/02/01 06:00:00'  # 6AM local time

    obsstart_dt_naive = dt.datetime.strptime(obsstart,'%Y/%m/%d %H:%M:%S')
    # or, e.g.:
    obsstart_dt_naive = ephem.Date(obsstart).datetime()
    obsstart_dt_naive = obs.construct_datetime(obsstart,'dt')
    # -->  datetime.datetime(2025, 1, 1, 18, 0)
    # There is no tzinfo there...

    # Apply the timezone, for example to US Pacific time:
    obsstart_local = obs.dt_naive_to_dt_aware(obsstart_dt_naive,'US/Pacific')
    # --> datetime.datetime(2025, 1, 1, 18, 0, tzinfo=<DstTzInfo 'US/Pacific' PST-1 day, 
    #                                                                       16:00:00 STD>)

Calculate time in UTC from a timezone-aware dt object

.. code:: python

    obsstart_utc = obs.local_to_utc(obsstart_local)
    # --> datetime.datetime(2025, 1, 2, 2, 0, tzinfo=<UTC>)

Alternatively, if your datetime object is already timezone-aware (has
tzinfo attached), use the following

.. code:: python

    # If it's tz-aware and the timezone is UTC
    obsstart_local = obs.utc_to_local(obsstart_utc,'US/Pacific')
    # --> datetime.datetime(2025, 1, 1, 18, 0, tzinfo=<DstTzInfo 'US/Pacific' PST-1 day, 
    #                                                                       16:00:00 STD>)

    # Or, to convert from one timezone to another (in this case the local Pacific time 
    # to Eastern time):
    obs.local_to_local(obsstart_local,'US/Eastern')
    # --> datetime.datetime(2025, 1, 1, 21, 0, tzinfo=<DstTzInfo 'US/Eastern' EST-1 day, 
    #                                                                       19:00:00 STD>)

Convert your ephem.Date objects into timezone-aware dt.datetime objects
similarly:

.. code:: python

    obsstart_local_aware = obs.dt_naive_to_dt_aware( ephem.Date(obsstart).datetime() , 
        'Atlantic/Canary' )

Or the reverse: convert from a timezone-aware dt.datetime to ephem.Date,
which does not store timezone information:

.. code:: python

    obs.dtaware_to_ephem(obsstart_local)
    # 45657.583333333336

You can now supply these timezone-aware datetime objects to any
functions that accept datetime format, or convert the tz-aware dt to
ephem.Date format as shown above, if you prefer to supply local times
instead of UTC.



Making Observers timezone-aware
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Observers created using obsplanning.create\_ephem\_observer() actually
use a slightly modified ("decorated" in python parlance) version of the
standard ephem.Observer, now including an additional optional attribute
called timezone. This information is used in several plotting functions
for displaying local time. Timezones can be included on creation as
follows:

.. code:: python

    # Default case, with timezone set to None.
    wht = obs.create_ephem_observer('WHT', '-17:52:53.8', '28:45:37.7', 2344, timezone=None) 

    # Manually specifying the timezone, if it is known by the user.
    # In the case of the WHT used in the above examples, the timezone is 'Atlantic/Canary'
    wht = obs.create_ephem_observer('WHT', '-17:52:53.8', '28:45:37.7', 2344, 
        timezone='Atlantic/Canary')

    # Automatically determine the timezone, using the latitude & longitude
    wht = obs.create_ephem_observer('WHT', '-17:52:53.8', '28:45:37.7', 2344, 
        timezone='calculate')

As seen in the last example above, the timezone can be calculated
automatically from the Observer's coordinates using tzwhere. Users can
determine the timezone for an Observer that is already defined (with or
without timezone already set) with, e.g.:

.. code:: python

    obs.autocalculate_observer_timezone(wht)
    # --> 'Atlantic/Canary'

    # Set the timezone manually after the fact like this:
    wht.timezone='Atlantic/Canary'

The autocalculate\_observer\_timezone() function also works on standard
ephem.Observers that do not have the .timezone attribute. (Though
setting the .timezone attribute will only work for Observers created
with obs.create\_ephem\_observer function or the
obs.Observer\_with\_timezone class. )



Extra timezone utilities:
^^^^^^^^^^^^^^^^^^^^^^^^^

Calculate the UTC offset of a particular timezone at the current time
(also accounts for daylight savings).

.. code:: python

    dt.datetime.utcnow()  
    #--> datetime.datetime(2021, 9, 26, 1, 32, 23, 635631)

    obs.calculate_current_utcoffset('America/Chicago')  
    #--> -5.0

Check e.g. `this Wikipedia
page <https://en.wikipedia.org/wiki/List_of_tz_database_time_zones>`__
for a list of timezone database names for many regions around the world.
You can view the set of all timezones in pytz with
pytz.all\_timezones\_set , or if you know the UTC offset in hours,
obsplanning can give you a list of the pytz timezones there.

.. code:: python

    obs.pytz_timezones_from_utc_offset(-8, common_only=True)

    # ['America/Los_Angeles', 'America/Tijuana', 'America/Vancouver', 'Canada/Pacific',
    #  'Pacific/Pitcairn', 'US/Pacific']

Check if a datetime object is already tz-aware

.. code:: python

    dt_naive = dt.datetime.strptime('2021/10/31 23:59:59','%Y/%m/%d %H:%M:%S')
    obs.is_dt_tzaware(dt_naive) #--> False 
    obs.is_dt_tzaware(dt_naive.replace(tzinfo=pytz.UTC)) #--> True

You can return just the timezone name string of a tz-aware datetime with

.. code:: python

    obsstart_local.tzinfo.zone
    # --> 'US/Pacific'


--------------


Astronomical targets
~~~~~~~~~~~~~~~~~~~~

Astronomical sources or targets are handled internally in obsplanning
with ephem.FixedBody objects. Let's take, for example, the Crab Nebula
or M1, with RA,DEC = [05:34:31.94, 22:00:52.2]. 
Create the target object easily from RA and DEC coordinates like so:

.. code:: python

    crab = obs.create_ephem_target('Crab Nebula','05:34:31.94','22:00:52.2') #'M1'

As with the Observer objects, the internally stored parameters can be
accessed and printed to human-readable formats:

.. code:: python

    crab.name
    # --> Crab Nebula

    crab.ra,crab.dec
    # --> (1.4653339885465477, 0.3844759277966574)  [In radians]

    print(crab.ra,crab.dec)
    # --> 5:35:49.79 22:01:43.9

Once the target and observer objects have been created, the ephemeris
for the target can be computed, after the observer has been updated with
the desired date & time:

.. code:: python

    wht.date='2025/01/01 23:59:59'
    crab.compute(wht)
    print('M1 altitude on %s is %.2f deg'%(wht.date, crab.alt*180./np.pi))
    # M1 altitude on 2025/1/1 23:59:59 is 83.27 deg

A general report of useful ephemeris data can be printed to screen with
obs.ephemeris\_report :

.. code:: python

    obs.ephemeris_report(crab, wht, '2025/01/01 23:59:00')
    #  Target rises at 2025/1/1 17:07:01 with azimuth 64.32 deg, sets at 2025/1/1 06:57:00 
    #                                                               with azimuth 295.68 deg
    #  Target transits at 2025/1/1 00:03:59 with altitude 83.27 deg
    #  Target rises during this night
    #  Target is not circumpolar
    #  For local time of 2025/1/1 23:59:59, sidereal time (LST) is 5:35:59.88



Converting coordinates
~~~~~~~~~~~~~~~~~~~~~~

The equatorial coordinates are already accessible with target.ra and
target.dec, but an ephem.Equatorial class also exists, which is useful
for calculations at a specific epoch:

.. code:: python

    crab_coords_equatorial = ephem.Equatorial(crab, epoch=ephem.J2000)
    print(crab_coords_equatorial.ra,crab_coords_equatorial.dec)
    # 5:34:31.94 22:00:52.2

To convert a target's coordinates from Equatorial (RA/DEC) to Ecliptic
(Lon/Lat) or Galactic (Lon/Lat):

.. code:: python

    # Convert to Ecliptic
    crab_coords_ecliptic = ephem.Ecliptic(crab, epoch=ephem.J2000)
    print(crab_coords_ecliptic.lon,crab_coords_ecliptic.lat)
    # 84:05:51.0 -1:17:40.0

    # Convert to Galactic
    crab_coords_galactic = ephem.Galactic(crab, epoch=ephem.J2000)
    print(crab_coords_galactic.lon,crab_coords_galactic.lat)
    # 184:33:26.8 -5:47:03.7

Coordinate objects like ephem.Equatorial or ephem.Galactic can be made
from the target as shown above, or directly from other coordinate
objects:

.. code:: python

    crab_coords_galactic = ephem.Galactic(crab_coords_ecliptic, epoch=ephem.J2000)
    print(crab_coords_galactic.lon,crab_coords_galactic.lat)
    # 184:33:26.8 -5:47:03.7

NOTE: the specified epoch is interpreted with ephem.Date(), so
'YYYY/MM/DD HH:MM:SS' strings and datetime objects can be used. Floats
are interpreted as the native ephem.Date values of Dublin Julian Days --
so do not use 'epoch=2000' for a J2000 epoch, as it will be interpreted
as 2000 days since the 31st of December, 1899.

.. code:: python

    crab_coords_2000 = ephem.Equatorial(crab, epoch=ephem.J2000)
    crab_coords_2000.epoch
    # --> 36525.0

    print(crab_coords_2000.epoch)
    # --> 2000/1/1 12:00:00

    ### This will give incorrect results:
    crab_coords_wrong = ephem.Equatorial(crab, epoch=2000.)
    print(crab_coords_wrong.epoch)
    # -->  1905/6/23 12:00:00

Calculate coordinates at an arbitrary date, such as epoch=2050:

.. code:: python

    crab_coords_2050 = ephem.Equatorial(crab, epoch='2050/01/01 12:00:00')
    print(' Coordinates at epoch=%s : RA = %s , DEC = %s'%(crab_coords_2050.epoch, 
            crab_coords_2050.ra, crab_coords_2050.dec))
    # Coordinates at epoch=2050/1/1 12:00:00 : RA = 5:37:32.60 , DEC = 22:02:36.8


--------------


Various calculations and tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Building on the basic functionality outlined above, let's explore some
of the tools for producing useful information for observations.

Again, let's take the example of observing the Crab Nebula from the WHT.
This time, we will plan observations for the night of Jan 1, 2025,
starting 30 minutes after sunset and ending 30 minutes before sunrise.

.. code:: python


    wht = obs.create_ephem_observer('WHT', '-17 52 53.8', '28 45 37.7', 2344)
    crab = obs.create_ephem_target('Crab Nebula','05:34:31.94','22:00:52.2') 

    sunset, twi_civil, twi_naut, twi_astro = obs.calculate_twilight_times(wht, 
        '2025/01/01 23:59:00')

    #Sun up/down +/-30min
    obsstart=sunset[0]+30.*ephem.minute
    obsend=sunset[1]-30.*ephem.minute 
    print('Start at %s, end at %s'%(ephem.date(obsstart),ephem.date(obsend)))
    # Start at 2025/1/1 18:54:27, end at 2025/1/2 07:36:20

Calculate the rise, set, and transit times of the target, from the
viewpoint of the specified observatory.

.. code:: python

    # Transit time (when it passes through the meridian / peak altitude)
    crab_transit = obs.calculate_transit_time_single(crab, wht, '2025/01/01 23:59:59', 
        return_fmt='str') 
    # --> '2025/1/2 00:00:03'  [using the default transit mode='nearest']

    # Rise & set times
    crab_RStimes = obs.calculate_rise_set_times_single(crab, wht, '2025/01/1 23:59:59', 
        return_fmt='str') 
    #--> ['2025/1/1 17:07:01', '2025/1/2 06:53:04']

Calculate values for the target's altitude and azimuth over the course
of the observations, from the viewpoint of the observatory. Here the
alt/az values are calculated at 200 intervals between the start and end
times.

.. code:: python

    m1_alts,m1_azs = obs.compute_target_altaz(crab, wht, obsstart, obsend, nsteps=200)

    #Then sec(z) airmass is easily computed from altitudes:
    m1_airmass=obs.alt2airmass(m1_alts)

These can be used to make plots of the target altitude vs time -- the
classic visibility plot. Examples of this are shown in the next tutorial
page.

The altitudes of the moon and Sun can also be calculated in a similar
way:

.. code:: python

    moon_alts,moon_azs = obs.compute_target_altaz( ephem.Moon(), wht, obsstart, obsend, 
        nsteps=200)
    sun_alts,sun_azs = obs.compute_target_altaz( ephem.Sun(), wht, obsstart, obsend, 
        nsteps=200)

Conversion to sidereal time is straightforward:

.. code:: python

    times_sidereal = obs.compute_sidereal_times(wht, obsstart, obsend, nsteps=200) 
    #LST, as ephem.Angle

Calculate the angular separation or distance on the sky from a target
and the Moon, for the specified time.

.. code:: python

    ### Separation from moon on given date:
    moonsep_start = obs.moonsep_single(crab,wht,obsstart)  #140.34 deg
    moonsep_end = obs.moonsep_single(crab,wht,obsend)      #131.74 deg

Separation from Sun is also of particular interest for daytime
observations (e.g., in the radio or submm).

.. code:: python

    sunsep_start = obs.sunsep_single(crab,wht,obsstart)  #162.79 deg

It's also straightforward to calculate the angular separation or
distance from any other fixed sky object. Useful for finding the nearest
flux calibrator, or on the next science target in a list to get a sense
for slew times, etc... In this example, the separation on the sky
between galaxies NGC 1052 and NGC 3079 is calculated.

.. code:: python

    ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751')
    ngc3079=obs.create_ephem_target('NGC3079','10:01:57.80','55:40:47.24') 

    obs.skysep_fixed_single(ngc1052,ngc3079)  #--> 108.13847548432832 [degrees]

This general sky separation function can also be used for separaation
from the Sun/moon, but you would first need to instantiate them with
specified time. The moonsep\_single and sunsep\_single functions are
recommended indead, as they include this step.

If you have a list of, e.g., potential calibrator targets and want to
determine which of them is closest to your science target, this can be
determined easily like in the following example that calculates the
nearest of a set of standard calibrators to NGC 1052.

.. code:: python

    obs.nearest_from_target_list(ngc1052, [obs.SRC_3C84,obs.SRC_3C286,obs.SRC_3C273], 
        verbose=True)
    # Angular separations on sky from NGC1052:
    #         3C84 = 50.55 deg
    #        3C286 = 152.38 deg
    #        3C273 = 146.57 deg

    # --> '3C84'

As seen in the example above, obsplanning has several common radio
calibrator objects pre-defined. Further discussion of radio-oriented
tools in obsplanning are covered in a later tutorial.
