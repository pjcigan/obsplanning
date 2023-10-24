
# Introduction to working with pyephem and obsplanning objects

A brief overview of how the various observer and target objects are created and used in various computations. Most of obsplanning's inner workings rely on the pyephem package - in particular, the telescope/observatory locations, sky targets, and many of the dates & times are instantiated as ephem Observer(), FixedBody() and Date() objects.

Let's start by importing various packages used below.

```python
import numpy as np
import ephem
import datetime as dt
from astropy.time import Time
import pytz

import obsplanning as obs
```


### Observers (telscope sites)

To start things off, here's a look at how telescope or observatory objects are created in obsplanning.  As an example, we will plan a set of observations of the Crab Nebula from the William Herschel Telescope on lovely La Palma.  The WHT is located at -17:52:53.8 East in longitude, +28:45:37.7 in latitude, and at an elevation of 2344 meters:
```python
# The basic format:
# observer_object = obs.create_ephem_observer( Name, longitude[+E], latitude, altitude[meters] )

wht = obs.create_ephem_observer('WHT', '-17 52 53.8', '28 45 37.7', 2344)

#Equivalent ways to create this object using different coordinate formats:
wht = obs.create_ephem_observer('WHT', '-17:52:53.8', '28:45:37.7', 2344)
wht = obs.create_ephem_observer('WHT', -17.88161, 28.760472, 2344)
wht = obs.create_ephem_observer('WHT', -0.31209297, 0.50196493, 2344, decimal_format='rad')
```
String sexagesimal coordinates are automatically parsed by pyephem, delimited by spaces or colons. Decimal coordinates can also be used. Pyephem Observer objects internally store these coordinates in radians, but as degrees are far more commonly used for decimal representation of RA/DEC coordinates, this is the default for float input in obsplanning.  The keyword decimal\_format can be set to 'radians' for native pyephem format if desired.  

The numerical values in a pyephem object can be accessed in different ways/formats, depending on how they are called.  Simply calling the object parameter (e.g., wht.lat for the latitude of our WHT observer above) returns its internal numerical value.  In the case of Observer() objects, the .lat and .lon are stored as radians.  If instead the user calls it inside of print(), it prints the values in a more normal human-readable format, such as sexagesimal coordinates for Observers and sky targets.  
```python
wht.lat,wht.lon          
#-->  (0.501964934706148, -0.3120929894500905) [in radians]

print(wht.lat,wht.lon)   
#-->  28:45:37.7 -17:52:53.8
```



### Dates

Internal calculations involving dates & times are mainly performed on ephem.Date() objects.  The common datetime.datetime format can also be used for input, as well as strings formatted as 'YYYY/MM/DD HH:MM:SS.s' which are automatically parsed by pyephem.  Internally, the date info in an ephem.Date object is stored as a floating point number, which is the Dublin Julian Date (days since 1899 Dec. 31 at noon).  The date/time can then be returned from the object in a variety of formats.

```python
obsdate_string='2025/01/01 23:59:00'
obsdate_ephem=ephem.Date(obsdate_string)

obsdate_ephem
# --> 45657.49930555555  [In Dublin Julian days]

obsdate_ephem.datetime()
# --> datetime.datetime(2025, 1, 1, 23, 58, 59, 999999)

print(obsdate_ephem)  #or, str(obsdate_ephem)
# --> 2025/1/1 23:59:00
```

If you want to convert times to Modified Julian Days, it's easy to do so from datetime, ephem.Date, or string format:
```python
obs.MJD(obsdate_ephem)
# --> 60676.9993055556

obs.MJD(dt.datetime(2025,1,1,23,59,0))
# --> 60676.9993055556

obs.MJD('2025/01/01 23:59:00')
# --> 60676.9993055556

### And similarly for regular Julian Days...
obs.JD('2025/01/01 23:59:00')
# --> 2460677.4993055556
```

Note that all input times are generally assumed to be in UTC.  Local time can be used in many cases with timezone-aware dt.datetime objects, which is discussed below.


To perform calculations of things that are dependent on time at a specific location (such as az/alt, position of Sun & Moon in sky from the observer...), you would typically set the Observer.date before the body.compute(Observer) calculations are performed.
```python
wht.date='2025/01/01 23:59:59'

ephem.Sun(wht).rise_time
# --> 45656.83757346579

print(ephem.Sun(wht).rise_time)
# --> 2025/1/1 08:06:06
```
This process is handled automatically by the relevant calculation functions in obsplanning.  


Once an Observer is defined, and a date/time is applied to it, several useful pieces of information can be computed, such as sunrise and sunset, various stages of twilight, phases of the moon, etc.
```python
### Some calculations for the Sun:

# Compute sunset/sunrise, civil twilights, nautical twilights, and astronomical twilights:
# Twilight definitions: Civil = Sun is 6 degrees below horizon, Nautical = -12 deg, Astronomical = -18 deg

sunset, twi_civil, twi_naut, twi_astro = obs.calculate_twilight_times(wht, '2025/01/01 23:59:00', verbose=True)
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

moonrise = obs.calculate_moon_times(wht,'2025/01/01 23:59:00',outtype='datetime', verbose=True);
#  Previous moonrise : 2025/1/1 09:36:32
#  Next moonset :      2025/1/1 20:03:24

obs.compute_moonphase(obsdate_string,return_fmt='perc')
# --> 5.068473554934938  [percent]

#return_fmt can be 'perc' for percent, 'frac' for fraction, or 'name' for colloquial name
obs.compute_moonphase('2025/01/01 23:59:00',return_fmt='name')
# --> 'Waxing Crescent'
obs.compute_moonphase('2024/12/30 23:59:00',return_fmt='name')
# --> 'New'
```



#### Local time & timezone-aware datetimes

The dates and times internally stored in ephem.Date objects are not timezone-aware.  That is, input times are assumed to be UTC.  Local timezone information _can_ be incorporated for input by using pytz with datetime.datetime. A standard datetime.datetime object defaults to UTC, but timezones can be applied with either pytz builtin defined objects or a standard Olson database timezone name (e.g., 'US/Mountain' or 'America/Chicago'), which will be resolved by pytz.  
```python
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
# --> datetime.datetime(2025, 1, 1, 18, 0, tzinfo=<DstTzInfo 'US/Pacific' PST-1 day, 16:00:00 STD>)
```

Calculate time in UTC from a timezone-aware dt object
```python
obsstart_utc = obs.local_to_utc(obsstart_local)
# --> datetime.datetime(2025, 1, 2, 2, 0, tzinfo=<UTC>)
```

Alternatively, if your datetime object is already timezone-aware (has tzinfo attached), use the following
```python
# If it's tz-aware and the timezone is UTC
obsstart_local = obs.utc_to_local(obsstart_utc,'US/Pacific')
# --> datetime.datetime(2025, 1, 1, 18, 0, tzinfo=<DstTzInfo 'US/Pacific' PST-1 day, 16:00:00 STD>)

# Or, to convert from one timezone to another (in this case the local Pacific time to Eastern time):
obs.local_to_local(obsstart_local,'US/Eastern')
# --> datetime.datetime(2025, 1, 1, 21, 0, tzinfo=<DstTzInfo 'US/Eastern' EST-1 day, 19:00:00 STD>)
```

Convert your ephem.Date objects into timezone-aware dt.datetime objects similarly:
```python
obsstart_local_aware = obs.dt_naive_to_dt_aware( ephem.Date(obsstart).datetime() , 'Atlantic/Canary' )
```

Or the reverse: convert from a timezone-aware dt.datetime to ephem.Date, which does not store timezone information:
```python
obs.dtaware_to_ephem(obsstart_local)
# 45657.583333333336
```


You can now supply these timezone-aware datetime objects to any functions that accept datetime format, or convert the tz-aware dt to ephem.Date format as shown above, if you prefer to supply local times instead of UTC.


#### Making Observers timezone-aware

Observers created using obsplanning.create\_ephem\_observer() actually use a slightly modified ("decorated" in python parlance) version of the standard ephem.Observer, now including an additional optional attribute called timezone.  This information is used in several plotting functions for displaying local time. Timezones can be included on creation as follows:
```python
# Default case, with timezone set to None.
wht = obs.create_ephem_observer('WHT', '-17:52:53.8', '28:45:37.7', 2344, timezone=None)

# Manually specifying the timezone, if it is known by the user.
# In the case of the WHT used in the above examples, the timezone is 'Atlantic/Canary'
wht = obs.create_ephem_observer('WHT', '-17:52:53.8', '28:45:37.7', 2344, timezone='Atlantic/Canary')

# Automatically determine the timezone, using the latitude & longitude
wht = obs.create_ephem_observer('WHT', '-17:52:53.8', '28:45:37.7', 2344, timezone='calculate')
```

As seen in the last example above, the timezone can be calculated automatically from the Observer's coordinates using tzwhere.  Users can determine the timezone for an Observer that is already defined (with or without timezone already set) with, e.g.:
```python
obs.autocalculate_observer_timezone(wht)
# --> 'Atlantic/Canary'

# Set the timezone manually after the fact like this:
wht.timezone='Atlantic/Canary'
```
The autocalculate\_observer\_timezone() function also works on standard ephem.Observers that do not have the .timezone attribute.  (Though setting the .timezone attribute will only work for Observers created with obs.create\_ephem\_observer function or the obs.Observer\_with\_timezone class. )



#### Extra timezone utilities:

Calculate the UTC offset of a particular timezone at the current time (also accounts for daylight savings).
```python
dt.datetime.utcnow()  
#--> datetime.datetime(2021, 9, 26, 1, 32, 23, 635631)

obs.calculate_current_utcoffset('America/Chicago')  
#--> -5.0
```

Check e.g. [this Wikipedia page](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones) for a list of timezone database names for many regions around the world.  You can view the set of all timezones in pytz with  pytz.all\_timezones\_set , or if you know the UTC offset in hours, obsplanning can give you a list of the pytz timezones there.
```python
obs.pytz_timezones_from_utc_offset(-8, common_only=True)

# ['America/Los_Angeles', 'America/Tijuana', 'America/Vancouver', 'Canada/Pacific',
#  'Pacific/Pitcairn', 'US/Pacific']
```


Check if a datetime object is already tz-aware
```python
dt_naive = dt.datetime.strptime('2021/10/31 23:59:59','%Y/%m/%d %H:%M:%S')
obs.is_dt_tzaware(dt_naive) #--> False
obs.is_dt_tzaware(dt_naive.replace(tzinfo=pytz.UTC)) #--> True
```

You can return just the timezone name string of a tz-aware datetime with
```python
obsstart_local.tzinfo.zone
# --> 'US/Pacific'
```



### Astronomical targets

Astronomical sources or targets are handled internally in obsplanning with ephem.FixedBody objects.  Let's take, for example, the Crab Nebula or M1, with RA,DEC = \[05:34:31.94, 22:00:52.2\].  Create the target object easily from RA and DEC coordinates like so:

```python
crab = obs.create_ephem_target('Crab Nebula','05:34:31.94','22:00:52.2') #'M1'

```

As with the Observer objects, the internally stored parameters can be accessed and printed to human-readable formats:
```python
crab.name
# --> Crab Nebula

crab.ra,crab.dec
# --> (1.4653339885465477, 0.3844759277966574)  [In radians]

print(crab.ra,crab.dec)
# --> 5:35:49.79 22:01:43.9
```



Once the target and observer objects have been created, the ephemeris for the target can be computed, after the observer has been updated with the desired date & time:
```python
wht.date='2025/01/01 23:59:59'
crab.compute(wht)
print('M1 altitude on %s is %.2f deg'%(wht.date, crab.alt*180./np.pi))
# M1 altitude on 2025/1/1 23:59:59 is 83.27 deg

```

A general report of useful ephemeris data can be printed to screen with obs.ephemeris\_report :
```python
obs.ephemeris_report(crab, wht, '2025/01/01 23:59:00')
#  Target rises at 2025/1/1 17:07:01 with azimuth 64.32 deg, sets at 2025/1/1 06:57:00 with azimuth 295.68 deg
#  Target transits at 2025/1/1 00:03:59 with altitude 83.27 deg
#  Target rises during this night
#  Target is not circumpolar
#  For local time of 2025/1/1 23:59:59, sidereal time (LST) is 5:35:59.88
```


#### Querying Coordinates

What if you don't know your source coordinates, or you want to check against 'standard' catalog positions?  obsplanning has a convenience function ```query_object_coords_simbad()``` to quickly do just that by querying the Simbad service over the internet.  By default, the queried coordinates are returned in decimal format (for use in calculations), but you can specify to return them as their native sexagesimal strings.

```python
# Look up the Simbad coordinates for the Crab Nebula
# -- by default the coordinates are returned in ICRS decimal format.

obs.query_object_coords_simbad('M1') #, return_fmt='dec')
# --> [83.62875, 22.014722222222222]

# You can also specify to return the coordinates in string sexagesimal format,
# which is useful for printing or comparing to tables, etc.

obs.query_object_coords_simbad('M1', return_fmt='sex')
# --> ['05 34 30.9', '+22 00 53']

# Example using keyword args that get passed to Simbad.query_object :
# query objects named m1 through m9 using wildcard, print out (verbose)
# details to terminal, and return 3rd entry from the resulting table
# (use_entry = 2, since python is zero-indexed)

query_object_coords_simbad("m [1-9]", wildcard=True, verbose=True, use_entry=2)

# MAIN_ID       RA         DEC     RA_PREC DEC_PREC ... COO_QUAL COO_WAVELENGTH     COO_BIBCODE     SCRIPT_NUMBER_ID
#            "h:m:s"     "d:m:s"                    ...                                                             
#--------- ----------- ----------- ------- -------- ... -------- -------------- ------------------- ----------------
#    M   1  05 34 30.9   +22 00 53       5        5 ...        E              R 1995AuJPh..48..143S                1
#    M   2 21 33 27.02 -00 49 23.7       6        6 ...        D              O 2010AJ....140.1830G                1
#    M   3 13 42 11.62 +28 22 38.2       6        6 ...        D              O 2010AJ....140.1830G                1
#    M   4 16 23 35.22 -26 31 32.7       6        6 ...        D              O 2010AJ....140.1830G                1
#    M   5 15 18 33.22 +02 04 51.7       6        6 ...        D              O 2010AJ....140.1830G                1
#NGC  6405  17 40 16.6   -32 14 31       5        5 ...        D              O 2021A&A...647A..19T                1
#NGC  6475  17 53 47.3   -34 50 28       5        5 ...        D              O 2021A&A...647A..19T                1
#    M   8    18 03 37    -24 23.2       4        4 ...        E                                                   1
#    M   9 17 19 11.78 -18 30 58.5       6        6 ...        D                2002MNRAS.332..441F                1
## --> Output:
# [205.54841666666667, 28.377277777777778]
```

At the moment, only Simbad queries by target name are implemented here; in the future other query services will be added. Under the hood, this is simply calling ```astroquery.simbad.Simbad.query_object(stringname, **kwargs)```, so you can alternatively call astroquery manually with other services like NED or many others of your choice.  



##### Converting Coordinates


The equatorial coordinates are already accessible with target.ra and target.dec (or target.g_ra, target.g_dec), but an ephem.Equatorial class also exists, which is useful for calculations at a specific epoch:
```python
crab_coords_equatorial = ephem.Equatorial(crab, epoch=ephem.J2000)
print(crab_coords_equatorial.ra,crab_coords_equatorial.dec)
# 5:34:31.94 22:00:52.2

```

Note that the 'astrometric' RA and DEC coordinates (not dependent on epoch) should be accessed with target.a_ra and target.a_dec, instead.  See the [pyephem documentation](https://rhodesmill.org/pyephem/radec) for more details.


To convert a target's coordinates from Equatorial (RA/DEC) to Ecliptic (Lon/Lat) or Galactic (Lon/Lat):
```python
# Convert to Ecliptic
crab_coords_ecliptic = ephem.Ecliptic(crab, epoch=ephem.J2000)
print(crab_coords_ecliptic.lon,crab_coords_ecliptic.lat)
# 84:05:51.0 -1:17:40.0

# Convert to Galactic
crab_coords_galactic = ephem.Galactic(crab, epoch=ephem.J2000)
print(crab_coords_galactic.lon,crab_coords_galactic.lat)
# 184:33:26.8 -5:47:03.7

```

Coordinate objects like ephem.Equatorial or ephem.Galactic can be made from the target as shown above, or directly from other coordinate objects:

```python
crab_coords_galactic = ephem.Galactic(crab_coords_ecliptic, epoch=ephem.J2000)
print(crab_coords_galactic.lon,crab_coords_galactic.lat)
# 184:33:26.8 -5:47:03.7

```


NOTE: the specified epoch is interpreted with ephem.Date(), so 'YYYY/MM/DD HH:MM:SS' strings and datetime objects can be used.  Floats are interpreted as the native ephem.Date values of Dublin Julian Days -- so do not use 'epoch=2000' for a J2000 epoch, as it will be interpreted as 2000 days since the 31st of December, 1899.

```python
crab_coords_2000 = ephem.Equatorial(crab, epoch=ephem.J2000)
crab_coords_2000.epoch
# --> 36525.0

print(crab_coords_2000.epoch)
# --> 2000/1/1 12:00:00

### This will give incorrect results:
crab_coords_wrong = ephem.Equatorial(crab, epoch=2000.)
print(crab_coords_wrong.epoch)
# -->  1905/6/23 12:00:00

```

Calculate coordinates at an arbitrary date, such as epoch=2050:
```python
crab_coords_2050 = ephem.Equatorial(crab, epoch='2050/01/01 12:00:00')
print(' Coordinates at epoch=%s : RA = %s , DEC = %s'%(crab_coords_2050.epoch,
        crab_coords_2050.ra, crab_coords_2050.dec))
# Coordinates at epoch=2050/1/1 12:00:00 : RA = 5:37:32.60 , DEC = 22:02:36.8

```


Obsplanning has a convenience function to output the RA,DEC coordinates from an ephem target source object, called ```eph2c()```.  This is useful for including compactly in functions to print or convert output.  Options for output style are sexagesimal (calling ```obs.dec2sex``` under the hood, and optionally taking its keyword args), degrees, hours (for RA), radians, or astropy.coordinates.SkyCoord.  
```python
obs.eph2c(crab) #,style='sex')
# ([5, 34, 31.93999999999869], [22, 0, 52.20000000000624])
obs.eph2c(crab, style='sex', as_string=True, decimal_places=3)
# ['05:34:31.940', '22:00:52.200']
obs.eph2c(crab, style='deg')
# [83.63308333333333, 22.0145]
obs.eph2c(crab, style='hour')
# [5.575538888888889, 22.0145]
obs.eph2c(crab, style='rad')
# [1.4596726677614609, 0.3842255081802917]
obs.eph2c(crab, style='skycoord')
# <SkyCoord (ICRS): (ra, dec) in deg
#    (83.63308333, 22.0145)>
```
By default this returns the astrometric/absolute coordinates stored in the object as target.a_ra and target.a_dec.  However, there is an option to return apparent coords (target.ra, target.dec) by setting ```apparent=True```, if the ephem object ephemeris info has already been updated with .compute()




### Various calculations and tools

Building on the basic functionality outlined above, let's explore some of the tools for producing useful information for observations.

Again, let's take the example of observing the Crab Nebula from the WHT.  This time, we will plan observations for the night of Jan 1, 2025, starting 30 minutes after sunset and ending 30 minutes before sunrise.
```python

wht = obs.create_ephem_observer('WHT', '-17 52 53.8', '28 45 37.7', 2344)
crab = obs.create_ephem_target('Crab Nebula','05:34:31.94','22:00:52.2')

sunset, twi_civil, twi_naut, twi_astro = obs.calculate_twilight_times(wht, '2025/01/01 23:59:00')

#Sun up/down +/-30min
obsstart=sunset[0]+30.*ephem.minute
obsend=sunset[1]-30.*ephem.minute
print('Start at %s, end at %s'%(ephem.date(obsstart),ephem.date(obsend)))
# Start at 2025/1/1 18:54:27, end at 2025/1/2 07:36:20
```

##### Source Elevation (Rise, Set, Transit times)

Calculate the rise, set, and transit times of the target, from the viewpoint of the specified observatory.
```python
# Transit time (when it passes through the meridian / peak altitude)
crab_transit = obs.calculate_transit_time_single(crab, wht, '2025/01/01 23:59:59', return_fmt='str')
# --> '2025/1/2 00:00:03'  [using the default transit mode='nearest']

# Rise & set times
crab_RStimes = obs.calculate_rise_set_times_single(crab, wht, '2025/01/1 23:59:59', return_fmt='str')
#--> ['2025/1/1 17:07:01', '2025/1/2 06:53:04']
```

Calculate values for the target's altitude and azimuth over the course of the observations, from the viewpoint of the observatory.  Here the alt/az values are calculated at 200 intervals between the start and end times.
```python
m1_alts,m1_azs = obs.compute_target_altaz(crab, wht, obsstart, obsend, nsteps=200)

#Then sec(z) airmass is easily computed from altitudes:
m1_airmass=obs.alt2airmass(m1_alts)
```
These can be used to make plots of the target altitude vs time -- the classic visibility plot.  Examples of this are shown in the next tutorial page.   


The altitudes of the moon and Sun can also be calculated in a similar way:
```python
moon_alts,moon_azs = obs.compute_target_altaz( ephem.Moon(), wht, obsstart, obsend, nsteps=200)
sun_alts,sun_azs = obs.compute_target_altaz( ephem.Sun(), wht, obsstart, obsend, nsteps=200)
```

Conversion to sidereal time is straightforward:
```python
times_sidereal = obs.compute_sidereal_times(wht, obsstart, obsend, nsteps=200) #LST, as ephem.Angle
```


##### Separations Between Sources

Calculate the angular separation or distance on the sky from a target and the Moon, for the specified time.
```python
### Separation from moon on given date:
moonsep_start = obs.moonsep_single(crab,wht,obsstart)  #140.34 deg
moonsep_end = obs.moonsep_single(crab,wht,obsend)      #131.74 deg
```

Separation from Sun is also of particular interest for daytime observations (e.g., in the radio or submm).
```python
sunsep_start = obs.sunsep_single(crab,wht,obsstart)  #162.79 deg
```

There are convenience functions to print daily Sun and Moon separations for your target, for
a specific time every N days.  This can be useful to help determine optimal observing days in a month or other time period.  Here is an example of printing the separations every 7 days for the month of 2023 October.  
```python
# Every 7 days between 10/1 and 10/31: October 1,8,15,22,29
# Moon separations, every 7 days at midnight:
daily_moonseps(crab, wht, '2023/10/01 00:00:00', '2023/10/31 00:00:00', every_N_days=7)
#Crab Nebula
#  On 2023/10/1 00:00:00, Moon separation = 54.9 deg
#  On 2023/10/8 00:00:00, Moon separation = 37.1 deg
#  On 2023/10/15 00:00:00, Moon separation = 119.7 deg
#  On 2023/10/22 00:00:00, Moon separation = 148.2 deg
#  On 2023/10/29 00:00:00, Moon separation = 47.1 deg
## np.array([ 54.85133653,  37.10377776, 119.71043295, 148.18651131,
#        47.05053793])

# Separations from the Sun at noon, every 14 days, observing from the MK VLBA
# station (which is built in to obsplanning)
daily_sunseps(crab, obs.vlbaMK, '2023/10/01 12:00:00', '2023/10/31 12:00:00', every_N_days=14)
#Crab Nebula
#  On 2023/10/1 12:00:00, Sun separation = 103.6 deg
#  On 2023/10/15 12:00:00, Sun separation = 117.4 deg
#  On 2023/10/29 12:00:00, Sun separation = 131.4 deg
##
# np.array([103.62409241, 117.43859323, 131.36060867])


```


It's also straightforward to calculate the angular separation or distance from any other fixed sky object.  Useful for finding the nearest flux calibrator, or on the next science target in a list to get a sense for slew times, etc...  In this example, the separation on the sky between galaxies NGC 1052 and NGC 3079 is calculated.
```python
ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751')
ngc3079=obs.create_ephem_target('NGC3079','10:01:57.80','55:40:47.24')

obs.skysep_fixed_single(ngc1052,ngc3079)  #--> 108.13847548432832 [degrees]
```
This general sky separation function can also be used for separation from the Sun/moon, but you would first need to instantiate those objects with a specified time. The moonsep\_single and sunsep\_single functions are recommended instead, as they include this step.

It's possible to return the component longitude and latitude values as well as the total separation, if those are needed.  In this case, you also need to specify the frame in which you want them calculated: default 'equatorial' (for RA,DEC), 'galactic' (for l,b), or 'ecliptic' (for lon,lat).  You also need to specify whether you want the components calculated as 'cartesian' (longitude separation following lines of constant latitude) which is the default (and also what is given when returncomponents=True), or calculated as 'spherical' -- which gives the components in spherically orthogonal directions (i.e. the 'longitude' will not be along lines of constant latitude).  
\[Note that the total angular separation will be the same regardless of the frame.\]
```python
src1 = obs.create_ephem_target('src1','01:00:00.0','-30:00:00.0') #[15.0,-30.0] in decimal
src2 = obs.create_ephem_target('src2','23:00:00.0','-30:00:00.0') #[345.0,-30.0] in decimal

obs.skysep_fixed_single(src1,src2)
#--> 25.906047546458453     (smaller than 15+15=30, because of the cos(DEC) term)

obs.skysep_fixed_single(src1,src2, returncomponents='cartesian', componentframe='equatorial')  
#--> (25.906049857216924, -25.905079284444753, 0.0)
#    (total separation, d_RA component, d_DEC component)

obs.skysep_fixed_single(src1,src2, returncomponents='cartesian', componentframe='galactic')
#--> (25.906049857216924, 39.67849268770435, 21.143725639068673)
#    (total, d_l, d_b)      
#   ==> though the individual components seem high, look at their Galactic coords:
print( np.degrees([ephem.Galactic(src1).lon,ephem.Galactic(src1).lat]) )
print( np.degrees([ephem.Galactic(src2).lon,ephem.Galactic(src2).lat]) )
#   [270.22743381 -86.56783073]
#   [ 19.60461978 -65.4241051 ]

obs.skysep_fixed_single(src1,src2, returncomponents='cartesian', componentframe='ecliptic')
#--> (25.906049857216924, -24.381436465113175, 11.53326328663634)
#    (total, d_lon, d_lat)
print( np.degrees([ephem.Ecliptic(src1).lon,ephem.Ecliptic(src1).lat]) )
print( np.degrees([ephem.Ecliptic(src2).lon,ephem.Ecliptic(src2).lat]) )
#   [  0.4629637  -33.22308652]
#   [334.19174077 -21.68982324]

```
There are also helper functions ```vincenty_sphere()``` and ```angulardistance()``` for computing separations with inputs given as simple floats rather than ephem objects.  See their API entries for more details.  


##### Searching for sources within a radius

There are two simple cone search functions in obsplanning, for determining which sources from an array fall within a specified separation radius on the sky.  ```sources_within_radius()``` will take a list (or tuple, array...) of ephem targets, and a reference position (which can be an ephem target, an astropy SkyCoord, or a list of [RA,DEC] values that can be parsed by ephem).  The reference position can be an ephem.Sun or Moon instance, as long as it has been instantiated and updated with time and Observer using .compute(). You can request output as the list of ephem sources themselves, just the names, the indices for the valid sources, or the numpy mask that would return the valid list -- which can be useful in scripting.  Here is an example for a handful of bright Messier sources in the Virgo cluster.
```python
# Let obsplanning automatically query coordinates for some targets
virgo_targets = [obs.create_ephem_target(n,*obs.query_object_coords_simbad(n))
             for n in ['m87', 'm85', 'm60', 'm49', 'm90', 'm98']]

for t in virgo_targets:
    print('  %4s : [ %s, %s]'%(t.name, t.a_ra, t.a_dec) )
#   m87 : [ 12:30:49.42, 12:23:28.0]
#   m85 : [ 12:25:24.05, 18:11:27.9]
#   m60 : [ 12:43:39.97, 11:33:09.7]
#   m49 : [ 12:29:46.80, 8:00:01.0]
#   m90 : [ 12:36:49.80, 13:09:46.5]
#   m98 : [ 12:13:48.29, 14:54:02.0]

# Which of these sources are within 3 degrees of the Virgo cluster center?
obs.sources_within_radius(virgo_targets, ['12:27:00','12:43:00'], 3.0,
    direction='inside', return_format='names')
#--> np.array(['m87', 'm90'])

# Which of these Virgo sources are further than 4 degrees from M87 (which is the
#  zeroth-index of virgo_targets)?  Use direction='outside'
obs.sources_within_radius(virgo_targets[1:], virgo_targets[0], 4.0,
    direction='outside', return_format='names')
#--> np.array(['m85', 'm49', 'm98'])

## The first example again (M87 and M90 within 3 deg of the center), but now
#  instead of just the names, return the indices, the numpy mask, and the array
#  of ephem targets themseleves
obs.sources_within_radius(virgo_targets, ['12:27:00','12:43:00'], 3.0,
    direction='inside', return_format='ind')
#--> np.array([0, 4])

obs.sources_within_radius(virgo_targets, ['12:27:00','12:43:00'], 3.0,
    direction='inside', return_format='mask')
#--> np.array([ True, False, False, False,  True, False])

obs.sources_within_radius(virgo_targets, ['12:27:00','12:43:00'], 3.0,
    direction='inside', return_format='targets')
#--> np.array([<ephem.FixedBody 'm87' at 0x7ff64d50d870>,
#              <ephem.FixedBody 'm90' at 0x7ff64d50ddf0>], dtype=object)
```

The other function for testing separations from an array of sources is ```srctable_within_radius()```.  This function is intended to work with more typical catalogs of sources, and so the input it takes is a pandas or astropy table (or similar recarray-style object that can be referenced by column names) of source coordinates.  Here is a quick example using the same Virgo galaxies above with pandas:
```python
import pandas as import pd

# generate a simple pandas table with columns ['name','RA','DEC'] from the
# ephem target list above.  It's possible to set all columns in one go initially,
# then set dtypes later (because of mixing string and float dtypes in the list
# comprehension), but simpler just to initialize it with the names column first,
# then set the RA,DEC float columns after.

virgo_table = pd.DataFrame([t.name for t in virgo_targets], columns=['name',])
virgo_table[['RA','DEC']] = [[t.ra*180/np.pi, t.dec*180/np.pi] for t in virgo_targets]

print(virgo_table)
#  name          RA        DEC
#0  m87  188.000690  12.261681
#1  m85  186.644222  18.060981
#2  m60  191.209931  11.424425
#3  m49  187.741156   7.870973
#4  m90  189.501168  13.033926
#5  m98  183.747979  14.770090

# Once again, find sources that are within 3 degrees of the cluster center.  Now,
# specify the column names for the RA and DEC values.
obs.srctable_within_radius(virgo_table, ['12:27:00','12:43:00'], 3.0,
    RAlabel='RA', DEClabel='DEC', direction='inside', return_format='sources')
#  name          RA        DEC
#0  m87  188.000690  12.261681
#4  m90  189.501168  13.033926
##--> Now the output is a table slice

## Other equivalent formats that would work for the reference position coordinates:

#float degrees  [186.75, 12.716666666666667]
ref_crd = [186.75, 12.716666666666667]
obs.srctable_within_radius(virgo_table, ref_crd, 3.0)

#SkyCoord       <SkyCoord (ICRS): (ra, dec) in deg
#                   (12.45, 12.71666667)>
ref_crd = obs.coordinates.SkyCoord('12:27:00', '12:43:00', unit='deg')
obs.srctable_within_radius(virgo_table, ref_crd, 3.0)
```
This function may be useful for simple selections of sources from lists such as the Gaia or ICRF3 catalogs.


##### Nearest source from a list

If you have a list of, e.g., potential calibrator targets and want to determine which of them is closest to your science target, this can be determined easily like in the following example that calculates the nearest of a set of standard calibrators to NGC 1052.  
```python
obs.nearest_from_target_list(ngc1052, [obs.SRC_3C84,obs.SRC_3C286,obs.SRC_3C273], verbose=True)
# Angular separations on sky from NGC1052:
#         3C84 = 50.55 deg
#        3C286 = 152.38 deg
#        3C273 = 146.57 deg

# --> '3C84'
```
As seen in the example above, obsplanning has several common radio calibrator objects pre-defined.  Further discussion of radio-oriented tools in obsplanning are covered in a later tutorial.  



##### Calculating Optimal Slew Ordering

Telescopes take time to change their pointing on the sky, and for those that have slow motors, this can be a significant addition to your overhead.  If you need to observe several sources over the course of one observing session, taking care to observe them in a speedy order can save you valuable time that can be better used integrating on your science targets.  This compounds if you need to make several passes at targets in a loop. You can use obsplanning to determine the optimal order in which to observe a list of targets, to make efficient use of your telescope time.  

Let's look at an example of observing a handful of Messier objects in a similar RA range.  Let's use the automatic simbad identifier query function to get the coordinates. (Though for your own observations you probably have much more precise values!)
```python
messier_list = ['m51', 'm101', 'm102', 'm104', 'm87', 'm82']
messier_targets = [obs.create_ephem_target(n,*obs.query_object_coords_simbad(n))
                 for n in messier_list]

for t in messier_targets:
  print('  %4s : [ %s, %s]'%(t.name, t.a_ra, t.a_dec) )
#   m51 : [ 13:29:52.70, 47:11:42.9]
#  m101 : [ 14:03:12.58, 54:20:55.5]
#  m102 : [ 15:06:29.56, 55:45:47.9]
#  m104 : [ 12:39:59.43, -11:37:23.0]
#   m87 : [ 12:30:49.42, 12:23:28.0]
#   m82 : [ 9:55:52.43, 69:40:46.9]

##Using the built-in ephem target .a_ra,.a_dec printout as shown here is
#  convenient, though manual formatting may be required to display more
#  precision, etc.  
```


The function ```calc_optimal_slew_loop()``` will take a list of ephem targets, and calculate the optimal slew order, with a variety of options. Setting ```verbose=True``` will print out every (unique) permutation if you care to see that nitty gritty detail.  Setting ```sortloops=True``` will sort those printed verbose outputs by increasing cumulative slew.

You can force it to start from a particular target with ```set_first=<target>```, which may be useful for bookending with calibrators, or a source that sets early.  By default it returns the list of names as strings, but you can have it return the list of cumulative slew separation angles/times or the list of ephem targets directly with the ```return_format``` keyword.  One assumption for the optimization is that targets will be observed in a repeated loop, and so the separation between the permutation's last and first entries are also considered.  To calculate just a single pass, without looping, you can set ```repeat_loop=False```.  

```python
# Optimize purely by angular separation (telescope slew motor speeds not considered)

obs.calc_optimal_slew_loop(messier_targets, verbose=True, sortloops=True)

#Permutations (repeating the loop)
#  ['m51', 'm101', 'm102', 'm82', 'm87', 'm104']: cumulative distance = 199.4 deg
#  ['m51', 'm104', 'm87', 'm82', 'm102', 'm101']: cumulative distance = 199.4 deg
#  ['m51', 'm101', 'm102', 'm82', 'm104', 'm87']: cumulative distance = 200.4 deg
# ...
#  ['m51', 'm102', 'm87', 'm101', 'm104', 'm82']: cumulative distance = 304.2 deg

## --> returns :
#['m51', 'm101', 'm102', 'm82', 'm87', 'm104']


#Best order for a single pass, without looping back around to the first source
obs.calc_optimal_slew_loop(messier_targets, verbose=True, sortloops=True,
    repeat_loop=False)

# Permutations (single pass)
#  ['m104', 'm87', 'm51', 'm101', 'm102', 'm82']: cumulative distance = 114.3 deg
#  ['m82', 'm102', 'm101', 'm51', 'm87', 'm104']: cumulative distance = 114.3 deg
#  ['m82', 'm101', 'm102', 'm51', 'm87', 'm104']: cumulative distance = 118.4 deg
# ...
```
The optimal slew order from this group of sources has a cumulative value of 199.4 degrees, while the anti-optimal order takes 304.2 degrees of slewing.  This difference of ~100 degrees may not be much if you are going through them in one pass, but if you need to do multiple loops, that can add up.  And if your sources are spread farther across the sky, this can result in significant differences in time.

Note that this algorithm is not optimized for large N, it can be extremely slow for large numbers of input sources!  (Since number of permutations increases rapidly with N...  Groups of 7 or fewer should be almost instantaneous, with dramatic compute requirements above that.)

Wrapped versions of permutations are automatically dropped from the verbose output by default. For example, if [A,B,C] exists, then would drop [B,C,A] and [C,A,B] ( but not [A,C,B] ).  This may make it seem like all the returned permutations are starting from the same target.  But you can turn that off to see ALL the permutations, including wrapped duplicates, if ```drop_wrap_repeats=False```.  

If you know the nominal telescope motor slew rate (or rates, for both Az,El axes), you can input those here to estimate the actual amount of time these loops will take.  Specify that you want to optimize by time with ```optimize_by='time'```, and set the actual slew speeds in degrees/minute with ```AZ_deg_min``` and ```EL_deg_min```.  Some telescopes may have different slew rates along the different axes, and this will affect the optimal ordering for targets.
```python
#Here, perform calculations for the same group of targets, but now specifying
# slew speeds of 90 deg/min in azimuth, and 30 deg/min in elevation (the nominal
# rates for the VLBA dishes).
obs.calc_optimal_slew_loop(messier_targets, verbose=True, sortloops=True,
    optimize_by='time', AZ_deg_min=90., EL_deg_min=30.)

#Permutations (repeating the loop)
#  ['m51', 'm104', 'm87', 'm102', 'm82', 'm101']: 5.42 min.  (211.8 deg)
#  ['m51', 'm101', 'm82', 'm102', 'm87', 'm104']: 5.42 min.  (211.8 deg)
#  ['m51', 'm87', 'm104', 'm102', 'm82', 'm101']: 5.42 min.  (210.1 deg)
# ...
#  ['m51', 'm82', 'm87', 'm101', 'm104', 'm102']: 8.79 min.  (301.5 deg)
#  ['m51', 'm82', 'm87', 'm102', 'm104', 'm101']: 8.79 min.  (300.2 deg)

```
For this particular group of targets and slew rates, optimizing can gain you over 3 minutes (per loop), or 30 seconds of extra integration time on source.
