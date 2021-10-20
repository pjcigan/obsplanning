### ObsPlanning
### v1.0
### written by Phil Cigan
__author__ = "Phil Cigan <pcigan@gmu.edu>"
__version__ = "1.0.1"


"""
Functions to aid in planning and plotting observations.
e.g., calculations for:
    - UTC, GMST, LST, LMST... 
    - Moon/Sun Separation
    - Parallactic Angle
    - Altitude/Airmass of target
    - etc...

Some online observing tools: 
Text data tools  http://www.briancasey.org/artifacts/astro/
JSkyCalc         http://www.dartmouth.edu/~physics/labs/skycalc/flyer.html
python obstools  https://pythonhosted.org/Astropysics/coremods/obstools.html#astropysics.obstools.Site.nightPlot
visplot          http://www.not.iac.es/observing/forms/visplot/




"""

### Import generally useful packages

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
#import time
import datetime as dt

hourformat = mdates.DateFormatter('%H:%M')


##############################################################################
# Import the packages necessary for finding coordinates and making
# coordinate transformations

import astropy.units as u
from astropy.time import Time
#from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle, Latitude, Longitude
#from astropy.coordinates import get_sun, get_moon
import astropy.io.fits as pyfits
import ephem
import pytz
from tzwhere import tzwhere
#import pygeodesy

#For making finder plots:
import os
from astroquery.skyview import SkyView
from astroquery.sdss import SDSS
from astroquery.simbad import Simbad
from astropy import coordinates
from astropy.wcs import WCS 
import astropy.units as u
import multicolorfits as mcf
from matplotlib import rcParams
import matplotlib.patches as patches
import matplotlib.patheffects as PathEffects

from scipy import interpolate #--> Currently only needed for calculate_antenna_visibility_limits()
from tqdm import tqdm


### numpy, as of vers. 19.0, raises a VisibleDeprecationWarning if a function creates an array from
#   'ragged'/non-uniform sequences (e.g., np.array([ [1,2],[3,4,5] ]) ), and warns that you must specify
#   dtype='object'.  
# Some of the imported packages use this behavior, so is not possible to to this on the fly.  
# I am opting here to suppress this particular warning type, to avoid these repeated warnings that do 
# not currently impact the functionality of obsplanning.
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


##############################################################################


##### ------------ General utility functions: angle conversions, etc ------------ #####

def alt2airmass(altitude):
    """Secant(z) = 1/cos(z). 
    z = 90-altitude(elevation)   
    
    Parameters
    ----------
    altitude: float (or numpy.array of floats)
        The altitude/elevation in degrees
    
    Returns
    -------
    airmass : float (or numpy.array of floats)
        The airmass at the specified altitude
    """
    return 1./np.cos( (90-altitude)*np.pi/180 )

def airmass2alt(airmass):
    """Secant(z) = 1/cos(z). 
    z = 90-altitude(elevation)   
    
    Parameters
    ----------
    airmass: float (or numpy.array of floats)
        Airmass along the telescope pointing LOS
    
    Returns
    -------
    altitude : float (or numpy.array of floats)
        The altitude/elevation in degrees corresponding to the airmass
    """
    return 90.-(np.arccos(1./airmass)*180./np.pi)

def wrap_24hr(timein,component=False):
    """
    Wraps a time (float) in hours to a 24-hour time.  Can be used, e.g., to 
    determine the hour of day after adding/subtracting some number of hours.
    
    Parameters
    ----------
    timein: float
        Decimal number of hours
    component : bool
        If True, input timein in H:M:S components will be converted to decimal
    
    Returns
    -------
    timeout : float
        The number of hours on a 0--24 hour scale
    
    Examples
    --------
    # e.g.: 49.5 hours is 1.5 hours (after 2x24); -2 corresponds to 22.0 \n
    obs.wrap_24hr(49.5) #--> 1.5 \n
    obs.wrap_24hr(-2)   #--> 22.0 \n
    obs.wrap_24hr(some_amount_of_seconds/3600.) \n
    obs.wrap_24hr(LMST_reference_hour+offset_hours) \n
    obs.wrap_24hr('49:30:00',component=True)  #--> 1.5
    """
    if component==True: timein=dms2deg(timein) #Not hour2deg, as we don't want to convert to full 360 deg range
    if (timein<0.0): n_hours=int(timein/24.)-1; timeout=timein-n_hours*24. #Equivalent to timein % 24.
    elif (timein>=24.0): n_hours=int(timein/24.); timeout=timein-n_hours*24.
    else: timeout=timein
    return timeout

def wrap_360(valin):
    """
    Wraps a value (float) to a 0--360 degree range.  
    
    Parameters
    ----------
    valin: float
        Input value in degrees
    
    Returns
    -------
    valout : float
    
    Example
    -------
    # e.g., 370 degrees is equivalent to 10 degrees.  \n
    obs.wrap_360(370)  #--> 10.
    
    Note
    ----
    Essentially equivalent to calculating as  (valin % 360) .
    """
    if (valin<0.0): n_revs=int(valin/360.)-1; valout=valin-n_revs*360.
    elif (valin>=360.0): n_revs=int(valin/360.); valout=valin-n_revs*360.
    else: valout=valin
    return valout

def wrap_pm180(valin):
    """
    Wraps a value (float) to a -180 to 180 degree range.  
    
    Parameters
    ----------
    valin: float
        Input value in degrees
    
    Returns
    -------
    valout : float
    
    Example
    -------
    # e.g., 200 degrees corresponds to -160 degrees when limited to [-180,180] \n  
    obs.wrap_pm180(200)  #--> -160.0
    
    Note
    ----
    Similar to (valin % 180.).
    """
    val_red=valin%360-360
    valout=val_red+[360. if val_red<=-180 else 0.][0]
    return valout

def deg2hour(valin):
    """
    Converts decimal degrees to a HMS list of [ Hours (int), Minutes (int), Seconds (decimal)].
    
    Parameters
    ----------
    valin: float
        Input value in degrees
    
    Returns
    -------
    HMS : list
        [ Hours (int), Minutes (int), Seconds (decimal)]
    
    Example
    -------
    # e.g., 180 degrees corresponds to 12 hours, 0 min, 0 sec \n
    obs.deg2hour(180.0)  #--> [12, 0, 0]  
    """
    rmins,rsec=divmod(24./360*valin*3600,60)
    rh,rmins=divmod(rmins,60)
    return [int(rh),int(rmins),rsec]

def hour2deg(valin):
    """
    Converts hours or HMS input to decimal degrees. 
    
    Parameters
    ----------
    valin: float
        Input value in HMS. Can be either: \n
        - a string delimeted by : or spaces \n
        - a list of [H,M,S] numbers (floats or ints)
    
    Returns
    -------
    valout : float
        Degrees corresponding to the HMS value
    
    Examples
    --------
    # e.g., '09:30:15' corresponds to 142.5625 deg \n
    obs.hour2deg('09:30:15')      #--> 142.5625 \n
    obs.hour2deg('09 30 15.000')  #--> 142.5625 \n
    obs.hour2deg([9,30,15.])      #--> 142.5625
    """
    if type(valin)==str:
        if ':' in valin: ra=[float(val)*360./24 for val in valin.split(':')]
        else: ra=[float(val)*360./24 for val in valin.split(' ')]
    else: 
        ra=list(valin); 
        for i in [0,1,2]: ra[i]*=360./24 
    valout=ra[0]+ra[1]/60.+ra[2]/3600.
    return valout

def deg2dms(valin):
    """
    Converts decimal degrees to a list of [ Degrees (int), Minutes (int), Seconds (decimal)]
    
    Parameters
    ----------
    valin: float
        Input value in degrees
    
    Returns
    -------
    DMS : list
        [ Degrees (int), Minutes (int), Seconds (decimal)]
    
    Example
    -------
    # e.g., 123.456789 degrees corresponds to 123 deg, 27 min, 24.4404 sec \n
    obs.deg2dms(123.456789)  #--> [123, 27, 24.440400000002]  
    """
    ddeg=int(valin)
    dmins=int(abs(valin-ddeg)*60)
    dsec=(abs((valin-ddeg)*60)-dmins)*60
    return [int(ddeg),int(dmins),dsec]

def dms2deg(valin):
    """
    Converts DMS input to decimal degrees.
    Input can be either a string delimeted by : or spaces, or a list of [D,M,S] numbers.
    
    Parameters
    ----------
    valin: float
        Input value in DMS. Can be either: \n
        - a string delimeted by : or spaces \n
        - a list of [D,M,S] numbers (floats or ints) \n
    
    Returns
    -------
    valout : float
        Degrees corresponding to the DMS value
    
    Examples
    --------
    # e.g., '-78:12:34.56' corresponds to -77.7904 deg \n
    obs.dms2deg('-78:12:34.56')  #--> -77.79039999999999 \n
    obs.dms2deg('-78 12 34.56')  #--> -77.79039999999999 \n
    obs.dms2deg([-78,12,34.56])  #--> -77.79039999999999
    """
    if type(valin)==str:
        if ':' in valin: ra=[float(val) for val in valin.split(':')]
        else: ra=[float(val) for val in valin.split(' ')]
    else: ra=valin
    valout=ra[0]+ra[1]/60.+ra[2]/3600.
    return valout

def dec2sex(longin,latin,as_string=False,decimal_places=2,str_format=':',RAhours=True,order='radec'):
    """
    Convert from decimal coordinate pairs to sexagesimal format.
    
    Parameters
    ----------
    longin : float
        Longitude coordinate in degrees (e.g., Long on Earth or Right Ascension on sky)
    latin : float
        Latitude coordinate in degrees (e.g., Lat on Earth or Declination on sky ) 
    as_string : bool
        False (default) = return as lists of floats
        True = return as a string
    decimal_places : int 
        Return the result with this number of decimal places 
    str_format : str
        The format in which to return the coords when as_string is True. \n
        Options are ':', ' ', 'DMS', 'dms', 'HMS', 'hms', or user supplied list \n
        of six delimiter strings: \n
           ':' : separate all coordinates by colons.  e.g., ['05:19:00.54','-23:07:31.10'] \n
           ' ' : separate all coordinates by spaces.  e.g., ['05 19 00.54','-23 07 31.10'] \n
           'DMS','dms': separate RA&DEC coords by letters (upper or lower case)  \n
                        e.g., ['05d19m00.54s','-23d07m31.10s'] \n
           'HMS','hms': same as DMS/dms, but RA uses hours instead of degrees    \n
                        e.g., ['05h19m00.54s','-23d07m31.10s'] \n
           user supplied list of 6 delimiters.  \n
            e.g. ['hr','min','sec','deg','min','sec']  -- ['05hr19min00.54sec',..] \n
    RAhours : bool 
        True (default) specifies that the longitude coordinate is in RA hours \n
            (and will be divided by 360/24=15) \n
        False specifies that the longitude coordinate is not in hours
    order : str
        The order of the coordinates: longitude then latitude or vice versa.  \n
        Options: 'lonlat','latlon','radec','decra' \n
        For example,  dec2sex(loncoord,latcoord)                 -- returns [lonstring,latstring]  \n
        or            dec2sex(latcoord,loncoord, order='latlon') -- returns [latstring,lonstring]  \n
    
    Returns
    -------
    LatMinSec : list, or string
        The sexagesimal versions of the input coorinates.
    
    Examples
    --------
    obs.dec2sex(146.50, -14.25)  #--> RAhours=True, so  ([9, 46, 0.0], [-14, 15, 0.0]) \n
    obs.dec2sex(146.50, -14.25, RAhours=False)     #--> ([146, 30, 0.0], [-14, 15, 0.0]) \n
    obs.dec2sex(0.,90., as_string=True) \n
    obs.dec2sex(0., 90., as_string=True, str_format='hms')  #--> ['00h00m00.00s', '90d00m00.00s'] \n
    obs.dec2sex(5.1111, 80.2222, as_string=True, decimal_places=8, str_format=':')  \n
        #-->  ['00:20:26.66400000', '80:13:19.92000000']
    """
    if RAhours is True: RAfactor=24./360.
    else: RAfactor=1.
    if 'lat' in order[:3].lower() or 'dec' in order[:3].lower():
        #If input coords are given in the order [lat,lon] or [dec,ra], switch them here to ensure longin matches with longitude
        tmplong=float(latin); tmplat=float(longin)
        longin=tmplong; latin=tmplat
    longmins,longsecs=divmod(RAfactor*longin*3600,60)
    longdegs,longmins=divmod(longmins,60)
    #latmins,latsecs=divmod(latin*3600,60)
    #latdegs,latmins=divmod(latmins,60)
    latdegs=int(latin)
    if latdegs==0 and latin<0: latdegs=np.NZERO #Case of -0 degrees
    latmins=int(abs(latin-latdegs)*60)
    latsecs=(abs((latin-latdegs)*60)-latmins)*60
    if as_string==True: 
        if latdegs==0 and latin<0: latdegstring='-00'
        else: latdegstring='{0:0>2d}'.format(latdegs)
        #if str_format==':': return ['{0}:{1}:{2:0>{4}.{3}f}'.format(int(longdegs),int(longmins),longsecs,decimal_places,decimal_places+3),'{0}:{1}:{2:0>{4}.{3}f}'.format(int(latdegs),int(latmins),latsecs,decimal_places,decimal_places+3)]
        if str_format==':': delimiters=[':',':','', ':',':','']
        elif str_format==' ': delimiters=[' ',' ','', ' ',' ','']
        elif str_format=='DMS': delimiters=['D','M','S','D','M','S']
        elif str_format=='dms': delimiters=['d','m','s','d','m','s']
        elif str_format=='HMS': delimiters=['H','M','S','D','M','S']
        elif str_format=='hms': delimiters=['h','m','s','d','m','s']
        else: 
            delimiters=str_format
            if len(delimiters)<6: 
                raise Exception('dec2sex(): user-supplied format must have six delimiters.\nSupplied str_format = %s'%str_format)
        lonstring='{0:0>2d}{5}{1:0>2d}{6}{2:0>{4}.{3}f}{7}'.format( int(longdegs),int(longmins),longsecs, decimal_places, decimal_places+3,  delimiters[0],delimiters[1],delimiters[2], )
        latstring='{0}{5}{1:0>2d}{6}{2:0>{4}.{3}f}{7}'.format( latdegstring,int(latmins),latsecs, decimal_places, decimal_places+3, delimiters[3],delimiters[4],delimiters[5], )
        if 'lat' in order[:3].lower() or 'dec' in order[:3].lower(): return [latstring,lonstring]
        else: return [lonstring,latstring]
    else: 
        if 'lat' in order[:3].lower() or 'dec' in order[:3].lower(): 
            return [latdegs,int(latmins),latsecs],[int(longdegs),int(longmins),longsecs]
        else: return [int(longdegs),int(longmins),longsecs],[latdegs,int(latmins),latsecs]

def sex2dec(longin,latin,RAhours=True,order='radec'):
    """
    Convert from sexagesimal coordinate pairs to decimal format.
    
    Parameters
    ----------
    longin : string, or array_like (list, tuple, numpy.array)
        Longitude coordinate in sexagesimal format (e.g., Long on Earth or Right Ascension on sky) 
    latin : string, or array_like (list, tuple, numpy.array)
        Latitude coordinate in sexagesimal format (e.g., Lat on Earth or Declination on sky ) 
    RAhours : bool 
        True (default) = the longitude coordinate will be RA hours (to be multiplied by 360/24=15) \n
        False = the longitude coordinate will NOT be in RA hours
    order : str
        The order of the coordinates: longitude then latitude or vice versa.   \n
        Options: 'lonlat','latlon','radec','decra' \n
        For example,  sex2dec(lonstring,latstring)                 -- returns [loncoord,latcoord] \n
        or            sex2dec(latstring,lonstring, order='latlon') -- returns [latcoord,loncoord] \n
    
    Returns
    -------
    [latout,longout] : list
        The decimal versions of the input coorinates.
    
    Examples
    --------
    obs.sex2dec('23:59:59.9999','-30:15:15.0')    #--> [359.9999995833333, -30.254166666666666] \n
    obs.sex2dec('23h59m59.9999s','-30d15m15.0s')  #--> [359.9999995833333, -30.254166666666666] \n
    obs.sex2dec([23,59,59.9999],[-30,15,15.0])    #--> [359.9999995833333, -30.254166666666666]
    """
    if RAhours is True or 'h' in longin.lower(): RAfactor=360./24.
    else: RAfactor=1.
    
    ### Check if input was given in list or numpy.array format, and convert to string
    lonstring=str(longin); latstring=str(latin)
    if '[' in lonstring or '(' in lonstring: 
        lonstring=str(list(longin)).replace('[','').replace(']','').replace(' ','').replace(',',':');
    if '[' in latstring or '(' in latstring: 
        latstring=str(list(latin)).replace('[','').replace(']','').replace(' ','').replace(',',':');
    
    if 'lat' in order[:3].lower() or 'dec' in order[:3].lower():
        #If input coords are given in the order [lat,lon] or [dec,ra], switch them here to ensure longin matches with longitude
        tmplong=str(latstring); tmplat=str(lonstring)
        lonstring=tmplong; latstring=tmplat
    if ':' in lonstring: 
        ### Format example: '09:45:42.05'
        lo=[float(val)*RAfactor for val in lonstring.split(':')]
    elif 'h' in lonstring.lower(): 
        ### Right ascension, Format example: '09h45m42.05s' or '09H45M42.05S'
        lotmp=lonstring.lower().replace('h',':').replace('m',':').replace('s','')
        lo=[float(val)*RAfactor for val in lotmp.split(':')]
        #import re
        #ra=[float(val)*360./24 for val in re.split('h|m|s',rain.lower())[:3]]
    elif 'd' in lonstring.lower():
        ### Not right ascension, Format example: '135d45m42.05s' or '135D45M42.05S'
        lotmp=lonstring.lower().replace('d',':').replace('m',':').replace('s','')
        lo=[float(val)*RAfactor for val in lotmp.split(':')]
    else: 
        ### Format example: '09 45 42.05'
        lo=[float(val)*RAfactor for val in lonstring.split(' ')]
    longout=lo[0]+lo[1]/60.+lo[2]/3600.
    
    if ':' in latstring: 
        ### Format example: '-14:19:34.98'
        la=[float(val) for val in latstring.split(':')]
    elif 'd' in latstring.lower():
        ### Format example: '-14d19m34.98s' or '-14D19M34.98S'
        latmp=latstring.lower().replace('d',':').replace('m',':').replace('s','')
        la=[float(val) for val in latmp.split(':')]
    else: 
        ### Format example: '-14 19 34.98'
        la=[float(val) for val in latstring.split(' ')]
    if la[0]<0 or '-' in latin.split(' ')[0]: latout=la[0]-la[1]/60.-la[2]/3600.
    else: latout=la[0]+la[1]/60.+la[2]/3600.
    if 'lat' in order[:3].lower() or 'dec' in order[:3].lower(): return [latout,longout]
    else: return [longout,latout]

def angulardistance(coords1,coords2,pythag_approx=False,returncomponents=False):
    """
    Calculate the angular distance between [RA1,DEC1] and [RA2,DEC2] given in decimal (not sexagesimal). 
    
    Parameters
    ----------
    coords1 : float
        [RA,DEC] coordinates in degrees for the first object (decimal, not sexagesimal) 
    coords2 : float
        [RA,DEC] coordinates in degrees for the second object (decimal, not sexagesimal) 
    pythag_approx : bool
        If True, use the pythagorean approximation for the distance (warning! only approximately valid for very small distances, such as for small fits images.  Likely to be deprecated in a future release.)
    returncomponents : bool
        Setting to True also returns the separated RA and Dec distance components
    
    Returns 
    -------
    totalseparation_deg : float
        Separation in degrees
    
    Notes
    -----
    This assumes coords1 and coords2 are in the SAME coordinate system/equinox/etc.!   \n
    Don't mix ICRS and FK5 for example.  \n
    For small separations (<~1deg), can approximate the result with a small correction  \n
    to the pythagorean formula: \n
       delta_RA=(RA1-RA2)*cos(DEC_average) \n
       delta_DEC=(DEC1-DEC2) \n
       total_sep=sqrt(delta_RA^2+delta_DEC^2) \n
    Otherwise must use the full equation for arbitrary angular separation (central angle  \n
    of great-circle distance) -- for RA mapped to [0,2pi], DEC to [-pi/2,+pi/2]): \n
        ang.dist.=arccos( sin(DEC1)sin(DEC2)+cos(DEC1)COS(DEC2)cos(RA1-RA2) )
    
    Examples
    --------
    obs.angulardistance([146.4247680, -14.3262771], [150.4908485, 55.6797891])  #--> 70.08988039585651 
    """
    if pythag_approx is True: 
        #Reasonable approximation for separations (esp. in RA) of less than ~ 1deg.  
        #Fine for normal small FOV fits images.
        dRA=(coords2[0]-coords1[0])*np.cos(np.mean([coords1[1],coords2[1]])*np.pi/180)
        dDEC=(coords2[1]-coords1[1])
        totalseparation_deg=np.sqrt(dRA**2+dDEC**2)
    else: 
        c1rad=np.array(coords1)*np.pi/180; c2rad=np.array(coords2)*np.pi/180;
        totalseparation_deg = np.arccos( np.sin(c1rad[1])*np.sin(c2rad[1])+np.cos(c1rad[1])*np.cos(c2rad[1])*np.cos(c1rad[0]-c2rad[0]) )*180./np.pi
        #Quick calc individual delta components from full formula after locking the 
        #other axis to its max value.  ...not really as useful for large angles, these 
        #individual components were really just meant for calculating pixel widths...
        dDEC=np.arccos( np.sin(c1rad[1])*np.sin(c2rad[1])+np.cos(c1rad[1])*np.cos(c2rad[1])*1. )*180./np.pi * np.sign(coords2[1]-coords1[1]) #here, deltaRA=0
        maxdec=np.array([c1rad[1],c2rad[1]])[np.argmax(np.abs([c1rad[1],c2rad[1]]))]
        dRA=np.arccos( np.sin(maxdec)**2+np.cos(maxdec)**2*np.cos(c2rad[0]-c1rad[0]) )*180./np.pi *np.sign(coords2[0]-coords1[0])
    if returncomponents==True: return totalseparation_deg,dRA,dDEC
    else: return totalseparation_deg



##### ------------ Utility functions for dt.datetime and pyephem.Date objects ------------ #####

def construct_datetime(listin,dtformat='time',timezone=None):
    """
    Convenience utility to construct a dt.datetime() object from a list of components, 
    including an optional timezone.  
    
    Parameters
    ----------
    listin : list (or other array-like)
        List of [H,M,S(.SS)] (for 'time'), [Y,M,D] (for 'date'), or [Y,M,D,H,M,S(.S)] (datetime) \n
        Strings also accepted in format 'H:M:S(.s)', 'Y/M/D', or 'Y/M/D H:M:S(.s)'
    dtformat : str
        'time', 'date', or 'datetime'/'dt'.  Denotes desired format.
    timezone : pytz.timezone, or datetime timezone string 
        Accepts timezones in pytz.timezone format or datetime string format (standard
        Olson database names such as 'US/Eastern', which will then be formed with pytz)
    
    Returns
    -------
    dt_out : datetime object
        dt.time, dt.date, or dt.datetime object
    
    Examples
    --------
    obs.construct_datetime([2021,9,15,16,0,0.0],'dt',timezone='US/Pacific')
    obs.construct_datetime('2021/09/15 16:00:00','dt',timezone='US/Pacific') \n
    #-->  datetime.datetime(2021, 9, 15, 16, 0, tzinfo=<DstTzInfo 'US/Pacific' PDT-1 day, 17:00:00 DST>)
    """
    try: dt.datetime;
    except: import datetime as dt
    if dtformat=='time':
        if type(listin)==str: 
            try: dt_out=dt.datetime.strptime(listin,'%H:%M:%S').time()
            except: dt_out=dt.datetime.strptime(listin,'%H:%M:%S.%f').time()
        elif type(listin[-1])==float: dt_out=dt.datetime.strptime('{}:{}:{}'.format(*listin),'%H:%M:%S.%f').time()
        else: dt_out=dt.time(*listin)
    elif dtformat=='date': 
        if type(listin)==str: dt_out=dt.datetime.strptime(listin,'%Y/%m/%d').date()
        else: dt_out=dt.date(*listin)
    elif dtformat=='datetime' or dtformat=='dt':
        if type(listin)==str: 
            try: dt_out=dt.datetime.strptime(listin,'%Y/%m/%d %H:%M:%S')
            except: dt_out=dt.datetime.strptime(listin,'%Y/%m/%d %H:%M:%S.%f')
        elif type(listin[-1])==float: dt_out=dt.datetime.strptime('{}/{}/{} {}:{}:{}'.format(*listin),'%Y/%m/%d %H:%M:%S.%f')
        else: dt_out=dt.datetime(*listin)
    else: raise Exception("dtformat must be 'date', 'time', 'datetime' (or 'dt')")
    if timezone is not None: 
        try: pytz;
        except: import pytz
        if type(timezone)==str: timezone=pytz.timezone(timezone)
        #dt_out=dt_out.replace(tzinfo=timezone) #No, if local timezone is supplied, attach it to naive time with localize
        dt_out=timezone.localize(dt_out)
    return dt_out

def convert_ephem_datetime(ephem_date_in):
    """
    Converts naive UTC pyephem Date() object into formatted datetime.datetime 
    object with pytz UTC timezone
    
    Parameters
    ----------
    ephem_date_in : ephem.Date 
        ephem.Date() object with no timezone info (naive)
    
    Returns
    -------
    obstime_dt : datetime.datetime 
        datetime.datetime object with UTC pytz timezone
    """
    tmp_dt=ephem_date_in.datetime()
    obstime_dt = dt.datetime(tmp_dt.year, tmp_dt.month, tmp_dt.day, tmp_dt.hour, tmp_dt.minute, tmp_dt.second, tmp_dt.microsecond, tzinfo=pytz.timezone('UTC'))
    return obstime_dt

def convert_ephem_observer_datetime(obsin):
    """
    Converts pyephem Observer().date object into formatted datetime.datetime 
    object with pytz UTC timezone
    Not currently used; better to just use ephem.Date().datetime() 
    
    Parameters
    ----------
    obsin : ephem.Observer() 
    
    Returns
    -------
    obstime_dt : datetime.datetime 
        datetime.datetime object with UTC pytz timezone
    """
    tmp_dt=obsin.date.datetime()
    obstime_dt = dt.datetime(tmp_dt.year, tmp_dt.month, tmp_dt.day, tmp_dt.hour, tmp_dt.minute, tmp_dt.second, tmp_dt.microsecond, tzinfo=pytz.timezone('UTC'))
    return obstime_dt

def calculate_current_utcoffset(local_timezone):
    """
    Calculate the current UTC offset for a timezone string 
    (standard Olson database name, e.g. 'America/Chicago')
    
    Parameters
    ----------
    local_timezone : pytz.timezone, or string
        Local timezone, either as a pytz.timezone object, or a standard
        Olson database name string
    
    Returns
    -------
    utcoffset : float
        Number of hours the timezone is (currently) offset from UTC
    
    Example
    -------
    dt.datetime.utcnow()  #--> datetime.datetime(2021, 9, 9, 3, 20, 7, 791748)
    obs.calculate_current_utcoffset('America/Chicago')  #--> -5.0
    """
    if type(local_timezone) is str: local_timezone=pytz.timezone(local_timezone)
    utcoffset_delta=dt.datetime.now(local_timezone).utcoffset() #utcoffset() returns a dt.timedelta
    return utcoffset_delta.total_seconds()/3600.

def dt_naive_to_dt_aware(datetime_naive,local_timezone):
    """
    Convert a naive dt object (no tzinfo attached) to an aware dt object.
    
    Parameters
    ----------
    datetime_naive : datetime.datetime
        naive datetime object -- e.g. dt.datetime(2021, 4, 15, 23, 59, 59)
    local_timezone : pytz.timezone, or str
        pytz.timezone object, or string of a standard Olson database 
        timezone name (e.g. 'US/Pacific')
    
    Example
    ------- 
    obs.dt_naive_to_dt_aware(dt.datetime.strptime('2021/04/15 23:59:59','%Y/%m/%d %H:%M:%S'),'US/Pacific')
    """
    if type(local_timezone) is str: local_timezone=pytz.timezone(local_timezone)
    #datetime_aware=dt.datetime(datetime_naive.year, datetime_naive.month, datetime_naive.day, datetime_naive.hour, datetime_naive.minute, datetime_naive.second, datetime_naive.microsecond,local_timezone)
    ##Here: just use .localize with the naive time.  localize applies the tzinfo to the naive time, including DST info.
    datetime_aware=local_timezone.localize(datetime_naive)
    return datetime_aware
    
def calculate_dtnaive_utcoffset(datetime_naive,local_timezone):
    """
    Calculate the UTC offset for a specific time and timezone string
    
    Parameters
    ----------
    datetime_naive : datetime.datetime 
        naive dt.datetime format time in the local desired timezone \n
        If your dt.datetime object already has the correct tzinfo (it's "aware"), then you can  
        just use utcofffset=datetime_aware.utcoffset().total_seconds()/3600.
    local_timezone : pytz.timezone, or str
        pytz.timezone object, or string of a standard Olson database 
        timezone name (e.g. 'US/Pacific')
    
    Returns
    -------
    utc_offset : float
        Number of hours the supplied datetime object with timezone is offset from UTC
    """
    #
    """
    if tz_string is not '': 
        try: datetime_aware=pytz.timezone(tz_string).localize(datetime_in) #localize is preferable to replace, for DST etc.
        except: datetime_aware=datetime_in.replace(tzinfo=pytz.timezone(tz_string)) #but localize only works on naive (UTC) times
    else: datetime_aware=datetime_in
    #Could alternatively just manually construct a new dt.datetime with the tzinfo...
    datetime_aware=dt.datetime(dt.datetime.strptime(print(datetime_in),'%Y-%m-%d %H:%M:%S'),tzinfo=pytz.timezone(tz_string))
    dt_aware=dt_naive_to_dt_aware(datetime_in)
    """
    datetime_aware=dt_naive_to_dt_aware(datetime_naive,local_timezone)
    utcoffset_delta=datetime_aware.utcoffset() #utcoffset() returns a dt.timedelta
    return utcoffset_delta.total_seconds()/3600.

def calculate_dt_utcoffset(datetime_aware):
    """
    Calculate the UTC offset for a specific dt object that is timezone-aware 
    (tzinfo already included).
    
    Parameters
    ----------
    datetime_in : datetime.datetime
        dt.datetime format time in the local desired timezone. (If the tzinfo='UTC',
        the offset will naturally be 0)
    
    Returns
    -------
    utc_offset : float
        Number of hours the supplied aware datetime object is offset from UTC
    
    Example
    -------  
    obs.calculate_dt_utcoffset( obs.dt_naive_to_dt_aware(dt.datetime.now(),'US/Eastern') )
    """
    utcoffset_delta=datetime_aware.utcoffset() #utcoffset() returns a dt.timedelta
    return utcoffset_delta.total_seconds()/3600.
    
def localnaive_to_utc(datetime_local_naive,local_timezone):
    """
    Convert a naive (no tzinfo attached) datetime from the local timezone to UTC. 
    
    Parameters
    ----------
    datetime_local_naive : datetime.datetime
        datetime object with values reflecting the time in the local timezone
    local_timezone : pytz.timezone, or str
        pytz.timezone object, or string of a standard Olson database 
        timezone name (e.g. 'US/Eastern')
    
    Returns
    -------
    datetime_utc : datetime.datetime
    
    Example
    -------
    obs.localnaive_to_utc(dt.datetime.now(),'US/Mountain')
    """
    #
    """
    datetime_aware=dt.datetime(datetime_naive.year, datetime_naive.month, datetime_naive.day, datetime_naive.hour, datetime_naive.minute, datetime_naive.second, datetime_naive.microsecond, pytz.timezone(local_timezone))
    return pytz.utc.localize(datetime_local) #--> only valid if datetime_in is already in local time?
    return datetime_aware.replace(tzinfo=pytz.timezone(local_timezone)).astimezone(pytz.utc)
    #Could alternatively just subtract the utcoffset() from the local time...
    """
    datetime_aware=dt_naive_to_dt_aware(datetime_local_naive,local_timezone)
    #utcoffset_hours=datetime_aware.utcoffset().total_seconds()/3600.
    return datetime_aware.astimezone(pytz.utc)

def is_dt_tzaware(dt):
    """
    Simple check for whether a dt.datetime object is timezone-aware or not
    
    Parameters
    ----------
    dt : datetime.datetime
    
    Returns
    -------
    dt_is_tzaware : bool
    
    Examples
    --------
    dt_naive = dt.datetime.strptime('2021/10/31 23:59:59','%Y/%m/%d %H:%M:%S')\n
    obs.is_dt_tzaware(dt_naive) #--> False \n
    obs.is_dt_tzaware(dt_naive.replace(tzinfo=pytz.UTC)) #--> True
    """
    dt_is_tzaware = dt.tzinfo is not None and dt.tzinfo.utcoffset(dt) is not None
    return dt_is_tzaware

def local_to_utc(datetime_local_aware):
    """
    Convert an aware dt object (correct tzinfo already present) to UTC
    
    Parameters
    ----------
    datetime_local_aware : datetime.datetime
        Timezone-aware datetime object
    
    Returns
    -------
    datetime_utc : datetime.datetime
    
    Example
    ------- 
    #using pytz.localize to account for DST \n
    tz = pytz.timezone('US/Mountain') \n
    t = tz.localize(dt.datetime.strptime('2021/04/15 23:59:59','%Y/%m/%d %H:%M:%S')) \n
    obs.local_to_utc(t)
    """
    #datetime_local = dt.datetime(datetime_in.year, datetime_in.month, datetime_in.day,datetime_in.hour, datetime_in.minute, datetime_in.second, datetime_in.microsecond, pytz.timezone(local_timezone))
    #return pytz.utc.localize(datetime_local) #--> only valid if datetime_in is already in local time?
    #Could alternatively just subtract the utcoffset() from the local time...
    #return datetime_local.replace(tzinfo=pytz.timezone(local_timezone)).astimezone(pytz.utc)
    if not is_dt_tzaware(datetime_local_aware):
        raise Exception('  local_to_utc(): input datetime object has no tzinfo!')
    return datetime_local_aware.astimezone(pytz.utc)

def local_to_local(dt_aware_in,timezone_out):
    """
    Converts/localizes a tz-aware datetime to another specified timezone. 
    
    Parameters
    ----------
    dt_aware_in : datetime.datetime
        Timezone-aware dt.datetime object (tzinfo is attached)
    timezone_out : pytz.timezone or str
        The desired output timezone, in either pytz.timezone format or a string from the Olson timezone name database
    
    Returns
    -------
    dt_aware_out : datetime.datetime
        A timezone-aware datetime object in the new timezone
    
    Examples
    --------
    dt_pacific = dt.datetime(2025,12,31,23,59,59,tzinfo=pytz.timezone('US/Pacific')) \n
    dt_eastern = obs.local_to_local(dt_pacific,'US/Eastern') \n
    dt_utc = obs.local_to_local(dt_pacific,pytz.UTC)
    """    
    if not is_dt_tzaware(dt_aware_in):
        raise Exception('  local_to_local(): input datetime object has no tzinfo!')    
    
    if type(timezone_out) is str: timezone_out=pytz.timezone(timezone_out)
    
    dt_aware_out = dt_aware_in.astimezone(timezone_out)
    return dt_aware_out

def utc_to_local(utc_dt,local_timezone,verbose=True):
    """
    Convert a datetime object in UTC (or other timezone) to a dt in a specified timezone
    
    Parameters
    ----------
    utc_dt : datetime.datetime
        dt.datetime object (naive or aware)
    local_timezone : pytz.timezone, or str
        pytz.timezone object, or string of a standard Olson database 
        timezone name (e.g. 'US/Eastern')
    verbose : bool
        True = print warning in case supplied dt_utc is NaN
    
    Returns
    -------
    datetime_local : datetime.datetime
        Timezone-aware dt object in a localized timezone
    
    Example
    ------- 
    obs.utc_to_local(dt.datetime.utcnow(),'US/Eastern')
    """
    if type(local_timezone) is str: local_timezone=pytz.timezone(local_timezone) 
    try:
        if utc_dt.tzinfo is None: utc_dt=pytz.utc.localize(utc_dt)#, is_dst=None)
    except:
        if np.isnan(utc_dt): 
            if verbose==True: print('utc_to_local Warning: input utc_dt was NaN')
            return np.nan
    #local_dt = utc_dt.replace(tzinfo=pytz.utc).astimezone(local_timezone)
    #return local_timezone.normalize(local_dt) # .normalize might be unnecessary
    #utc_naive=dt.datetime(utc_dt.year,utc_dt.month,utc_dt.day,utc_dt.hour,utc_dt.minute, utc_dt.second,utc_dt.microsecond)
    #datetime_local=local_timezone.localize(utc_naive) #--> No, this just applies the tzinfo to the time, not converting.
    datetime_local=utc_dt.astimezone(local_timezone)
    return datetime_local

def dtaware_to_ephem(dt_in):
    """
    Converts timezone-aware dt object into ephem.Date() [UTC time] object
    
    Parameters
    ----------
    dt_in : datetime.datetime
    
    Returns
    -------
    ephem_out : ephem.Date
    """
    return ephem.Date(local_to_utc(dt_in))


def create_local_time_object(time_in, timezone, string_format='%Y/%m/%d %H:%M:%S', return_fmt='dt'):
    """
    Take an input naive time (no timezone explicitly attached yet) in the desired local time, and use the supplied timezone to return a time object in the specified format: datatime, ephem.Date, or string representation.  First converts input time from specified input format to tz-aware dt object, then returns in the specified output format. Timezone-aware dt objects are now also accepted as input, and will be localized from their attached tzinfo to the specified timezone.
    
    Parameters
    ----------
    time_in : ephem.Date(), datetime.datetime, MJD (float), or str formatted as 'YYYY/MM/DD HH:MM:SS' 
        The desired clock time in the local timezone. Can be input as a string formatted for ephem.Date, or as a naive datetime.datetime object (with no timezone attached yet).  If a tz-aware dt.datetime is used, it will use its attached tzinfo and the input timezone to localize.  Users may also input ephem.Date objects, that were previously created from a string corresponding to a local time (and not UTC as pyephem assumes), because this tz-unaware pyephem time will merely be converted to a naive datetime object. 
    timezone : str, or pytz.timezone
        The local timezone to use.  Olson database names or pytz timezones are accepted. 
    string_format : str
        The format of time_in when it's given as a string. Passed to dt.datetime.strptime to parse the input values.  Default is '%Y/%m/%d %H:%M:%S', which corresponds to input as 'YYYY/MM/DD hh:mm:ss'
    return_fmt : str ['dt', 'MJD', or 'str']
        The format in which the start time will be returned. Options are: \n
            'dt' or 'datetime' : dt.datetime \n
            'str' : return a string of the date, using the format code specified in string_format to convert from dt.datetime using dt.datetime.strftime(). \n
            'ephem' : ephem.Date format.  
    
    Returns
    -------
    local_time_formatted : datetime.datetime, ephem.Date, float, or str
        Returns the time in format specified by return_fmt -- datetime, ephem.Date, MJD (float), or string representation
    
    Examples
    --------
    obs.create_local_time_object('2021/10/31 23:59:59','US/Pacific') \n
    obs.create_local_time_object( dt.datetime(2021,10,31,23,59,59, \ \n
        tzinfo=pytz.timezone('US/Eastern')),'US/Pacific') 
    """
    if type(timezone) is str: timezone=pytz.timezone(timezone) 
    
    ### First, convert input time from whatever format it's in to a dt.datetime whose values
    #   represent the time in the local timezone.  (if tz-aware dt was input, localize it now)
    if type(time_in) is dt.datetime : 
        if is_dt_tzaware(time_in):  
            #In case a tz-aware dt is input (e.g., could be from a different timezone), 
            #localize it to get the local time values.  Yes, it will actually be aware 
            #at this point, but later when it's localized again it won't change.
            dt_local = time_in.astimezone(timezone) 
            dt_naive = dt_local.replace(tzinfo=None)
        else: 
            dt_naive = time_in #The input dt.datetime is already naive in this case
    elif type(time_in) is float :
        #Assumes float inputs are MJD.
        #MJD = JD − 2400000.5 , Dublin Julian Day (i.e. in pyephem) = JD − 2415020 
        #--> Dublin Julian Day = MJD - (2415020 - 2400000.5) = MJD - 15019.5
        dt_naive = ephem.Date(time_in-15019.5).datetime()
    elif type(time_in) is str : 
        dt_naive=dt.datetime.strptime(time_in,string_format)
    elif type(time_in) is ephem.Date: 
        #Note: if input is given as ephem.Date format, it is assumed that the value
        # returned by ephem.Date().datetime() will represent the hours/minutes/sec in local time
        dt_naive=time_in.datetime()
    else: raise Exception('  obs.create_local_time_object(): invalid input time_in of type %s.  Must be datetime or string (or ephem.Date)'%(type(time_in)))
    
    dt_aware = dt_naive_to_dt_aware(dt_naive,timezone)
    
    if 'dt' in return_fmt.lower() or 'datetime' in return_fmt.lower(): return dt_aware
    elif 'str' in return_fmt.lower(): return dt_aware.strftime(string_format)
    elif 'ephem' in return_fmt.lower(): return ephem.Date(dt_aware.astimezone(pytz.utc))
    else : 
        raise Exception("create_local_time_object(): invalid return_fmt (%s), must be one of: ['dt', 'str', 'ephem']")
    

def JD(time_in, string_format='%Y/%m/%d %H:%M:%S'):
    """
    Converts a date/time to Julian Date
    
    Parameters
    ----------
    time_in : datetime.datetime, ephem.Date, or string
        The input time.  If given as a string, specify the format using string_format.
    string_format : str
        The format of time_in when it's given as a string.  Passed to dt.datetime.strptime to parse the input values.
    
    Returns
    -------
    jd : float
    
    Notes
    -----
    Equivalent output to pyephem.julian_date \n
    Time since noon Universal Time on Monday, January 1, 4713 BC. \n
    Algorithm from  L. E. Doggett, Ch. 12, "Calendars", p. 606, in Seidelmann 1992  
    """
    if type(time_in)==dt.datetime : datetime_in=time_in
    elif type(time_in)==ephem.Date: datetime_in=time_in.datetime()
    elif type(time_in)==str : datetime_in=dt.datetime.strptime(time_in,string_format)
    else: raise Exception('  obs.JD(): invalid input time_in of type %s.  Must be datetime, ephem.Date, or string'%(type(time_in)))
    
    yr,mo,d,h,m,s = datetime_in.year, datetime_in.month, datetime_in.day, datetime_in.hour, datetime_in.minute, datetime_in.second
    jd = 367.*yr - int((7*(yr+int((mo+9)/12)))/4.) + int((275.*mo)/9.) + d  + 1721013.5 + ((((s/60.)+m)/60+h)/24.)
    return jd

def MJD(time_in, string_format='%Y/%m/%d %H:%M:%S'):
    """
    Converts a date/time to Modified Julian Date
    
    Parameters
    ----------
    time_in : datetime.datetime, ephem.Date, or string
        The input time.  If given as a string, specify the format using string_format.
    string_format : str
        The format of time_in when it's given as a string.  Passed to dt.datetime.strptime to parse the values.
    
    Returns
    -------
    mjd : float
        The Modified Julian Date
    
    Notes
    -----
    Definition of MJD relative to Julian Date:  MJD = JD - 2400000.5
    """
    if type(time_in)==dt.datetime : datetime_in=time_in
    elif type(time_in)==ephem.Date: datetime_in=time_in.datetime()
    elif type(time_in)==str : datetime_in=dt.datetime.strptime(time_in,string_format)
    else: raise Exception('  obs.MJD(): invalid input time_in of type %s.  Must be datetime, ephem.Date, or string'%(type(time_in)))
    return JD(datetime_in)-2400000.5


def GMST(datetime_in_UTC,dms=True):
    """
    For an input UT time, output the GMST.  
    
    Parameters
    ----------
    datetime_in_UTC : datetime.datetime
        dt object in UTC
    dms : bool
        Whether to return in sexagesimal components or decimal format. \n
        True (default) returns GMST in [int degrees, int minutes, float seconds] \n
        False returns GMST in decimal degrees
    
    Returns
    -------
    GMST_24 : list, or float 
        GMST time in 24-hr format
    """
    if datetime_in_UTC.tzinfo==None: raise Exception('datetime_in_local must include pytz.timezone')
    
    #Greenwich Mean Sidereal Time (in s at UT1=0) = 24110.54841 + 8640184.812866*T +
    #  + 0.093104*T^2 - 0.0000062*T^3   ;  T = (JD-2451545.0)/36525.  
    # T is in Julian centuries from 2000 Jan. 1 12h UT1 
    # D is days from 2000 January 1, 12h UT, which is Julian date 2451545.0
    #D=ephem.julian_date(datetime_in_UTC) #Julian Date
    #T=(D-2451545.)/36525. #Centuries since year 2000
    #GMST_raw=24110.54841/3600+8640184.812866/3600*T+0.093104/3600*(T**2)+6.2e-6/3600*(T**3) #off by noticeable number of seconds
    
    ### Using navy.mil definition, GMST in hours -- http://aa.usno.navy.mil/faq/docs/GAST.php
    datetime_prevmidnight=dt.datetime(datetime_in_UTC.year,datetime_in_UTC.month,datetime_in_UTC.day,0,0,0,0, pytz.utc)
    JD_0=ephem.julian_date(datetime_prevmidnight) #Julian date of Previous midnight
    D=(ephem.julian_date(datetime_in_UTC)-2451545.) 
    D_0=(ephem.julian_date(datetime_prevmidnight)-2451545.) 
    H=(D-D_0)*24. #The number of UT hours past the previous UT midnight
    T=D/36525. #Centuries since year 2000. T should be input time, not previous midnight.
    GMST_raw=6.697374558+0.06570982441908*D_0+1.00273790935*H+0.000026*(T**2) 
    ##Approximation accurate to 0.1 second per century:
    #GMST_raw=18.697374558+24.06570982441908*D #Hours --> Gives same as raw term above including H, but that's also off...
    
    GMST_24=wrap_24hr(GMST_raw) #float GMST hour in 24 format
    #GMST_out=wrap_24hr(GMST_24+ut*1.002737909)
    if dms==True: GMST_24=deg2dms(GMST_24) #[int degrees, int minutes, float seconds]
    return GMST_24 #decimal degrees

def LMST(GMST_in_hours,longitude_east_deg):
    """
    Compute LMST for GMST and longitude. \n
    LMST = GMST + (observer's east longitude)
    
    Parameters
    ----------
    GMST_in_hours : list 
        GMST in [int h, int m, float s]
    longitude_east_deg : float 
        East longitude, in degrees
    
    Returns
    -------
    LMST_hours : float
    """
    ##longitude_east_hour=np.array(deg2hour(longitude_east_deg))
    ##LMST_hours=GMST_in_hours+longitude_east_hour #--> minutes and seconds can go over 60...
    
    #longitude_east_hours=deg2hour(longitude_east_deg)
    #LMST_hours=np.array(GMST_in_hours)+np.array(longitude_east_hours)
    #LMST_hours[0]=wrap_24hr(LMST_hours[0])
    #return [int(LMST_hours[0]),int(LMST_hours[1]),LMST_hours[2]]
    
    GMST_in_deg=hour2deg(GMST_in_hours) #Full 360 deg range
    LMST_deg=wrap_360(GMST_in_deg+longitude_east_deg)
    LMST_hours=deg2hour(LMST_deg)
    return LMST_hours

def HA(LMST_hms,RA_in_decimal):
    """
    Compute hour angle from LMST and RA. \n
    Hour Angle = L(M)ST - Right Ascension
    
    Parameters
    ----------
    LMST_hms : list
        LMST in [H,M,S] format
    RA_in_decimal : float
        Right Ascension in decimal format
        
    Returns
    -------
    HA_hms : list
        Hour angle, in [H,M,S] format
    """
    LMST_decimal=LMST_hms[0]/24.+LMST_hms[1]/60.+LMST_hms[2]/3600.
    HA_decimal=LMST_in_decimal-RA_in_decimal
    return deg2hour(HA_decimal)

def LST_from_local(datetime_in_local,longitude_east_deg):
    """
    Compute Local Sidereal Time from timezone-aware datetime. 
    See http://www.stargazing.net/kepler/altaz.html
    
    Parameters
    ----------
    datetime_in_local : datetime.datetime
        Timezone-aware dt.datetime object (pytz.timezone included)
    longitude_east_deg : float 
        East longitude, in degrees
    
    Returns
    -------
    LST_hms : list
        Local Sidereal Time in [H,M,S format]
    """
    UT=datetime_in_local.astimezone(pytz.utc)
    UT_prevmidnight=dt.datetime(datetime_in_local.year,datetime_in_local.month,datetime_in_local.day,0,0,0,0, datetime_in_local.tzinfo).astimezone(pytz.utc)
    H=(JD(UT)-JD(UT_prevmidnight))*24. #The number of hours past the previous midnight
    
    #d = days from J2000, including fraction of a day.  (JD for 1200 hrs UT on Jan 1st 2000 AD is 2451545.0)
    d =JD(UT)-2451545.
    
    LST_deg = 100.46 + 0.985647 * d + longitude_east_deg + 15*H   #In degrees
    LST_hms = deg2hour(wrap_360(LST_deg))
    return LST_hms


def pytz_timezones_from_utc_offset(tz_offset, common_only=True):
    """
    Determine timezone strings corresponding to the given timezone (UTC) offset
    
    Parameters
    ----------
    tz_offset : int, or float 
        Hours of offset from UTC
    common_only : bool
        Whether to only return common zone names (True) or all zone names (False)
    
    Returns
    -------
    results : list
        List of Olson database timezone name strings
    
    Examples
    --------
    obs.pytz_timezones_from_utc_offset(-7)
    obs.pytz_timezones_from_utc_offset(-7.62, common_only=False)
    """
    #pick one of the timezone collections (All possible vs only the common zones)
    timezones = pytz.common_timezones if common_only else pytz.all_timezones

    # convert the float hours offset to a timedelta
    offset_days, offset_seconds = 0, int(tz_offset * 3600)
    if offset_seconds < 0:
        offset_days = -1
        offset_seconds += 24 * 3600
    desired_delta = dt.timedelta(offset_days, offset_seconds)

    # Loop through the timezones and find any with matching offsets
    null_delta = dt.timedelta(0, 0)
    results = []
    for tz_name in timezones:
        tz = pytz.timezone(tz_name)
        non_dst_offset = getattr(tz, '_transition_info', [[null_delta]])[-1]
        if desired_delta == non_dst_offset[0]:
            results.append(tz_name)

    return results

def date_linspace(start_dt,end_dt,n_elements):
    """
    From start & end dt times, generate an array of datetime objects, similar to np.linspace. 
    Useful for, e.g., making an array of times for calculating ephemeris.
    
    Parameters
    ----------
    start_dt : datetime.datetime
    end_dt : datetime.datetime
    n_elements : int
        number of elements in the final array, including start/end values (as in np.linspace usage)
    
    Returns
    -------
    dt_array : list
        Datetime objects in regular increments between start/end point 
    
    Example
    -------
    a=dt.datetime.strptime('2020/09/15 00:00:00','%Y/%m/%d %H:%M:%S')
    b=dt.datetime.strptime('2020/09/15 23:59:59','%Y/%m/%d %H:%M:%S')
    obs.date_linspace(a,b,10)
    """
    delta_dt=(end_dt-start_dt)/(n_elements-1)
    increments=range(0,n_elements)*np.array([delta_dt]*(n_elements))
    return start_dt+increments

def create_obstime_array(timestart,timeend,timezone_string='UTC',output_as_utc=False,n_steps=100):
    """
    Takes start & end times in either dt or string format, to create an array of 
    datetime times. Does this by converting the input date/time string to datetime 
    format, then calculating the increments with obs.date_linspace()
    
    Parameters
    ----------
    timestart : str, or datetime.datetime
        Start time, either in string format as 'YYYY/MM/DD HH:MM:SS' or in naive datetime format
    timeend : str, or datetime.datetime
        End time, either in string format as 'YYYY/MM/DD HH:MM:SS' or in naive datetime format
    timezone_string : str
        Standard timezone string for use with pytz. 'UTC', 'GMT', 'US/Mountain', 'Atlantic/Canary', etc..
    output_as_utc : bool
        Set to True to receive output in UTC
    n_steps : int
        number of elements in the final array, including start/end values (as in np.linspace usage)
    
    Returns
    -------
    times_arr : list
        Datetime objects in regular increments between start/end point 
    
    Example
    -------
    obs.create_obstime_array('2020/09/15 00:00:00','2020/09/15 23:59:59',timezone_string='UTC',n_steps=10)
    """
    #tstart_utc=ephem.Date(timestart)
    #tend_utc=ephem.Date(timeend)
    #tstart_utc=dt.datetime.strptime(str(timestart),'%Y/%m/%d %H:%M:%S')#.time()
    #tend_utc=dt.datetime.strptime(str(timeend),'%Y/%m/%d %H:%M:%S')#.time()
    
    if type(timestart) is dt.datetime: tstart_local=timestart
    else: tstart_local=dt.datetime.strptime(timestart,'%Y/%m/%d %H:%M:%S')
    if type(timeend) is dt.datetime: tend_local=timeend
    else: tend_local=dt.datetime.strptime(timeend,'%Y/%m/%d %H:%M:%S')
    # Apply timezone indicator.  If UTC, just use that as the 'local' timezone.
    #tstart_local=tstart_local.replace(tzinfo=pytz.timezone(timezone_string)) #No, use localize to attach local timezone to naive.
    #tend_local=tend_local.replace(tzinfo=pytz.timezone(timezone_string))  
    tstart_local=pytz.timezone(timezone_string).localize(tstart_local) 
    tend_local=pytz.timezone(timezone_string).localize(tend_local) 
    
    times_arr=date_linspace(tstart_local,tend_local,n_steps)
    
    if output_as_utc is True:
        #Default is to return the array as times in the local timezone.  
        #This step converts if desired output is UTC 
        for i in range(len(times_arr)): 
            times_arr[i]=times_arr[i].astimezone(pytz.utc)
    
    return times_arr



##### ------------ Utility functions for pyephem objects and calculations ------------ #####

def create_ephem_target(namestring, RA, DEC, decimal_format='deg'):
    """
    Convenience function to create an ephem.FixedBody() object for a sky source target
    
    Parameters
    ----------
    namestring : str
        The name/identifier to attach to the object (e.g., 'NGC1275')
    RA : string, or float
        Right Ascension coordinate. If string, must be in sexagesimal H:M:S form (separated by : ). 
        If float, must be in RADIANS (which is what pyephem assumes internally) 
    DEC : string, or float
        Declination coordinate. If string, must be in sexagesimal D:M:S form (separated by : ). 
        If float, must be in RADIANS (which is what pyephem assumes internally) 
    decimal_format : str
        Specifies the format of the input RA & DEC when they are given as floats -- 
        'deg'/'degrees' if the inputs are in degrees, or 'rad'/'radians' if the inputs are in 
        radians (which are the native float input format for pyephem FixedBody type)
    
    Returns
    -------
    ephemobj : ephem.FixedBody
        A pyephem FixedBody() instance
    
    Examples
    --------
    ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751') \n
    # obs.sex2dec('02:41:04.7985','-08:15:20.751')   #--> [40.269994, -8.255764] \n
    ngc1052=obs.create_ephem_target('NGC1052', 40.269994 * np.pi/180, -8.255764 * np.pi/180) \n
    #in radians: np.array([40.269994, -8.255764])*np.pi/180  #--> [ 0.70284399, -0.14409026] \n
    ngc1052=obs.create_ephem_target('NGC1052', 0.70284399, -0.14409026, decimal_format='rad')
    
    Notes
    -----
    This FixedBody object's internal values can be accessed and displayed in a number of ways. By default, calling a pyephem object property returns its numerical value, while using print() on that property returns a more common human-readable format, such as sexagesimal coordinates.\n
    \n
    ngc1052.ra,ngc1052.dec   #-->  (0.7075167104622412, -0.142444748709817) [radians] \n
    print(ngc1052.ra,ngc1052.dec)   #-->  2:42:09.05 -8:09:41.3 \n
    \n
    Information depending on time, such as altitude & azimuth from an Observer, can be computed using the FixedBody.compute(Observer) method. See the pyephem documentation for more details. \n
    https://rhodesmill.org/pyephem/quick.html#bodies\n
    \n
    Input coordinates are assumed to be apparent topocentric position.  \n
    See https://rhodesmill.org/pyephem/radec for more details on pyephem RA & DEC coordinates.
    """
    ### If the input coordinates are floats in degrees, convert to radians for input to pyephem
    if type(RA) not in [str,np.str_] and 'deg' in decimal_format.lower():
        RA=float(RA)*np.pi/180
    if type(DEC) not in [str,np.str_] and 'deg' in decimal_format.lower():
        DEC=float(DEC)*np.pi/180
    
    ephemobj=ephem.FixedBody() #MUST create blank FixedBody() object first and add ._ra & ._dec afterwards!
    ephemobj.name=namestring 
    #ngc1052._ra=ephem.hours(40.269994*np.pi/180); ngc1052._dec=ephem.degrees(-8.255764*np.pi/180); ngc1052.name='NGC1052' 
    ephemobj._ra=ephem.hours(RA); ephemobj._dec=ephem.degrees(DEC); 
    ephemobj.compute() #Initialize RA,DEC with an empty .compute() 
    return ephemobj

class Observer_with_timezone(ephem.Observer):
    """
    A lightly-decorated version of a pyephem Observer object, which now includes a timezone attribute.  
    This timezone can be set to a string name (Olson database) of the timezone corresponding to the
    Observer's location (longitude/latitude).
    """
    def __init__(self):
        super().__init__()
        self.timezone=None
    
    ### Just add the auto-calculate-timezone function here as a class method?
    #def calculate_observer_timezone()...

#class Date_with_timezone(ephem.Date):
#    """
#    A lightly-decorated version of a pyephem.Date object, which now includes a timezone attribute.
#    This timezone can be set to a string name (Olson database) of the timezone corresponding to the
#    desired local time. 
#    
#    NOTE: The timezone does not affect the stored time value -- the numerical value already 
#    contained in the Date object (in Dublin Julian Days) still represents the time in UTC.
#    The timezone is used for localization to a specific timezone, and conversions for plotting 
#    functions.  
#    """
#    def __init__(self):
#        super().__init__()
#        self.timezone=None
#        ### If input time is a tz-aware datetime, convert to UTC for storage with ephem.Date, and
#        #   also also initialize with the specified timezone. [pytz.UTC.zone --> 'UTC']
#        #...
#    
#    ### Add some convenience methods here?  e.g., to localize (return local dt) to the specified timezone?
#    #def localize()...


def create_ephem_observer(namestring, longitude, latitude, elevation, decimal_format='deg', timezone=None):
    """
    Convenience function to create a (lightly modified) ephem.Observer() object for a telescope location.  Atimezone attribute has been added to Observers created with this function.
    
    Parameters
    ----------
    namestring : str
        The name/identifier to attach to the observatory (e.g., 'JCMT', 'GBT', ...)
    longitude : str or float 
        Pyephem longitudes are +E, or (increasing) East coords (pygeodesy may print coords as +West, if the +E are negative).  String coordinates can be delimited by spaces or colons. Note that decimal coordinates are accepted by pyephem objects, but it expects decimal coordinates to be in radians, not degrees.
    latitude : str or float 
        String coordinates can be delimited by spaces or colons. Note that decimal coordinates are accepted by pyephem objects, but it expects decimal coordinates to be in radians, not degrees.
    elevation : int or float
        The elevation of the telescope / observer site, in meters
    decimal_format : str
        Specifies the format of the input longitude & latitude when they are given as floats -- 
        'deg'/'degrees' if the inputs are in degrees, or 'rad'/'radians' if the inputs are in 
        radians (which are the native float input format for pyephem Observer type)
    timezone : str, None
        The name of the timezone (from the Olson database) corresponding to the Observer's location.  The timezone can be calculated automatically, using tzwhere, from the latitude and longitude, though this option can taka a few seconds.
    
    Returns
    -------
    tmpobserver : ephem.Observer
        A pyephem Observer() instance -- lightly modified with Observer_with_timezone() to include a timezone attribute
    
    Examples
    --------
    wht = obs.create_ephem_observer('WHT', '-17 52 53.8', '28 45 37.7', 2344)  \n
    wht = obs.create_ephem_observer('WHT', '-17:52:53.8', '28:45:37.7', 2344)  \n
    wht = obs.create_ephem_observer('WHT', -17.88161, 28.760472, 2344) \n
    wht = obs.create_ephem_observer('WHT', -0.31209297, 0.50196493, 2344, decimal_format='rad') \n
    \n
    ### Attach timezone info by specifying the name manually, and by calculating from coordinates: \n
    wht = obs.create_ephem_observer('WHT', '-17 52 53.8', '28 45 37.7', 2344, timezone='Atlantic/Canary')\n
    wht = obs.create_ephem_observer('WHT', '-17 52 53.8', '28 45 37.7', 2344, timezone='calculate')
    
    Notes
    -----
    This Observer object's internal values can be accessed and displayed in a number of ways. By default, calling a pyephem object property returns its numerical value, while using print() on that property returns a more common human-readable format, such as sexagesimal coordinates.\n
    \n
    wht.lat,wht.lon          #-->  (0.501964934706148, -0.3120929894500905) [radians]\n
    print(wht.lat,wht.lon)   #-->  28:45:37.7 -17:52:53.8 \n
    \n
    See the pyephem documentation for more details. \n
    https://rhodesmill.org/pyephem/quick.html#observers
    """
    ### If the input coordinates are floats in degrees, convert to radians for input to pyephem
    if type(longitude) is not str and 'deg' in decimal_format.lower():
        longitude=float(longitude)*np.pi/180
    if type(latitude) is not str and 'deg' in decimal_format.lower():
        latitude=float(latitude)*np.pi/180
    
    #tmpobserver=ephem.Observer(); #Must be initialized with blank parentheses
    tmpobserver=Observer_with_timezone()
    tmpobserver.name=namestring;
    tmpobserver.lon=longitude; 
    tmpobserver.lat=latitude; 
    tmpobserver.elevation=elevation; 
    if timezone is not None: 
        if True in [i in timezone.lower() for i in ['auto','calc']]:
            tmpobserver.timezone = autocalculate_observer_timezone(tmpobserver)
        else: tmpobserver.timezone = timezone
    return tmpobserver

def autocalculate_observer_timezone(observer):
    """
    Determine the local timezone from a pyephem Observer object, using tzwhere.  
    
    Parameters
    ----------
    observer : ephem.Observer
        The pyephem Observer object which includes the site latitude & longitude
    
    Returns
    -------
    timezone : str
        The name of the local timezone (from the Olson database)
    
    Example
    -------
    wht = obs.create_ephem_observer('WHT', '-17 52 53.8', '28 45 37.7', 2344) \n
    obs.autocalculate_observer_timezone(wht) \n
    # --> 'Atlantic/Canary'
    """
    try: 
        timezone = tzwhere.tzwhere().tzNameAt( observer.lat*180/np.pi, wrap_pm180(observer.lon*180/np.pi) )
    except: 
        timezone = tzwhere.tzwhere(forceTZ=True).tzNameAt( observer.lat*180/np.pi, wrap_pm180(observer.lon*180/np.pi), forceTZ=True)
    return timezone

def tz_from_observer(observer):
    """
    Determine the local timezone from a pyephem Observer object - if the timezone attribute is already set (not None), use that.  If it's set to None, calculate the local timezone with autocalculate_observer_timezone().  
    
    Parameters
    ----------
    observer : ephem.Observer
        The pyephem Observer object which includes the site latitude & longitude
    
    Returns
    -------
    timezone : str
        The name of the local timezone (from the Olson database)
    """
    #Use tzwhere to compute the timezone based on the observer lat/lon (input in degrees)
    try: 
        if observer.timezone is not None: timezone = observer.timezone
        else: timezone = autocalculate_observer_timezone(observer)
    except: 
        #If the observer was instantiated as a standard ephem.Observer instead of an 
        # obsplanning.Observer_with_timezone class, perform the automatic calculation.
        timezone = autocalculate_observer_timezone(observer)
    return timezone

def ephemeris_report(target,observer,obstime):
    """
    For a specified Observer, astronomical target, and observation time, print to screen various ephemeris information about the target.
    
    Parameters
    ----------
    target : ephem.FixedBody(), ephem.Moon(), or ephem.Sun()
        pyephem target source
    observer : ephem.Observer() 
        The observer/telescope/station as a pyephem Observer object
    obstime : ephem.Date(), or string in format 'YYYY/MM/DD HH:MM:SS' 
        A relevant time for the observations, from which to calculate the rise/set times
    
    """
    tmp_obs=observer.copy(); tmp_obs.date=ephem.Date(obstime)
    tmp_tar=target.copy(); tmp_tar.compute(tmp_obs)
    
    # rise/set time, rise/set azimuth
    print('  Target rises at %s with azimuth %.2f deg, sets at %s with azimuth %.2f deg'%(tmp_tar.rise_time, tmp_tar.rise_az/ephem.degree, tmp_tar.set_time, tmp_tar.set_az/ephem.degree))
    
    # Transit time and altitude
    print('  Target transits at %s with altitude %.2f deg'%(tmp_tar.transit_time,tmp_tar.transit_alt/ephem.degree))
    
    #Whether target NEVER comes up during the night of observations?
    print('  Target %s during this night'%['does not rise' if tmp_tar.neverup else 'rises'][0]); 
    
    # whether target is circumpolar
    print('  Target is %scircumpolar'%('' if tmp_tar.circumpolar==True else 'not '))

    #Observer sidereal time / LST for its current local time...
    print('  For local time of %s, sidereal time (LST) is %s'%(tmp_obs.date,tmp_obs.sidereal_time()))


def calculate_twilight_times(obsframe,startdate,verbose=False):
    """
    Calculate the sunset/sunrise and various twilight times for a classical night observing session.  The expected startdate (date/time) is expected to be in the night, between sunset and sunrise, to calculate previous_setting and next_rising  
    
    Parameters
    ----------
    obsframe : ephem.Observer()
        pyephem Observer() frame
    startdate : str, formatted as e.g., 'YYYY/MM/DD 23:00:00'
         The starting date/time (the evening at the beginning observations, for classical night sessions) 
    verbose : bool
        Set to True to print the sunrise/sunset/twilights info to screen
    
    Returns
    -------
    sunsetrise,t_civil,t_naut,t_astro : np.array,np.array,np.array,np.array
        Numpy arrays of the sunset/sunrise, Civil twilights, Nautical twilights, and Astronomical twilights.  All values are in pyephem.Date format.  
     
    Notes
    -----
    Definitions for the various levels of twilight: \n
    Civil: for horizon = -6 degrees, Nautical: -12 deg, Astronomical: -18 deg
    
    Example
    -------
    # Print times for the Crab Nebula observed from WHT on 2025/01/01 at 23:59:00 \n
    wht = ephem.Observer()
    wht.lat='28 45 37.7'; wht.lon='-17 52 53.8'; wht.elevation=2344
    obstime='2025/01/01 23:59:00'
    sunset,twi_civil,twi_naut,twi_astro = obs.calculate_twilight_times(wht,obstime)
    """
    tmpobsframe=obsframe.copy()
    tmpobsframe.date=startdate
    
    #Sunset/Dawn
    sunset=[tmpobsframe.previous_setting(ephem.Sun(), use_center=True), tmpobsframe.next_rising(ephem.Sun(), use_center=True)]
    #Civil Twilight
    tmpobsframe.horizon='-6'
    t_civil=[tmpobsframe.previous_setting(ephem.Sun(), use_center=True), tmpobsframe.next_rising(ephem.Sun(), use_center=True)]
    #Nautical Twilight
    tmpobsframe.horizon='-12'
    t_naut=[tmpobsframe.previous_setting(ephem.Sun(), use_center=True), tmpobsframe.next_rising(ephem.Sun(), use_center=True)]
    #Astronomical Twilight
    tmpobsframe.horizon='-18'
    t_astro=[tmpobsframe.previous_setting(ephem.Sun(), use_center=True), tmpobsframe.next_rising(ephem.Sun(), use_center=True)]
    
    if verbose==True:
        print('  Sunset :   %s\n  Sunrise :  %s'%(sunset[0],sunset[1]))
        print('  Twilights')
        print('  Civil:        previous start at  %s ,  next end at  %s'%(t_civil[0],t_civil[1]))
        print('  Nautical:     previous start at  %s ,  next end at  %s'%(t_naut[0],t_naut[1]))
        print('  Astronomical: previous start at  %s ,  next end at  %s'%(t_astro[0],t_astro[1]))
    
    return np.array(sunset),np.array(t_civil),np.array(t_naut),np.array(t_astro)

def calculate_target_darktime_singlenight(target,observer,night_time,darklevel='astro'):
    """
    Calculate the available dark time for a given target on a given night.  Dark time here is defined from a specified level: between sunset/sunrise, civil, nautical, or astronomical twilight.  Does not consider light from the moon, only the Sun. 
    
    Parameters
    ----------
    target : ephem.FixedBody(), ephem.Moon(), or ephem.Sun()
        pyephem target source
    observer : ephem.Observer() 
        The observer/telescope/station as a pyephem Observer object
    night_time : ephem.Date(), dt.datettime, or string in format 'YYYY/MM/DD HH:MM:SS' 
        A reference time during the night for the observations. A safe bet is midnight of the desired day, since the calculation is for the previous rise and next set relative to that time. ** Important: do not give a time before sunset or after sunrise, or it will be off by a day.
    darklevel : str
        The level of darkness to consider.  Should be one of 'sunrise'/'sunset', or one of the twilight levels 'civil', 'nautical', 'astronomical'
    
    Returns
    -------
    darktime_singlenight : float
        The number of hours of darktime in the specified night
    
    Examples
    --------
    crab = obs.create_ephem_target('Crab Nebula','05:34:31.94','22:00:52.2') \n
    wht = obs.create_ephem_observer('WHT', '-17 52 53.8', '28 45 37.7', 2344) \n
    obs.calculate_target_darktime_singlenight(crab, wht, '2025/10/31 23:59:59') \n
    # --> 8.867340195807628  [hours] \n
    obs.calculate_target_darktime_singlenight(ephem.Jupiter(), wht, '2025/10/31 23:59:59') \n
    # --> 6.645513419411145  [hours]
    """
    #if target is never up on this night, return 0 hours
    if is_target_never_up(target,observer,night_time): return 0.
    
    R,S=calculate_rise_set_times_single(target,observer,night_time,return_fmt='ephem')
    sunset, twi_civil, twi_naut, twi_astro = calculate_twilight_times(observer,night_time)
    if 'sun' in darklevel.lower(): reftimes=sunset
    elif 'civ' in darklevel.lower(): reftimes=twi_civil
    elif 'naut' in darklevel.lower(): reftimes=twi_naut
    elif 'astro' in darklevel.lower(): reftimes=twi_astro
    else: raise Exception("calculate_target_darktime_single: Invalid input for darklevel, Options: ['sunset','civil,'nautical','astronomical']")
    
    #In the case that target is circumpolar/always up, return the full night's darktime
    if is_target_always_up(target,observer,night_time): return reftimes[1]-reftimes[0]
    
    if R <= reftimes[0]:
        if S < reftimes[0]: darktime_singlenight = 0. #Rises and sets before dark, no dark time
        elif S <= reftimes[1]: darktime_singlenight = S-reftimes[0]  #sets before end of dark
        else: darktime_singlenight = reftimes[1]-reftimes[0]  #Up all through the night
    elif R >= reftimes[0]:
        if S<=reftimes[1]: darktime_singlenight = S-R #rises and sets within night
        else: darktime_singlenight = reftimes[1]-R #rises in night, sets after light
    return darktime_singlenight*24.


def calculate_moon_times(obsframe,startdate,outtype='dec',verbose=False):
    """
    Calculates the moonrise/moonset -- takes into account when sunset/sunrise occur
    
    Parameters
    ----------
    obsframe : ephem.Observer()
        pyephem Observer() frame
    startdate : str, formatted as e.g., 'YYYY/MM/DD 23:00:00'
        The starting date/time (the evening at the beginning observations, for classical night sessions) 
    outtype : str ('decimal', 'tuple', or 'dt')
        Desired output format: \n
        'decimal' ephem.Date() format (e.g., 43489.99976096985) \n
        'tuple' (e.g., (2019, 1, 26, 11, 59, 39.34779482893646) ) \n  
        'dt' datetime format  (e.g., [datetime.datetime(2021, 4, 15, 14, 56, 33), ...]  )
    verbose : bool
        Set to True to print the previous moonrise / next moonset info to screen
    
    Returns
    -------
    moonriseset : list 
        The [moonrise,moonset] times in the specified format
    """
    tmpobsframe=obsframe.copy()
    tmpobsframe.date=startdate
    
    ### Sunset/sunrise
    #suntimes=[tmpobsframe.previous_setting(ephem.Sun(),use_center=True), tmpobsframe.next_rising(ephem.Sun(),use_center=True)]
    
    ### Moonrise/Moonset
    moonrise=tmpobsframe.previous_rising(ephem.Moon(),use_center=True)
    if moonrise.tuple()[2]<tmpobsframe.date.tuple()[2]:
        moonrise=tmpobsframe.next_rising(ephem.Moon(),use_center=True)
    
    tmpobsframe2=obsframe.copy()
    tmpobsframe2.date=moonrise
    moonset=tmpobsframe2.next_setting(ephem.Moon(),use_center=True)
    
    if verbose == True:
        print('  Previous moonrise : %s\n  Next moonset :      %s'%(moonrise,moonset))
    
    if 'tup' in outtype.lower(): moonrise=moonrise.tuple(); moonset=moonset.tuple() 
    elif 'dt' in outtype.lower() or 'dat' in outtype.lower(): 
        #datetime
        moonrise=dt.datetime.strptime(str(moonrise),'%Y/%m/%d %H:%M:%S')#.time()
        moonset=dt.datetime.strptime(str(moonset),'%Y/%m/%d %H:%M:%S')#.time()
    #else: ... already in ephem.Date decimal format
    
    return [moonrise,moonset]

def calculate_target_times(target,observer,obstime,outtype='dec'):
    """
    Compute a target's rise and set times, from a specified observer site and for a given obstime
    
    Parameters
    ----------
    target : ephem.FixedBody(), ephem.Moon(), or ephem.Sun()
        pyephem target source
    observer : ephem.Observer() 
        The observer/telescope/station as a pyephem Observer object
    obstime : ephem.Date(), or string in format 'YYYY/MM/DD HH:MM:SS' 
        A relevant time for the observations, from which to calculate the rise/set times
    outtype : str ('decimal', 'tuple', or 'dt')
        Desired output format: \n
        'decimal' ephem.Date() format (e.g., 43489.99976096985) \n
        'tuple' (e.g., (2019, 1, 26, 11, 59, 39.34779482893646) ) \n  
        'dt' datetime format  (e.g., [datetime.datetime(2021, 4, 15, 14, 56, 33),...]  )
    
    Returns
    -------
    target_riseset : list
        [target_rise,target_set] times in the specified format
    """
    tmp_obs=observer.copy(); tmp_obs.date=ephem.Date(obstime)
    tmp_tar=target.copy(); tmp_tar.compute(tmp_obs)
    target_rise=tmp_obs.next_rising(target,use_center=True) 
    target_set=tmp_obs.next_setting(target,use_center=True) 
    if 'tup' in outtype.lower(): target_rise=target_rise.tuple(); target_set=target_set.tuple() 
    elif 'dt' in outtype.lower() or 'dat' in outtype.lower(): 
        #datetime
        target_rise=dt.datetime.strptime(str(target_rise),'%Y/%m/%d %H:%M:%S')#.time()
        target_set=dt.datetime.strptime(str(target_set),'%Y/%m/%d %H:%M:%S')#.time()
    #else: ... already in ephem.Date decimal format
    return [target_rise,target_set]

def compute_target_altaz_single(target,observer,obstime):
    """
    Compute a target's altitude and azimuth at a single obstime, using the supplied observer location
    
    Parameters
    ----------
    target : ephem.FixedBody(), ephem.Moon(), or ephem.Sun()
        pyephem target source
    observer : ephem.Observer() 
        The observer/telescope/station as a pyephem Observer() object
    obstime : ephem.Date(), or string in format 'YYYY/MM/DD HH:MM:SS' 
        The desired time at which to calculate the observer alt/az
    
    Returns
    -------
    target_altaz : np.array
        Array of the target altitude and azimuth, each in degrees
    """
    tmp_obs=observer.copy(); tmp_obs.date=ephem.Date(obstime)
    tmp_tar=target.copy(); tmp_tar.compute(tmp_obs)
    target_altaz=np.array([tmp_tar.alt,tmp_tar.az])*180./np.pi #convert from radians to degrees
    return target_altaz #in degrees

def compute_target_altaz(target,observer,t1,t2,nsteps=1000):
    """
    Compute a range of target altitude and azimuth values, using the supplied observer location
    
    Parameters
    ----------
    target : ephem.FixedBody(), ephem.Moon(), or ephem.Sun()
        pyephem target source
    observer : ephem.Observer() 
        The observer/telescope/station as a pyephem Observer() object
    t1 : ephem.Date() or datetime.datetime()
        The start time for the array calculations
    t2 : ephem.Date() or datetime.datetime()
        The end time for the array calculations
    nsteps : int
        The number of elements in the final array, inclusive of t1 and t2
    
    Returns
    -------
    target_altaz : np.array  (target_altaz.shape = [2,nsteps])
        The array of altitude/azimuth calculations spanning t1 to t2, over nsteps.  
    
    Examples
    --------
    tstart=ephem.Date('2025/01/01 00:00:00'); tend=ephem.Date('2025/01/01 23:59:59'); \n
    ngc3079=obs.create_ephem_target('NGC3079','10:01:57.80','55:40:47.24') \n
    obs.compute_target_altaz(ngc3079,obs.vlbaLA,tstart,tend,nsteps=100) \n
    obs.compute_target_altaz(ephem.Sun(),obs.vlbaLA,tstart,tend,nsteps=100) 
    """
    if type(t1) is str: t1=ephem.Date(t1)
    if type(t2) is str: t2=ephem.Date(t2)
    tmp_obs=observer.copy(); #tmp_obs.date=ephem.Date(obstime)
    tmp_tar=target.copy(); tmp_tar.compute(tmp_obs)
    target_alt=[]; target_az=[]
    t_all=np.linspace(t1,t2,nsteps)
    for t in t_all:
        tmp_obs.date=t
        tmp_tar.compute(tmp_obs)
        target_alt.append(tmp_tar.alt); target_az.append(tmp_tar.az)
    target_altaz=np.array([target_alt,target_az])*180./np.pi #convert from radians to degrees
    return target_altaz #in degrees

def compute_moonphase(obstime,return_fmt='perc'):
    """
    Compute the phase of the moon, as a fraction or percent, at a specific time.
    
    Parameters
    ----------
    obstime : emphem.Date(), or str in 'YYYY/MM/DD hh:mm:ss' format
        date/time at which to calculate the moon phase
    return_fmt : str  ('perc', 'frac', or 'name')
        'perc' (default) to return moon phase as a percent, or 'frac' to return fraction.  'name' returns the colloquial name of the moon phase, such as 'Waxing Gibbous'
    
    Returns
    -------
    moonphase : float
        The fractional or percent phase of the moon
    """
    moon=ephem.Moon(); moon.compute(obstime)
    moonphase=moon.moon_phase
    if 'frac' in return_fmt.lower(): return moonphase
    elif 'per' in return_fmt.lower(): return moonphase*100
    elif 'name' in return_fmt.lower(): 
        #Definitions used for moon phase names: [min. fraction, max. fraction]
        #Full = 100% illuminated, approximated here as above 99%
        #New = 0% illuminated, approximated here as below 1%
        #name_dict={'Full':[0.99,1.], 'New':[0.0,0.01], 'Gibbous':[0.5,0.99], 'Crescent':[0.01,0.5]}
        
        #First, determine the principal lunar phase (new/full/gibbous/crescent)
        #If gibbous or crescent, compute the phase 1 minute later to determine if it's waxing or waning
        if moonphase<0.01: phasename='New'
        elif moonphase>0.99: phasename='Full'
        else:
            if moonphase<0.5: principalname='Cresent'
            else: principalname='Gibbous' #(moonphase>= 0.5)
            #Now compute the phase 1 minute (1/60/24 days) later
            moon.compute(ephem.Date(obstime)+1./60/24)
            phase_1minlater=moon.moon_phase
            print(phase_1minlater,type(phase_1minlater))
            if phase_1minlater-moonphase>0: phasename='Waxing '+principalname
            else: phasename='Waning '+principalname
        return phasename
    else: raise Exception('compute_moonphase(): return_fmt = "%s" not understood. Options are ["perc","frac"]'%(return_fmt))

# Idea for calling 'name' of moon phase, e.g. 'Waxing crescent', etc....  https://stackoverflow.com/a/26727658

def compute_moon_tracks(observer,obsstart,obsend,nsteps=100):
    """
    Compute the altitudes & azimuths (visibility tracks) for the Moon from a 
    given observer(telescope/station)
    
    Parameters
    ----------
    observer : ephem.Observer() 
        The observer/telescope/station as a pyephem Observer() object
    obsstart : ephem.Date() or datetime.datetime()
        The start time for the array calculations
    obsend : ephem.Date() or datetime.datetime()
        The end time for the array calculations
    nsteps : int
        The number of elements in the final array, inclusive of obsstart and obsend
    
    Returns
    -------
    moon_altaz : np.array  (moon_altaz.shape = [2,nsteps])
        The array of altitude/azimuth calculations spanning obsstart to obsend, over nsteps.  
    """
    
    """
    ### Using astropy:
    from astropy.coordinates import get_moon, SkyCoord, EarthLocation, AltAz
    wht = EarthLocation(lat=28.760472*u.deg, lon=-17.881611*u.deg, height=2344*u.m)
    delta_midnight_1day = np.linspace(-12, 12, 1000)*u.hour
    times_1day = midnight + delta_midnight_1day
    frame_1day = AltAz(obstime=times_1day, location=wht)
    moon_1day = get_moon(times_1day)
    moonaltazs_1day = moon_1day.transform_to(frame_1day)
    """
    
    return compute_target_altaz(ephem.Moon(),observer,obsstart,obsend,nsteps=nsteps)

def compute_sun_tracks(observer,obsstart,obsend,nsteps=100):
    """
    Compute the altitudes & azimuths (visibility tracks) for the Sun from a 
    given observer(telescope/station)
    
    Parameters
    ----------
    observer : ephem.Observer() 
        The observer/telescope/station as a pyephem Observer() object
    obsstart : ephem.Date() or datetime.datetime()
        The start time for the array calculations
    obsend : ephem.Date() or datetime.datetime()
        The end time for the array calculations
    nsteps : int
        The number of elements in the final array, inclusive of obsstart and obsend
    
    Returns
    -------
    sun_altaz : np.array  (sun_altaz.shape = [2,nsteps])
        The array of altitude/azimuth calculations spanning obsstart to obsend, over nsteps.  
    """
    
    """
    ### Using astropy:
    from astropy.coordinates import get_sun, SkyCoord, EarthLocation, AltAz
    wht = EarthLocation(lat=28.760472*u.deg, lon=-17.881611*u.deg, height=2344*u.m)
    delta_midnight_1day = np.linspace(-12, 12, 1000)*u.hour
    times_1day = midnight + delta_midnight_1day
    frame_1day = AltAz(obstime=times_1day, location=wht)
    sunaltazs_1day = get_sun(times_1day).transform_to(frame_1day)
    """
    
    return compute_target_altaz(ephem.Sun(),observer,obsstart,obsend,nsteps=nsteps)

def compute_sidereal_time(observer,t1,as_type='datetime'):
    """
    Compute the sidereal time corresponding to an input time t1, at a specified
    observer location.
    
    Parameters
    ----------
    observer : ephem.Observer() 
        The observer/telescope/station as a pyephem Observer() object
    t1 : ephem.Date() or datetime.datetime()
        The input time
    as_type : str -- 'datetime', 'dms', 'angle_rad', 'angle_deg', or 'string'
        The desired output format of the computed LST
    
    Returns
    -------
    LST_out : datetime.datetime, [D,M,S], float, or str
    """
    tmpobs=observer.copy()
    t1=ephem.Date(t1); #Ensure time is in ephem format, if input was datetime format
    if as_type=='datetime':
        #prev_midnight=ephem.Date(observer.date.datetime().replace(hour=0,minute=0,second=0,microsecond=0))
        #prev_midnight_dt=dt.datetime.strptime(str(prev_midnight),'%Y/%m/%d %H:%M:%S').date()
        startdate=t1.tuple()[:3]
        tmpobsstart=observer.copy()
        tmpobsstart.date=startdate; 
        starttime_lst=tmpobsstart.sidereal_time()
    tmpobs.date=t1
    LST_rad=tmpobs.sidereal_time() #pyephem sidereal_time is ephem.Angle format (radians)
    #lst_dt=dt.datetime.strptime(str(LST_rad),'%H:%M:%S.%f')
    if as_type=='datetime': 
        #LST_dt=dt.datetime.strptime(str(LST_rad),'%H:%M:%S.%f').replace(year=startdate[0], month=startdate[1], day=startdate[2]+[1 if LST_rad<starttime_lst else 0][0])
        LST_dt=dt.datetime.strptime(str(LST_rad),'%H:%M:%S.%f').replace(year=startdate[0], month=startdate[1], day=startdate[2])
        if starttime_lst<LST_rad: LST_dt+=dt.timedelta(days = -1.)
        LST_out=LST_dt
    elif as_type=='string': LST_out=str(LST_rad)
    elif as_type=='dms': LST_out=deg2hour(LST_rad*180./np.pi)
    elif as_type=='angle_rad': LST_out=LST_rad
    elif as_type=='angle_deg': LST_out=LST_rad*180./np.pi
    else: raise Exception("as_type must be in ['datetime', 'dms', 'angle_rad', 'angle_deg', 'string']")
    return LST_out #in degrees

def compute_sidereal_times(observer,t1,t2,nsteps=1000,as_type='datetime'):
    """
    Compute an array of sidereal times between times t1 and t2 over nsteps, 
    at a specified observer location.
    
    Parameters
    ----------
    observer : ephem.Observer() 
        The observer/telescope/station as a pyephem Observer() object
    t1 : ephem.Date() or datetime.datetime()
        The start time for the array calculations
    t2 : ephem.Date() or datetime.datetime()
        The end time for the array calculations
    nsteps : int
        The number of elements/steps in the array for calculating times
    as_type : str -- 'datetime', 'dms', 'angle_rad', 'angle_deg', or 'string'
        The desired output format of the computed LST values
    
    Returns
    -------
    LST_arr : list (with elements of type datetime.datetime, [D,M,S], float, or str)
    """
    tmpobs=observer.copy()
    t1=ephem.Date(t1); t2=ephem.Date(t2) #Ensure times are in ephem format, if input was datetime format
    LST_arr=[]
    t_all=np.linspace(t1,t2,nsteps) #local times
    if as_type=='datetime':
        #prev_midnight=ephem.Date(observer.date.datetime().replace(hour=0,minute=0,second=0,microsecond=0))
        #prev_midnight_dt=dt.datetime.strptime(str(prev_midnight),'%Y/%m/%d %H:%M:%S').date()
        startdate=t1.tuple()[:3]
        tmpobs.date=t1; 
        starttime_lst=tmpobs.sidereal_time()
    for t in t_all:
        tmpobs.date=t
        LST_rad=tmpobs.sidereal_time() #pyephem sidereal_time is ephem.Angle format (radians)
        #lst_dt=dt.datetime.strptime(str(LST_rad),'%H:%M:%S.%f')
        if as_type=='datetime': 
            LST_dt=dt.datetime.strptime(str(LST_rad),'%H:%M:%S.%f').replace(year=startdate[0], month=startdate[1], day=startdate[2]+[1 if LST_rad<starttime_lst else 0][0])
            LST_arr.append(LST_dt)
        elif as_type=='string': LST_arr.append(str(LST_rad))
        elif as_type=='dms': LST_arr.append(deg2hour(LST_rad*180./np.pi))
        elif as_type=='angle_rad': LST_arr.append(LST_rad)
        elif as_type=='angle_deg': LST_arr.append(LST_rad*180./np.pi)
        else: raise Exception("as_type must be in ['datetime', 'dms', 'angle_rad', 'angle_deg', 'string']")
    return LST_arr #in degrees


def calculate_transit_time_single(target, observer, approximate_time, mode='nearest', return_fmt='str'):
    """
    Calculates the transit time of a sky target for the supplied observer(telescope/station), for a single specified approximate date+time.  Use keyword "mode" to specify the transit previous to the input time, the next transit, or the nearest (default, whether it's before or after the input time)
    
    Parameters
    ----------
    target : ephem.FixedBody(), ephem.Sun(), or ephem.Moon()
        The source of interest on the sky
    observer : ephem.Observer() 
        Telescope/observer location from which the observation calculation will be made
    approximate_time : ephem.Date(), or str in 'YYYY/MM/DD HH:MM:SS.s' format
        approximate date/time (UTC, not local!) for which to compute the nearest transit time
    mode : str
        Specifies which transit time to return. Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    return_fmt : str ('str', 'ephem', or 'dt')
        The desired the output format.  'str' for string output, 'ephem' for ephem.Date format, 
        'dt' for datetime.
    
    Returns
    -------
    t_transit : str, ephem.Date, or datetime.datetime
    
    Example
    -------
    ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751') \n
    obs.calculate_transit_time_single(ngc1052,obs.vlbaBR,'2020/09/15 12:00:00',return_fmt='str') \n
    #--> '2020/9/15 11:01:22'
    """
    tmp_ant=observer.copy(); tmp_ant.date=ephem.Date(approximate_time)
    tmp_tar=target.copy(); tmp_tar.compute(tmp_ant)
    if tmp_tar.neverup: 
        try: print('Warning: Target %s is never up at %s for approximate_time = %s'%(target.name,observer.name,approximate_time))
        except: print('Warning: Target is never up at supplied observer at approximate_time = %s'%(approximate_time))
        return np.nan
    else: 
        transits=np.array([tmp_ant.previous_transit(tmp_tar),tmp_ant.next_transit(tmp_tar)])
        if 'prev' in mode.lower() or 'bef' in mode.lower(): t_transit=ephem.Date(transits[0])
        elif 'next' in mode.lower() or 'aft' in mode.lower(): t_transit=ephem.Date(transits[1])
        elif 'near' in mode.lower():
            t_transit=ephem.Date( transits[ np.argmin(np.abs(transits-ephem.Date(approximate_time))) ] )
        if 'eph' in return_fmt.lower(): 
            #return tmp_tar.transit_time #--> deprecated usage
            return t_transit
        elif True in [i in return_fmt.lower() for i in ['dt','datetime']]: 
            #return tmp_tar.transit_time.datetime() #--> deprecated usage
            return t_transit.datetime()
        else: 
            #return str(tmp_tar.transit_time) #--> deprecated usage
            return str(t_transit)

def calculate_rise_set_times_single(target, observer, time_in, mode='nearest', return_fmt='str', verbose=True):
    """
    Calculate a target source's rise/set times near a specified observation time. Depending on the specified 'mode', determine the previous rise and previous set, the next rise next set, or the rise/set corresponding to the nearest transit (default).
    
    Parameters
    ----------
    target = ephem.FixedBody(), ephem.Sun(), or ephem.Moon() 
        The sky source of interest
    observer : ephem.Observer() 
        Telescope/observer location from which the observation calculation will be made
    time_in : ephem.Date(), or str in 'YYYY/MM/DD HH:MM:SS.s' format
        approximate date/time (UTC, not local!) for which to compute the nearest transit time
    mode : str
        Specifies how to determine the rise/set times. Options are \n
        'previous'/'before' = calculate the rise/set times before the input time \n
        'next'/'after' = calculate the rise/set times after the input time \n
        'nearest'/'transit' = calculate the nearest transit to the input time, and the corresponding rise/set times
    return_fmt : str ('str', 'ephem', or 'dt')
        The desired the output format.  'str' for string output, 'ephem' for ephem.Date format, 
        'dt' for datetime.
    verbose : bool
        Set to True to print extra warnings to terminal
    
    Returns
    -------
    risesettimes : list
        [Previous rise time, Next set time]
    
    Examples
    --------
    ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751') \n
    obs.calculate_rise_set_times_single(ngc1052,obs.vlbaLA,'2025/09/15 12:00:00',return_fmt='str') \n
    #--> ['2025/9/15 04:30:38', '2025/9/15 15:47:02']
    """
    tmp_ant=observer.copy(); tmp_ant.date=ephem.Date(time_in)
    tmp_tar=target.copy(); tmp_tar.compute(tmp_ant)
    if tmp_tar.neverup: 
        if verbose==True: 
            try: print('Warning: Target %s is never up at %s for input time = %s'%(target.name,observer.name,time_in))
            except: print('Warning: Target is never up at supplied observer location at input time = %s'%(time_in))
        return [np.nan,np.nan]
    elif tmp_tar.circumpolar:
        if verbose==True:
            try: print('Warning: Target %s is always above horizon at %s for input time = %s'%(target.name,observer.name,time_in))
            except: print('Warning: Target is always above horizon at supplied observer location at input time = %s'%(time_in))
        return [np.nan,np.nan]
    else: 
        if 'prev' in mode.lower() or 'bef' in mode.lower():
            set_time=tmp_ant.previous_setting(tmp_tar)
            tmp_ant.date=set_time; tmp_tar.compute(tmp_ant)
            rise_time=tmp_ant.previous_rising(tmp_tar)
        elif 'next' in mode.lower() or 'aft' in mode.lower():
            rise_time=tmp_ant.next_rising(tmp_tar)
            tmp_ant.date=rise_time; tmp_tar.compute(tmp_ant)
            set_time=tmp_ant.next_setting(tmp_tar)
        elif 'near' in mode.lower() or 'trans' in mode.lower():
            nearest_transit = calculate_transit_time_single(target, observer, time_in, return_fmt='ephem')
            tmp_ant.date=nearest_transit; tmp_tar.compute(tmp_ant)
            rise_time=tmp_ant.previous_rising(tmp_tar)
            set_time=tmp_ant.next_setting(tmp_tar)
        else:
            raise Exception("calculate_rise_set_times_single(): invalid input for 'mode'. Use one of ['nearest','previous','next']")
        
        if 'eph' in return_fmt.lower(): return np.array([rise_time,set_time])
        elif True in [i in return_fmt.lower() for i in ['dt','datetime']]: 
            return np.array([rise_time.datetime(),set_time.datetime()])
        else: return [str(rise_time),str(set_time)]

def is_target_always_up(target,observer,approximate_time):
    """
    Checks whether a target is always above the observing horizon (circumpolar) for a specified observer on a specified date.
    
    Parameters
    ----------
    target = ephem.FixedBody(), ephem.Sun(), or ephem.Moon() 
        The sky source of interest
    observer : ephem.Observer() 
        Telescope/observer location from which the observation calculation will be made
    approximate_time : ephem.Date(), or str in 'YYYY/MM/DD HH:MM:SS.s' format
        approximate date/time (UTC, not local!) of the calculation
    
    Returns
    -------
    alwaysup : bool
        True if target is always up on that day, False if not.
    """
    tmp_ant=observer.copy(); tmp_ant.date=ephem.Date(approximate_time)
    tmp_tar=target.copy(); tmp_tar.compute(tmp_ant)
    return tmp_tar.circumpolar

def is_target_never_up(target,observer,approximate_time):
    """
    Checks whether a target is always below the observing horizon for a specified observer on a specified date.
    
    Parameters
    ----------
    target = ephem.FixedBody(), ephem.Sun(), or ephem.Moon() 
        The sky source of interest
    observer : ephem.Observer() 
        Telescope/observer location from which the observation calculation will be made
    approximate_time : ephem.Date(), or str in 'YYYY/MM/DD HH:MM:SS.s' format
        approximate date/time (UTC, not local!) of the calculation
    
    Returns
    -------
    alwaysup : bool
        True if target is always up on that day, False if not.
    """
    tmp_ant=observer.copy(); tmp_ant.date=ephem.Date(approximate_time)
    tmp_tar=target.copy(); tmp_tar.compute(tmp_ant)
    return tmp_tar.circumpolar

def calculate_antenna_visibility_limits(target,station,referenceephemtime,plusminusdays=1., elevation_limit_deg=15., interpsteps=100, alwaysup_fmt='nan', timeformat='ephem', LST_PT=False, verbose=False):
    """
    Calculate the previous and next time within e.g. 24hrs when a target reaches 
    a specified elevation for the given antenna \n
    e.g., 10deg or 15deg for VLBA stations...  \n
    (sometimes 'always up' targets might still go below 10deg)
    
    Parameters
    ----------
    target = ephem.FixedBody(), ephem.Sun(), or ephem.Moon() 
        The sky source of interest
    station : ephem.Observer() 
        Telescope/station/observer location from which the observation calculation will be made
    referenceephemtime : ephem.Date() 
        Reference observation time.  Can be any time, but expected use is ~ a (mean) transit time. 
    plusminusdays : float
        The number of days to span before/after the reference time
    elevation_limit_deg : float
        The minimum telescope elevation limit to be considered visible/up
    inerpsteps : int
        The number of steps over which to calculate elevations between the reference time and plusminusdays
    alwaysup_fmt : str ('nan','limits')
        The option for handling output when a target is 'always up' with respect to a specified elevation limit \n
        'nan'= return np.nan for the rise/set times \n
        'limits' = return the limits -- the beginning of the day (-plusminusdays) for previous rise/set, and the end of the day (+plusminusdays) for next rise/set
    timeformat : str  ('ephem', 'mjd', 'dt', or dt strptime format)
        String specifying the format for the returned times \n
        'ephem': return times as pyephem.Date objects (default) [NOTE: these floats are Dublin JD, not MJD] \n
        'mjd': return times as float Modified Julian Dates \n
        'dt': return times as datetime.datetime objects \n
        String with strptime % formatting: return string dates (e.g., '%Y/%m/%d %H:%M:%S')
    LST_PT: bool 
        True to return times in LST at Pietown instead of in the default UTC
    
    Returns
    -------
    previous_setrise,next_setrise : list, list 
        [previous set, previous rise], [next set, next rise]
    
    Example
    -------
    ngc3147=obs.create_ephem_target('NGC3147','10:16:53.65','73:24:02.7') \n
    prev_setrise,next_setrise = obs.calculate_antenna_visibility_limits(ngc3147,obs.vlbaSC, ephem.Date('2021/08/01 20:36:31'), elevation_limit_deg=15., verbose=False) 
    """
    try: interpolate
    except: from scipy import interpolate
    
    ### *** NOTE: Developed this way partly in preparation to also return elevation
    #   tracks, but to only return the rise/set times, it's probably far easier to
    #   create temporary station object copies, modify their horizon upwards, and use
    #   the pyephem builtins for e.g.  next_rising
    # tmpstation=station.copy()
    # tmpstation.date=referenceephemtime
    # tmpstation.horizon=str(elevation_limit_deg)
    # tmpstation.previous_setting(target, use_center=True)  
    #...
    
    ## Create arrays of test times and target altitudes
    tmpx_prev=np.linspace(referenceephemtime-plusminusdays,referenceephemtime,interpsteps)
    target_alts_prev=compute_target_altaz(target,station,referenceephemtime-plusminusdays, referenceephemtime, nsteps=interpsteps)[0]
    ## Find the roots (0-crossings) of the test altitude array minus the target value
    alt_matches_prev=interpolate.UnivariateSpline(tmpx_prev,target_alts_prev-elevation_limit_deg,s=0).roots()
    rising_indmask_prev=np.diff(target_alts_prev)>=0; setting_indmask_prev=np.diff(target_alts_prev)<0;
    try:
        alt_matches_setting_prev=setting_indmask_prev[[np.where(i<tmpx_prev)[0][0] for i in alt_matches_prev]]
        prev_set_alt_match=ephem.Date(alt_matches_prev[alt_matches_setting_prev][-1])
        if LST_PT==True:
            prev_set_dt=compute_sidereal_time(vlbaPT,prev_set_alt_match) 
            prev_set_alt_match=ephem.Date(prev_set_dt) #Back to ephem.Date format
    except: 
        if 'nan' in alwaysup_fmt.lower(): prev_set_alt_match=np.nan
        else: prev_set_alt_match=ephem.Date(tmpx_prev[0])
    try:
        alt_matches_rising_prev=rising_indmask_prev[[np.where(i<tmpx_prev)[0][0] for i in alt_matches_prev]]
        prev_rise_alt_match=ephem.Date(alt_matches_prev[alt_matches_rising_prev][-1])
        if LST_PT==True:
            prev_rise_dt=compute_sidereal_time(vlbaPT,prev_rise_alt_match)
            prev_rise_alt_match=ephem.Date(prev_rise_dt) #Back to ephem.Date format
    except:
        if 'nan' in alwaysup_fmt.lower(): prev_rise_alt_match=np.nan
        else: prev_rise_alt_match=ephem.Date(tmpx_prev[0])
    #testalt=interpolate.splev(prev_rise_alt_match,interpolate.splrep(tmpx_prev,target_alts_prev,k=3)) 
    
    ## Create arrays of test times and target altitudes
    tmpx_next=np.linspace(referenceephemtime,referenceephemtime+plusminusdays,interpsteps)
    target_alts_next=compute_target_altaz(target,station,referenceephemtime, referenceephemtime+plusminusdays, nsteps=interpsteps)[0]
    ## Find the roots (0-crossings) of the test altitude array minus the target value
    alt_matches_next=interpolate.UnivariateSpline(tmpx_next,target_alts_next-elevation_limit_deg,s=0).roots()
    rising_indmask_next=np.diff(target_alts_next)>=0; 
    setting_indmask_next=np.diff(target_alts_next)<0;
    try:
        alt_matches_setting_next=setting_indmask_next[[np.where(i<tmpx_next)[0][0] for i in alt_matches_next]]
        next_set_alt_match=ephem.Date(alt_matches_next[alt_matches_setting_next][0])
        if LST_PT==True:
            next_set_dt=compute_sidereal_time(vlbaPT,next_set_alt_match)
            next_set_alt_match=ephem.Date(next_set_dt) #Back to ephem.Date format
    except:
        if 'nan' in alwaysup_fmt.lower(): next_set_alt_match=np.nan
        else: next_set_alt_match=ephem.Date(tmpx_next[0])
    try: 
        alt_matches_rising_next=rising_indmask_next[[np.where(i<tmpx_next)[0][0] for i in alt_matches_next]]
        next_rise_alt_match=ephem.Date(alt_matches_next[alt_matches_rising_next][0])
        if LST_PT==True:
            next_rise_dt=compute_sidereal_time(vlbaPT,next_rise_alt_match)
            next_rise_alt_match=ephem.Date(next_rise_dt) #Back to ephem.Date format
    except: 
        if 'nan' in alwaysup_fmt.lower(): next_rise_alt_match=np.nan
        else: next_rise_alt_match=ephem.Date(tmpx_next[0])
    #testalt_next=interpolate.splev(next_set_alt_match,interpolate.splrep(tmpx_next,target_alts_next,k=3)) 
    
    ### Convert times from MJD to specified output format (actually, probably better to do it above)
    if np.sum(np.isnan([prev_set_alt_match,prev_rise_alt_match, next_set_alt_match,next_rise_alt_match]))==0:
        if timeformat.lower()=='dt' or timeformat.lower()=='datetime':
            prev_set_alt_match=prev_set_alt_match.datetime()
            prev_rise_alt_match=prev_rise_alt_match.datetime()
            next_set_alt_match=next_set_alt_match.datetime()
            next_rise_alt_match=next_rise_alt_match.datetime()
        elif '%' in timeformat: 
            prev_set_alt_match=prev_set_alt_match.datetime().strftime(timeformat)
            prev_rise_alt_match=prev_rise_alt_match.datetime().strftime(timeformat)
            next_set_alt_match=next_set_alt_match.datetime().strftime(timeformat)
            next_rise_alt_match=next_rise_alt_match.datetime().strftime(timeformat)
        elif 'mjd' in timeformat.lower():
            prev_set_alt_match=Time(prev_set_alt_match.datetime()).mjd
            prev_rise_alt_match=Time(prev_rise_alt_match.datetime()).mjd
            next_set_alt_match=Time(next_set_alt_match.datetime()).mjd
            next_rise_alt_match=Time(next_rise_alt_match.datetime()).mjd
        else: pass #Any other arguments, leave as ephem.Date format
    
    return [prev_set_alt_match,prev_rise_alt_match],[next_set_alt_match,next_rise_alt_match]

def calculate_vlbi_mean_transit_time(target,observer_list,approximate_time,weights=None, force_after_observer_i=None, wrap_date=False, mode='nearest'):
    """
    Calculates the mean transit time of a single target for the supplied antennae, for a specified approximate date+time
    
    Parameters
    ----------
    target : ephem.FixedBody(), ephem.Sun(), or ephem.Moon() 
        Sky target of interest
    observer_list : list (elements of ephem.Observer() type) 
        A list of the ephem.Observer() objects for all the stations to calculate for
    approximate_time : ephem.Date(), or str, formatted as 'YYYY/MM/DD HH:MM:SS.s'
        String UTC date/time, to use for calculating the transit times with ephem
    weights : np.array, list, or other array-like
        Optional (float) weights to give to certain stations.  (e.g. to prioritize certain baselines)   Larger values give higher importance. 
    force_after_observer_i : None or int  
        If set to an integer, force the computation of the next transit to be after 
        the transit time at this observer/station index in the supplied observer_list.  
        i.e., useful for forcing all times to be after the St. Croix station transit 
        for VLBA calculations.
    wrap_date : bool 
        Set to True to wrap the mean transit time to restrict it to the same day 
        as the supplied "approximate_time".  Useful when a calculation requires a specific date. 
        (e.g., if mean transit is 03:24 the following day, forces it to return 03:24 on the same day)
    mode : str
        Specifies which transit time to return. Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    
    Returns
    -------
    mean_transit_time : str
        The mean of the transit times of the target at the various stations.  The output is produced by calling the ephem.Date object as a string
    
    Notes
    -----
    Weighted mean = Sum(w_i*x_i)/Sum(w_i)
    
    Examples
    --------
    ngc3079=obs.create_ephem_target('NGC3079','10:01:57.80','55:40:47.24') \n
    obs.calculate_vlbi_mean_transit_time(ngc3079,[obs.vlbaBR,obs.vlbaMK,obs.vlbaSC],'2021/09/15 12:00:00') \n
    #To upweight MK and SC over BR, given them higher values in this list \n
    obs.calculate_vlbi_mean_transit_time(ngc3079,[obs.vlbaBR,obs.vlbaMK,obs.vlbaSC],'2021/09/15 12:00:00',weights=[1,5,2]) 
    """
    if weights is not None:
        if len(observer_list)!=len(weights): 
            raise Exception('weights must be None or an array the same length as observer_list array')
    
    transit_times=np.zeros(len(observer_list))
    for i in range(len(observer_list)):
        tmp_ant=observer_list[i].copy()
        tmp_ant.date=ephem.Date(approximate_time)
        tmp_tar=target.copy()
        tmp_tar.compute(tmp_ant)
        if tmp_tar.neverup: transit_times[i]=np.nan
        else: 
            #transit_times[i]=tmp_tar.transit_time #deprecated usage
            transit_times[i]=calculate_transit_time_single(tmp_tar,tmp_ant,approximate_time, mode=mode, return_fmt='ephem')
    
    ## If a specific reference observer is given, recalculate the next transit times after that one.
    if isinstance(force_after_observer_i, (int, np.integer)):
        force_earliest_time=transit_times[force_after_observer_i]
        for i in range(len(observer_list)):
            if i == force_after_observer_i: continue
            if transit_times[i] < force_earliest_time:
                tmp_ant=observer_list[i].copy()
                tmp_ant.date=ephem.Date(approximate_time)
                tmp_tar=target.copy()
                transit_times[i]=tmp_ant.next_transit(tmp_tar, start=force_earliest_time)
    
    nanmask=np.isnan(transit_times)
    if True in nanmask: 
        neveruplist=[]
        for i in range(len(observer_list)):
            try: neveruplist.append(observer_list[i].name)
            except: neveruplist.append('Observer_i='+str(i))
        #raise Warning('Target never up, for input date, at station(s):\n%s'%(neveruplist))
        print('\nWarning: Target never up, for input date %s, at station(s):\n%s\n'%(approximate_time,neveruplist))
    
    if weights is not None: 
        meantransit=ephem.Date(np.average(transit_times[~nanmask],weights=weights[~nanmask]))
    else: 
        meantransit=ephem.Date(np.nanmean(transit_times))
    if wrap_date==True and meantransit.datetime().day > ephem.Date(approximate_time).datetime().day:
        #This is just a quick & dirty hack.  Proper way is to re-calculate each on the 
        #transit times starting from the previous date if the mean really goes past midnight.
        meantransit-=1
    return str(ephem.Date(meantransit))

def calculate_targets_mean_transit_time(target_list, observer, approximate_time, weights=None, mode='nearest'):
    """
    Calculates the mean transit time of multiple targets at the supplied single observatory, for a specified approximate date+time
    
    Parameters
    ----------
    target_list : list [ elements being ephem.FixedBody(), ephem.Sun(), or ephem.Moon() ]
        Sky targets of interest
    observer : ephem.Observer() 
        The observatory in ephem.Observer() format
    approximate_time : str, formatted as 'YYYY/MM/DD HH:MM:SS.s'
        String UTC date/time, to use for calculating the transit times with ephem
    weights : np.array, list, or other array-like
        Optional (float) weights to give to certain stations.  (e.g. to prioritize certain baselines)   Larger values give higher importance. 
    mode : str
        Specifies which transit time to return. Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    
    Returns
    -------
    mean_transit_time : str
        The mean of the transit times of the mutliple specified targets
    
    Notes
    -----
    Weighted mean = Sum(w_i*x_i)/Sum(w_i)
    
    Examples
    --------
    mrk348=obs.create_ephem_target('Mrk348','00:48:47.14','31:57:25.1') \n
    ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751') \n
    ngc1068=obs.create_ephem_target('NGC1068','02:42:40.71','-00:00:47.81')  \n
    obs.calculate_targets_mean_transit_time([mrk348,ngc1052,ngc1068],obs.vlbaFD,'2025/09/15 12:00:00') \n
    ## Now, upweight ngc1052 and ngc1068 \n
    obs.calculate_targets_mean_transit_time([mrk348,ngc1052,ngc1068],obs.vlbaFD,'2025/09/15 12:00:00',weights=[1,5,2]) #--> '2025/9/15 09:46:05'
    """
    if weights is not None:
        if len(target_list)!=len(weights): 
            raise Exception('weights must be None or an array the same length as target_list array')
    
    tmp_ant=observer.copy() #Make a copy, so it doesn't modify the original object outside of this function...
    tmp_ant.date=ephem.Date(approximate_time)
    
    transit_times=np.zeros(len(target_list))
    for i in range(len(target_list)):
        tmp_tar=target_list[i].copy()
        tmp_tar.compute(tmp_ant)
        #transit_times[i]=tmp_tar.transit_time #deprecated use
        transit_times[i]=calculate_transit_time_single(tmp_tar,tmp_ant,approximate_time,mode=mode)
    
    if weights is not None: return str(ephem.Date(np.average(transit_times,weights=weights)))
    else: return str(ephem.Date(np.mean(transit_times)))

def calculate_VLBA_mean_transit_time(target, approximate_time, weights=None, mode='nearest'):
    """
    Calculates the mean transit time of a single target for the VLBA stations, for a specified approximate date+time. This is a convenience function which calls  calculate_vlbi_mean_transit_time()  for all VLBA antennae.
    
    Parameters
    ----------
    target : ephem.FixedBody() 
        Sky target of interest. 
    approximate_time : str, formatted as 'YYYY/MM/DD HH:MM:SS.s'
        String UTC date/time, to use for calculating the transit times with ephem
    weights : np.array, list, or other array-like
        Optional (float) weights to give to certain stations.  (e.g. to prioritize certain baselines)   Larger values give higher importance. The order of stations is alphabetical: [BR,FD,HN,KP,LA,MK,NL,OV,PT,SC].
    mode : str
        Specifies which transit time to return. Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    
    Returns
    -------
    mean_transit_time : str
        The mean of the transit times of the target at the various VLBA stations
    
    Notes
    -----
    Weighted mean = Sum(w_i*x_i)/Sum(w_i)
    
    Example
    -------
    mrk348=obs.create_ephem_target('Mrk348','00:48:47.14','31:57:25.1') \n
    obs.calculate_VLBA_mean_transit_time(mrk348,'2025/09/15 12:00:00')
    """
    return calculate_vlbi_mean_transit_time(target,[vlbaBR,vlbaFD,vlbaHN,vlbaKP,vlbaLA, vlbaMK,vlbaNL,vlbaOV,vlbaPT,vlbaSC], approximate_time, weights=weights, force_after_observer_i=9, wrap_date=True, mode=mode)

def calculate_optimal_VLBAstarttime(target,approximate_time,duration_hours,weights=None, return_fmt='%Y/%m/%d %H:%M:%S', LST_PT=False, mode='nearest'):
    """
    Calculate the optimal VLBA observation start time (to maximize station uptime), for a given target and desired observation duration. \n
    This is estimated from the mean transit time for the full array, minus half the specified duration. \n
    Calls  calculate_VLBA_mean_transit_time()
    
    Parameters
    ----------
    target : ephem.FixedBody() 
        Sky target of interest
    approximate_time : ephem.Date(), datetime.datetime, or str formatted as 'YYYY/MM/DD HH:MM:SS.s' 
        UTC date/time
    duration_hours : float
        The desired duration of the full VLBA observation session, in hours
    weights : array-like [list, np.array...]
        Optional (float) weights to give to certain stations.  (e.g. to prioritize certain baselines)   Larger values give higher importance. The order of stations is alphabetical: [BR,FD,HN,KP,LA,MK,NL,OV,PT,SC].
    return_fmt : str ['ephem', 'dt', or dt strftime format string]
        The format in which the start time will be returned. Options are: \n
             'ephem' : ephem.Date format \n
             'dt' or 'datetime' : dt.datetime \n
             strftime format code : return a string of the date, using the format code to convert from dt.datetime using dt.datetime.strftime().  Default is '%Y/%m/%d %H:%M:%S', which prints as 'YYYY/MM/DD hh:mm:ss'
    LST_PT : bool
        True returns the time in PieTown LST (e.g., used for VLBA dynamic scheduling) 
    mode : str
        Specifies which transit time to use. Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    
    Returns
    -------
    start_time : ephem.Date, datetime.datetime, or str
    
    Example
    -------
    # Observing quasar on 2025/09/15 for a duration of 5.5 hours \n
    ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751') \n
    obs.calculate_optimal_VLBAstarttime(ngc1052,'2025/09/15',5.5)   \n
    #--> '2025/09/15 07:19:27' \n
    obs.calculate_optimal_VLBAstarttime(ngc1052,'2025/09/15',5.5,LST_PT=True)  \n
    #--> '2025/09/14 23:45:01'
    """
    
    meantransit=calculate_VLBA_mean_transit_time(target,approximate_time,weights=weights,mode=mode)
    halfduration=dt.timedelta(hours=duration_hours/2.)
    start_dt=dt.datetime.strptime(meantransit,'%Y/%m/%d %H:%M:%S')-halfduration 
    #The calculated optimal start time start_dt is now in dt.datetime format
    
    ### Change start time to PieTown LST (timezone = ) if LST_PT set to True
    #   --> tzwhere.tzwhere().tzNameAt(obs.vlbaPT.lat*180/np.pi, obs.wrap_pm180(obs.vlbaPT.lon*180/np.pi)) 
    #       Returns tz string as 'America/Denver' ('US/Mountain' is now deprecated)
    #   At any rate, Mountain Time is UTC-6 or 7 hours, depending on Daylight Savings
    if LST_PT==True:
        start_dt=compute_sidereal_time(vlbaPT,start_dt,as_type='datetime')
        ## Note - NOT local time as done with, e.g.  utc_to_local(start_dt,'America/Denver')
    
    if 'ephem' in return_fmt.lower():
        #return an ephem.Date() object
        return ephem.Date(start_dt)
    elif True in [i in return_fmt.lower() for i in ['dt','datetime']]: 
        #return a dt.datetime object (already in that format, so just return it directly)
        return start_dt
    else:
        #return a string of the start time, following the formatting in return_fmt (default = '%Y/%m/%d %H:%M:%S' )
        try: return dt.datetime.strftime(start_dt,return_fmt)
        except: raise Exception("calculate_optimal_VLBAstarttime(): input return_fmt '%s' not recognized by dt.datetime.strftime()"%(return_fmt))


def compute_yearly_target_data(target, observer, obsyear, timezone='auto', time_of_obs='night', peak_alt=False, local=True, mode='nearest'):
    """
    Computes a target source's transit & rise/set times from a specific observer site, over the course of a calendar year.
    
    Parameters
    ----------
    target : an ephem.FixedBody() 
        Target observation source of interest
    observer : ephem.Observer() 
        Observatory/antenna/location object to use for the calculations
    obsyear : ephem.Date(), datetime.datetime, or str formatted as 'YYYY'
        Year in which to calculate the source transit & rise/set times
    timezone : str or None
        Timezone (standard Olson names) for local time on upper x-axis.  Options: \n
            'none' or None : Do not add local time info to the upper axis
            string of a standard Olson database name [e.g., 'US/Mountain' or 'America/Chicago'] --  use this directly for dt calculations \n
            'auto' or 'calculate' : compute the timezone from the observer lat/lon using module tzwhere.  Takes a while (but is convenient!) 
    time_of_obs : str  ['noon', 'night', 'midnight', 'peak' or 'HH:MM:SS'] 
        Denotes how to generate the times to use for each day. Options: \n
            'noon' = calculate the values at 12:00 local time each day \n
            'night' or 'midnight' = use 23:59:59 local time each day (for night obs.  00:00:00 would be for the previous night.) \n
            'peak' = return values during daily peak altitude \n
            #'exact' = use the exact timesteps as np.linspace returns them
    peak_alt : bool
        True also returns the daily peak altitudes
    local : bool 
        True returns times in the local observer site time. False returns in UTC time.
    mode : str
        Specifies which transit time to return. Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    
    Returns
    -------
    transit_times,rise_times,set_times,obstime_alts : np.array,np.array,np.array,np.array
    
    Example
    -------
    ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751') 
    TRS_alts_2021=obs.compute_yearly_target_data(ngc1052, obs.vlbaBR, '2021', local=True, timezone='US/Pacific')
    """
    if type(obsyear)==str: yearstring=str(obsyear)
    elif type(obsyear)==dt.datetime: yearstring=obsyear.year
    elif type(obsyear)==ephem.Date: yearstring=obsyear.datetime().year
    else: raise Exception('Input obsyear="%s" not understood. Must be of type str, dt.datetime, or ephem.Date'%(obsyear))
    
    if 'noon' in time_of_obs.lower(): timestring='12:00:00'; hoursdelta=12.
    elif 'night' in time_of_obs.lower(): timestring='23:59:59'; hoursdelta=23.99
    elif 'peak' in time_of_obs.lower(): timestring='00:00:00'; hoursdelta=0. #just use midnight for doing the peak altitude calcs for now, need to to update at some point
    else: 
        #User-specified time of day in format 'HH:MM:SS'
        timestring=str(time_of_obs); 
        hoursdelta=(dt.datetime.strptime(timestring,'%H:%M:%S')-dt.datetime.strptime('00:00:00','%H:%M:%S')).total_seconds()/3600.
    
    #if True in [i in timezone.lower() for i in ['auto','calc']]:
    #    #Use tzwhere to compute the timezone based on the observer lat/lon (input in degrees)
    #    try: timezone=tzwhere.tzwhere().tzNameAt(observer.lat*180/np.pi, wrap_pm180(observer.lon*180/np.pi))
    #    except: timezone=tzwhere.tzwhere(forceTZ=True).tzNameAt(observer.lat*180/np.pi, wrap_pm180(observer.lon*180/np.pi), forceTZ=True)
    timezone = tz_from_observer(observer)
    
    ### In order to deal with leap years, better to create start and end dates, and let datetime handle the number of days 
    #yearstart=pytz.timezone(timezone).localize(dt.datetime.strptime('%s/01/01 %s'%(yearstring,timestring),'%Y/%m/%d %H:%M:%S'))
    #yearend=pytz.timezone(timezone).localize(dt.datetime.strptime('%s/12/31 %s'%(yearstring,timestring),'%Y/%m/%d %H:%M:%S'))
    yearstart=pytz.timezone(timezone).localize(dt.datetime.strptime('%s/01/01 00:00:00'%(yearstring),'%Y/%m/%d %H:%M:%S'))
    yearend=pytz.timezone(timezone).localize(dt.datetime.strptime('%s/12/31 00:00:00'%(yearstring),'%Y/%m/%d %H:%M:%S'))
    ndays=(yearend-yearstart).days+1
    dt_array_local=np.array([yearstart+dt.timedelta(days=d) for d in list(range(ndays))]) #doing this explicitly as may want later
    dt_array_utc=np.array([d.astimezone(pytz.utc) for d in dt_array_local])
    
    transit_times=np.zeros(ndays).astype(dt.datetime) #UTC transit times
    rise_times=np.zeros(ndays).astype(dt.datetime) #UTC rise times
    set_times=np.zeros(ndays).astype(dt.datetime) #UTC set times
    if peak_alt==True: peak_alts=np.zeros(ndays)
    else: obstime_alts=np.zeros(ndays)
    for obstime_utc,i in zip(dt_array_utc+dt.timedelta(hoursdelta),list(range(ndays))):
        #These calculations assume time in UTC. So convert.
        #alt,az=compute_target_altaz_single(target,observer,obstime_utc) 
        transit_times[i]=calculate_transit_time_single(target, observer, ephem.Date(obstime_utc), return_fmt='dt', mode=mode)
        rise_times[i],set_times[i]=calculate_rise_set_times_single(target, observer, ephem.Date(obstime_utc), mode=mode, return_fmt='dt', verbose=False)
        if peak_alt==True: peak_alts[i]=compute_target_altaz_single(target, observer, ephem.Date(transit_times[i]) )[0]
        else: obstime_alts[i]=compute_target_altaz_single(target, observer, ephem.Date(obstime_utc) )[0]
    
    if local==True:
        transit_times_local=np.array([utc_to_local(d,timezone) for d in transit_times])
        rise_times_local=np.array([utc_to_local(d,timezone,verbose=False) for d in rise_times])
        set_times_local=np.array([utc_to_local(d,timezone,verbose=False) for d in set_times])
        if peak_alt==True: return transit_times_local,rise_times_local,set_times_local,peak_alts
        else: return transit_times_local,rise_times_local,set_times_local,obstime_alts
    else:
        if peak_alt==True: return transit_times,rise_times,set_times,peak_alts
        else: return transit_times,rise_times,set_times,obstime_alts

def calculate_VLBA_dynamic_starttime_range(target,approximate_time,duration_hours,weights=None, return_fmt='%Y/%m/%d %H:%M:%S', elevation_limit_deg=15., LST_PT=False, mode='nearest', plotresults=False):
    """
    Estimate sensible time range to start VLBA observation for a given target, based on the
    specified observation duration.  This is quite useful, for example, for helping to 
    quickly determine a sensible LST_PT range to put in your .key file for dynamic scheduling.  \n
    This process for estimating this range starts by first calculating the mean source transit 
    over all VLBA antennae, and subtracting half the duration as a starting point.  \n
    If not all stations remain up for the duration, it returns this mean transit minus half 
    the duration for both start & end times. \n
    If all stations are up for at least the full duration though, it returns the time range
    between MK rising and SC setting. \n 
    If all stations are always up (circumpolar source), it returns the start & end of the day.
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The sky source of interest
    approximate_time : ephem.Date(), datetime.datetime, or str formatted as 'YYYY/MM/DD HH:MM:SS.s' 
        UTC date/time, not local time
    duration_hours : float
        The desired duration of the full VLBA observation session, in hours
    weights : np.array, list, or other array-like
        Optional (float) weights to give to certain stations.  (e.g. to prioritize certain baselines)   Larger values give higher importance. The order of stations is alphabetical: [BR,FD,HN,KP,LA,MK,NL,OV,PT,SC].
    return_fmt : str  ['ephem', 'mjd', 'dt', or dt strftime format string]
        The format in which the start time will be returned. Options are: \n
            'ephem': return times as pyephem.Date objects [NOTE: these floats are Dublin JD, not MJD] \n
            'mjd': return times as float Modified Julian Dates \n
            'dt': return times as datetime.datetime objects \n
            String with strptime % formatting: return string dates (e.g., '%Y/%m/%d %H:%M:%S' [default]) 
    elevation_limit_deg : float [degrees]
        The elevation limit at the observatory to use for rise/set calculations, in degrees
    LST_PT : bool
        True returns the time in PieTown LST (e.g., used for VLBA dynamic scheduling) 
    mode : str
        Specifies which transit time to return. Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    plotresults : bool
        Set to True to plot MK/SC elevations and starttimes.  [NOT YET IMPLEMENTED]
    
    Returns
    ------- 
    dynamic_starttime_range : list  [dynamic_obs_range_start, dynamic_obs_range_end]
    
    Example
    ------- 
    mrk348=obs.create_ephem_target('Mrk348','00:48:47.14','31:57:25.1') \n
    dynamicobs_starttime_range = obs.calculate_VLBA_dynamic_starttime_range(mrk348,ephem.Date('2021/07/15 20:36:00'), 3.5, weights=None, return_fmt='%Y/%m/%d %H:%M:%S', elevation_limit_deg=10., LST_PT=True, plotresults=False)  \n
    #--> ['2021/07/14 21:57:04', '2021/07/15 00:22:53']
    """
    
    meantransit_str=calculate_VLBA_mean_transit_time(target,approximate_time,weights=weights,mode=mode)
    meantransit=ephem.Date(meantransit_str)
    halfduration=dt.timedelta(hours=duration_hours/2.)
    start_dt=dt.datetime.strptime(meantransit_str,'%Y/%m/%d %H:%M:%S')-halfduration
    
    #Optimal time (range) here means the target is visible from all stations.  
    #Range of times from when it rises at MK through when it sets at PT.
    #If that's not possible for the specified obs duration, just return the mean 
    #transit minus half the obs duration
    
    ## prev/next specific elevation set/rise times returned as [prev_set,prev_rise],[next_set,next_rise]
    setrisetimes_MK=calculate_antenna_visibility_limits(target,vlbaMK,meantransit, elevation_limit_deg=elevation_limit_deg, plusminusdays=1., interpsteps=100, verbose=False, alwaysup_fmt='limits', timeformat='ephem', LST_PT=LST_PT)
    setrisetimes_SC=calculate_antenna_visibility_limits(target,vlbaSC,meantransit, elevation_limit_deg=elevation_limit_deg, plusminusdays=1., interpsteps=100, verbose=False, alwaysup_fmt='limits', timeformat='ephem', LST_PT=LST_PT)
    
    allstationsup_hrs=(setrisetimes_SC[1][0]-setrisetimes_MK[0][1])*24
    
    if allstationsup_hrs>23.99:
        #If the target is always up / never sets, generate the obs day start/end times
        if LST_PT==True:
            start_dt_PT=compute_sidereal_time(vlbaPT,start_dt,as_type='datetime')
            ## Note - NOT to be confused with local time as done with, e.g.  
            ##  utc_to_local(start_dt,'America/Denver')
            currentdaystart_dt=start_dt_PT.replace(hour=0,minute=0,second=0)
            currentdayend_dt=start_dt_PT.replace(hour=23,minute=59,second=59)
        else: 
            currentdaystart_dt=start_dt.replace(hour=0,minute=0,second=0)
            currentdayend_dt=start_dt.replace(hour=23,minute=59,second=59)
        
        currentdaystart=ephem.Date(currentdaystart_dt)
        currentdayend=ephem.Date(currentdayend_dt)
    
    
    if allstationsup_hrs<duration_hours:
        meantransit_minus_half_duration=ephem.Date(calculate_optimal_VLBAstarttime(target, approximate_time, duration_hours, weights=weights, return_fmt=return_fmt, LST_PT=LST_PT, mode=mode)) 
        #--> Keep it in ephem.Date format for now; convert to return_fmt at the end.
        dynamic_starttime_range=[meantransit_minus_half_duration]*2
    else: 
        if allstationsup_hrs>23.99:
            dynamic_starttime_range=[currentdaystart,currentdayend]
        else: 
            dynamic_starttime_range=[setrisetimes_MK[0][1], ephem.Date(setrisetimes_SC[1][0]-duration_hours/24.)]
            #dynamic_starttime_range=[setrisetimes_MK[0][1], ephem.Date(np.nanmin([setrisetimes_SC[1][0]-duration_hours/24.,currentdayend])) ]
            #if setrisetimes_SC[1][0]==dt.datetime(start_dt.year,start_dt.month,start_dt.day,23,59,59):
    
    #print(dynamic_starttime_range)
    if return_fmt.lower()=='dt' or return_fmt.lower()=='datetime':
        for i in [0,1]: dynamic_starttime_range[i]=dynamic_starttime_range[i].datetime()
    elif '%' in return_fmt:
        for i in [0,1]: dynamic_starttime_range[i]=dynamic_starttime_range[i].datetime().strftime(return_fmt)
    elif 'mjd' in return_fmt.lower():
        for i in [0,1]: dynamic_starttime_range[i]=Time(dynamic_starttime_range[i].datetime()).mjd
    else: pass  #For any other inputs, leave as ephem.Date format
    
    return dynamic_starttime_range


def skysep_fixed_single(source1,source2):
    """
    Calculate the angular separation on the sky (in degrees), for fixed pyephem objects source1 and source2.  To calculate the angular distance between a source and the Sun or moon, use sunsep_single() or moonsep_single() functions instead.
    
    Parameters
    ----------
    source1 : ephem.FixedBody()
    source2 = ephem.FixedBody() 
    
    Returns
    -------
    angsep : float
        The angular separation between the sources, in degrees
    
    Example
    -------
    ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751') \n
    ngc3079=obs.create_ephem_target('NGC3079','10:01:57.80','55:40:47.24')  \n
    obs.skysep_fixed_single(ngc1052,ngc3079)  #--> 108.13770035565858  [degrees]
    """
    ### Manual calculation:
    #coords1_dec=[source1.ra*180/np.pi, source1.dec*180/np.pi] #decimal coordinates [in rad, conv. to degrees]
    #coords2_dec=[source2.ra*180/np.pi, source2.dec*180/np.pi] #decimal coordinates [in rad, conv. to degrees]
    #angsep=angulardistance(coords1_dec,coords2_dec) #,pythag_approx=False,returncomponents=False)
    
    ### Replaced manual calculations above with ephem builtin function for simplicity...
    angsep=ephem.separation(source1,source2)*180/np.pi #Separation [in rad] converted to degrees
    
    return angsep

##### Separation of target from Sun & Moon...

def moonsep_single(target,observer,obstime):
    """
    Calculate the angular separation between the Moon and a target, observed from a specified telescope/station and obsdate
    
    Parameters
    ----------
    target : ephem.FixedBody() 
        Target observation source
    observer : ephem.Observer() 
        Observatory/antenna location object to use for the calculations
    obstime : ephem.Date, datetime.datetime, str formatted as 'YYYY/MM/DD HH:MM:SS.s'
        The time at which the separation should be determined
    
    Returns
    -------
    moonsep : float
        The angular separation between the Moon and the source, in degrees
    """
    tmp_ant=observer.copy(); tmp_ant.date=ephem.Date(obstime)
    tmp_tar=target.copy(); tmp_tar.compute(tmp_ant)
    moon=ephem.Moon(); 
    moon.compute(tmp_ant); 
    moonsep=ephem.separation(moon,tmp_tar)*180/np.pi #Separation [in rad] converted to degrees
    return moonsep

def moonsep_timearray(target,observer,obstime_array):
    """
    Calculate an array of angular separations between the Moon and a target, observed from a specified telescope/station at multiple obsdate/times.  Useful, for example, for determining moon separation at each point in a night's observations.
    
    Parameters
    ----------
    target : ephem.FixedBody() 
        Target observation source
    observer : ephem.Observer() 
        Observatory/antenna location object to use for the calculations
    obstime_array : array-like [list,np.array...] 
        A series of date+times, in ephem.Date, dt.datetime, or string ('YYYY/MM/DD HH:MM:SS.s') 
    
    Returns
    -------
    moonseps : np.array
        The target separations from the Moon, in degrees, at each time in the series
    """
    tmp_ant=observer.copy(); tmp_ant.date=ephem.Date(obstime_array[0])   
    tmp_tar=target.copy(); tmp_tar.compute(tmp_ant) 
    moon=ephem.Moon();
    moonseps=np.zeros(len(obstime_array))
    for t in list(range(len(moonseps))):
        tmp_ant.date=obstime_array[t]
        moon.compute(tmp_ant);
        moonseps[t]=ephem.separation(moon,tmp_tar)*180/np.pi #Separation [in rad] converted to degrees
    return np.array(moonseps)

def moonsep_arrayminmax(target,observer_list,obstime,verbose=False,return_names=False):
    """
    Calculate the Moon separation from the target, observed from a specified obstime, and each supplied telescope/station, to determine the minimum/maximum separation 
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The target observation source
    observer_list : array-like [list, np.array, ...] 
        The list of ephem.Observer() antenna location objects to use for the calculations
    obstime : ephem.Date, datetime.datetime, str formatted as 'YYYY/MM/DD HH:MM:SS.s'
        The time at which the separation should be determined
    verbose : bool 
        True prints results to terminal
    return_names : bool
        Setting to True will also return the station names where the min/max occur
    
    Returns
    -------
    minmaxlist : list  [minsep,maxsep]
        The min/max target separation from the Moon, in degrees
    """
    moonseps=np.zeros(len(observer_list))
    for i in list(range(len(moonseps))):
        moonseps[i]=moonsep_single(target,observer_list[i],obstime)
    #The minimum & maximum separations
    minsep=np.nanmin(moonseps); maxsep=np.nanmax(moonseps)
    #The names of the stations/telescopes/observers corresponding to the min/max
    try: minname=observer_list[np.argmin(moonseps)].name
    except: minname= 'Observer_i='+str(np.argmin(moonseps))
    try: maxname=observer_list[np.argmax(moonseps)].name
    except: maxname= 'Observer_i='+str(np.argmax(moonseps))
    
    if verbose==True: print('  Min. separation from Moon = %.5f deg [%s], Max. sep. = %.5f deg [%s]'%(minsep,minname,maxsep,maxname))
    if return_names==True: return [minsep,maxsep],[minname,maxname]
    else: return [minsep,maxsep]

def moonsep_VLBAminmax(target,obstime,verbose=False,return_names=False):
    """
    Convenience function to call moonsep_arrayminmax() for the VLBA stations.
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The target observation source
    obstime : ephem.Date, datetime.datetime, str formatted as 'YYYY/MM/DD HH:MM:SS.s'
        The time at which the separation should be determined
    verbose : bool 
        True prints results to terminal
    return_names : bool
        Setting to True will also return the station names where the min/max occur
    
    Returns
    -------
    minmaxlist : list  [minsep,maxsep]
        The min/max target separation from the Moon, in degrees
    """
    return moonsep_arrayminmax(target,[vlbaBR,vlbaFD,vlbaHN,vlbaKP,vlbaLA,vlbaMK,vlbaNL,vlbaOV,vlbaPT,vlbaSC], obstime, verbose=verbose, return_names=return_names)

"""
For minimum safe separation of source from Sun for VLBA observations: see NRAO page
https://science.nrao.edu/facilities/vlba/docs/manuals/oss/obs-prep-ex
Freq (GHz)      0.33  0.61  1.6  2.3  5.0   8.4   15    22  43
Min.Sep (deg)   117   81    45   36   23    17    12    9   6
"""
def sunsep_single(target,observer,obstime):
    """
    Calculate the angular separation of the Sun from the target, observed from a specified telescope/station and obsdate+time
    
    Parameters
    ----------
    target : ephem.FixedBody() 
        Target observation source
    observer : ephem.Observer() 
        Observatory/antenna location object to use for the calculations
    obstime : ephem.Date, datetime.datetime, str formatted as 'YYYY/MM/DD HH:MM:SS.s'
        The time at which the separation should be determined
    
    Returns
    -------
    sunsep : float
        The angular separation between the Sun and the source, in degrees
    """
    tmp_ant=observer.copy(); tmp_ant.date=ephem.Date(obstime)
    tmp_tar=target.copy(); tmp_tar.compute(tmp_ant)
    sun=ephem.Sun(); 
    sun.compute(tmp_ant); 
    sunsep=ephem.separation(sun,tmp_tar)*180/np.pi #Separation [in rad] converted to degrees
    return sunsep

def sunsep_timearray(target,observer,obstime_array):
    """
    Calculate the angular sky separation of the Sun from the target, observed from a specified telescope/station at each obstime.  Useful, for example, for determining Sun separation at each point in an observing session.
    
    Parameters
    ----------
    target : ephem.FixedBody() 
        Target observation source
    observer : ephem.Observer() 
        Observatory/antenna location object to use for the calculations
    obstime_array : array-like [list,np.array...] 
        A series of date+times, in ephem.Date, dt.datetime, or string ('YYYY/MM/DD HH:MM:SS.s') 
    
    Returns
    -------
    sunseps : np.array
        The target separations from the Sun, in degrees, at each time in the series
    """
    tmp_ant=observer.copy(); tmp_ant.date=ephem.Date(obstime_array[0])   
    tmp_tar=target.copy(); tmp_tar.compute(tmp_ant) 
    sun=ephem.Sun();
    sunseps=np.zeros(len(obstime_array))
    for t in list(range(len(sunseps))):
        tmp_ant.date=obstime_array[t]
        sun.compute(tmp_ant);
        sunseps[t]=ephem.separation(sun,tmp_tar)*180/np.pi #Separation [in rad] converted to degrees
    return np.array(sunseps)

def sunsep_arrayminmax(target,observer_list,obstime,verbose=False,return_names=False):
    """
    Calculate the angular separation of the Sun from the target, observed from a specified obsdate, and each supplied telescope/station to get the minimum/maximum separation 
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The target observation source
    observer_list : array-like [list, np.array, ...] 
        The list of ephem.Observer() antenna location objects to use for the calculations
    obstime : ephem.Date, datetime.datetime, str formatted as 'YYYY/MM/DD HH:MM:SS.s'
        The time at which the separation should be determined
    verbose : bool 
        True prints results to terminal
    return_names : bool
        Setting to True will also return the station names where the min/max occur
    
    Returns
    -------
    minmaxlist : list  [minsep,maxsep]
        The min/max target separation from the Sun, in degrees
    """
    sunseps=np.zeros(len(observer_list))
    for i in list(range(len(sunseps))):
        sunseps[i]=sunsep_single(target,observer_list[i],obstime)
    #The minimum & maximum separations
    minsep=np.nanmin(sunseps); maxsep=np.nanmax(sunseps)
    #The names of the stations/telescopes/observers corresponding to the min/max
    try: minname=observer_list[np.argmin(sunseps)].name
    except: minname= 'Observer_i='+str(np.argmin(sunseps))
    try: maxname=observer_list[np.argmax(sunseps)].name
    except: maxname= 'Observer_i='+str(np.argmax(sunseps))
    
    if verbose==True: print('  Min. separation from Sun = %.5f deg [%s], Max. sep. = %.5f deg [%s]'%(minsep,minname,maxsep,maxname))
    if return_names==True: return [minsep,maxsep],[minname,maxname]
    else: return [minsep,maxsep]

def sunsep_VLBAminmax(target,obstime,verbose=False,return_names=False):
    """
    Convenience function to call sunsep_arrayminmax() for the VLBA stations.
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The target observation source
    obstime : ephem.Date, datetime.datetime, str formatted as 'YYYY/MM/DD HH:MM:SS.s'
        The time at which the separation should be determined
    verbose : bool 
        True prints results to terminal
    return_names : bool
        Setting to True will also return the station names where the min/max occur
    
    Returns
    -------
    minmaxlist : list  [minsep,maxsep]
        The min/max target separation from the Sun, in degrees
    """
    return sunsep_arrayminmax(target,[vlbaBR,vlbaFD,vlbaHN,vlbaKP,vlbaLA,vlbaMK,vlbaNL,vlbaOV,vlbaPT,vlbaSC], obstime, verbose=verbose, return_names=return_names)

def sunsep_allVLBA(target,obstime,verbose=False,return_names=False):
    """
    Calculates the angular separation from the sun for a given target and obsdate, for all VLBA stations.
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The target observation source
    obstime : ephem.Date, datetime.datetime, str formatted as 'YYYY/MM/DD HH:MM:SS.s'
        The time at which the separation should be determined
    verbose : bool 
        True prints results to terminal
    return_names : bool
        Setting to True will also return the station names where the min/max occur
    
    Returns
    -------
    sunseps : list
        List of the Sun separations (in degrees) as observed from all VLBA stations
    """
    sunseps=np.zeros(10)
    for i,station in zip(list(range(len(sunseps))),[vlbaBR,vlbaFD,vlbaHN,vlbaKP,vlbaLA,vlbaMK,vlbaNL,vlbaOV,vlbaPT,vlbaSC]):
        sunseps[i]=sunsep_single(target,station,obstime)
        if verbose==True:
            print('  Sun separation from %s = %f deg'%(station.name,sunseps[i]))
    if return_names==True: return sunseps,['BR','FD','HN','KP','LA','MK','NL','OV','PT','SC']
    else: return sunseps

def sunsep_VLBAminmax_timearray(target,obstime_array,verbose=False,return_names=False):
    """
    Calculate the min & max Sun separation for the VLBA stations, across an array of obstimes.
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The target observation source
    obstime_array : array-like [list,np.array...] 
        A series of date+times, in ephem.Date, dt.datetime, or string ('YYYY/MM/DD HH:MM:SS.s') 
    verbose : bool 
        True prints results to terminal
    return_names : bool
        Setting to True will also return the station names where the min/max occur
    
    Returns
    -------
    minmaxlist : list  [minsep,maxsep]
        The min/max target separation from the Sun from among the specified times, in degrees
    """
    sunseps=np.zeros((10,len(obstime_array)))
    for i,station in zip(list(range(10)),[vlbaBR,vlbaFD,vlbaHN,vlbaKP,vlbaLA,vlbaMK,vlbaNL,vlbaOV,vlbaPT,vlbaSC]):
        sunseps[i,:]=sunsep_timearray(target,station,obstime_array)
    
    #The minimum & maximum separations
    minsep=np.nanmin(sunseps); maxsep=np.nanmax(sunseps)
    #The names of the stations/telescopes/observers corresponding to the min/max
    minind=np.unravel_index(np.argmin(sunseps,axis=None), sunseps.shape); minstaind=minind[0]
    maxind=np.unravel_index(np.argmax(sunseps,axis=None), sunseps.shape); maxstaind=maxind[0]
    minname=[vlbaBR,vlbaFD,vlbaHN,vlbaKP,vlbaLA,vlbaMK,vlbaNL,vlbaOV,vlbaPT,vlbaSC][minstaind].name
    maxname=[vlbaBR,vlbaFD,vlbaHN,vlbaKP,vlbaLA,vlbaMK,vlbaNL,vlbaOV,vlbaPT,vlbaSC][maxstaind].name
    mintime=obstime_array[minind[1]]; maxtime=obstime_array[maxind[1]]; 
    
    if verbose==True: print('  Min. separation from Sun = %.5f deg [%s], at %s\n  Max. separation from Sun = %.5f deg [%s], at %s'%(minsep,minname,mintime,maxsep,maxname,maxtime))
    if return_names==True: return [minsep,maxsep],[minname,maxname],[mintime,maxtime]
    else: return [minsep,maxsep]


def optimal_visibility_date(target, observer, obsyear, extra_info=True, verbose=False, return_fmt='str', timezone='auto', time_of_obs='night', peak_alt=True, local=True):
    """
    Returns the optimal day/time to observe a specified target from a specified observer site, for the given year. (Based on the peak altitude at the time specified by time_of_obs.)
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The target observation source
    observer : ephem.Observer() 
        Observatory/antenna location object to use for the calculations
    obsyear : int, float, str, ephem.Date, or datetime.datetime
        The year to consider for the calculations.  Use 4-digit format 'YYYY' if given as a string.
    extra_info : bool 
        Set to True to also return a list of [rise time, set time, peak altitude] on that date
    verbose : bool 
        Set to True to also print results to screen
    return_fmt : str  ['string', 'dt', or 'ephem']
        The desired output format for the returned date.  
    timezone : str or None
        Timezone (standard Olson names) for local time on upper x-axis.  Options: \n
            'none' or None : Do not add local time info to the upper axis
            string of a standard Olson database name [e.g., 'US/Mountain' or 'America/Chicago'] :  use this directly for dt calculations \n
            'auto' or 'calculate' : compute the timezone from the observer lat/lon using module tzwhere.  Takes a while (but is convenient!)
    time_of_obs : str  ['noon', 'night', 'midnight', 'peak' or 'HH:MM:SS'] 
        Denotes how to generate the times to use for each day. Options: \n
            'middark' = calculate the altitudes at each night's middle-dark time \n
            'noon' = calculate the values at 12:00 local time each day \n
            'night' or 'midnight' = use 23:59:59 local time each day (for night obs.  00:00:00 would be for the previous night.) \n
            'peak' = return values during daily peak altitude \n
    
    Returns
    -------
    outdates : same as return_fmt, or np.array if extra_info=True
        Returns the optimal observing date/time in the format specified by return_fmt.  Optionally also returns a list of extra info comprised of [rise time, set time, peak altitude]
    
    Example
    -------
    ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751') \n
    obs.optimal_visibility_date(ngc1052,obs.vlbaBR,'2021',verbose=True,time_of_obs='12:00:00')  \n
    #--> '2021/09/19 03:46:38'
    """
    if time_of_obs=='middark': TRSP=compute_yearly_target_data(target, observer, obsyear, timezone=timezone, time_of_obs='23:59:59', peak_alt=True, local=local)
    else: TRSP=compute_yearly_target_data(target, observer, obsyear, timezone=timezone, time_of_obs=time_of_obs, peak_alt=True, local=local)
    #Caculate the daily alt values 
    if time_of_obs == 'peak':  pass
    elif time_of_obs == 'middark':
        suntimes=np.array([calculate_twilight_times(observer,t.replace(hour=12)) for t in TRSP[0]])[:,0].astype(ephem.Date)
        middarks=np.squeeze([ suntimes[:-1,1] + (suntimes[:-1,1]-suntimes[1:,0]) ])
        middarks=np.array([ephem.Date(m).datetime() for m in middarks])
        alts=np.array([compute_target_altaz_single(target,observer,m) for m in middarks])[:,0]
    else:
        if type(obsyear) in [str,float,int]: yearint=int(obsyear)
        elif type(obsyear) is dt.datetime: yearint=obsyear.year
        elif type(obsyear) is ephem.Date: yearint=obsyear.datetime().year
        else: raise Exception('Input obsyear="%s" not understood. Must be of type int, float, str, dt.datetime, or ephem.Date'%(obsyear))
        if 'night' in time_of_obs.lower(): 
            obsstart=dt.datetime(yearint,1,1,23,59,59); obsend=dt.datetime(yearint,12,31,23,59,59); 
        elif 'noon' in time_of_obs.lower(): 
            obsstart=dt.datetime(yearint,1,1,12,0,0); obsend=dt.datetime(yearint,12,31,12,0,0); 
        else: 
            obsstart=dt.datetime(yearint,1,1,int(time_of_obs[0:2]), int(time_of_obs[3:5]), int(time_of_obs[6:8]));
            obsend=dt.datetime(yearint,12,31,int(time_of_obs[0:2]), int(time_of_obs[3:5]), int(time_of_obs[6:8]));
        alts=compute_target_altaz(target,observer,ephem.Date(obsstart),ephem.Date(obsend),nsteps=365)[0]
    
    if time_of_obs=='peak':
        optind=np.argmax(TRSP[3]); optdate=TRSP[0][optind]
    elif time_of_obs=='middark': 
        optind=np.argmax(alts); optdate=ephem.Date(middarks[optind]).datetime()
    else: 
        optind=np.argmax(alts); 
        ts=create_obstime_array(obsstart, obsend, timezone_string='UTC', n_steps=365)
        optdate=ts[optind]
    
    sunsep=sunsep_single(target,observer,TRSP[0][optind])
    moonsep=moonsep_single(target,observer,TRSP[0][optind])
    moonperc=compute_moonphase(TRSP[0][optind],return_fmt='perc')
    
    
    if verbose==True:
        if target.name is '': target.name='Target'
        if observer.name is '': observer.name='Observer'
        print('\nOptimal observing date for %s, from %s, in year %s:'%(target.name,observer.name,obsyear))
        print('  %s  with transit occurring at %s %s time'%(TRSP[0][optind].strftime('%Y/%m/%d'),TRSP[0][optind].strftime('%H:%M:%S'),['local' if local==True else 'UTC'][0]))
        try: risestring=TRSP[1][optind].strftime('%H:%M:%S')
        except: risestring='Always down' #Case where TRSP[1][optind] is np.nan
        try: setstring=TRSP[2][optind].strftime('%H:%M:%S')
        except: setstring='Always up' #Case where TRSP[2][optind] is np.nan
        print('  On that date, rise time = %s, set time = %s, peak altitude = %.1f deg'%(risestring, setstring, TRSP[3][optind]))
        print('  At transit, separation from Sun = %i deg, Moon separation = %i deg, Moon = %i%% illuminated\n'%(sunsep,moonsep,moonperc))
    
    if 'str' in return_fmt.lower():
        #Return dates in string format
        ##outdates=[d.strftime('%Y/%m/%d %H:%M:%S') for d in np.array(TRSP)[:3,optind]] #--> raises exceptions when including NaNs (for dates when always up/down)
        outdates=[]
        for d in np.array(TRSP)[:3,optind]:
            try: outdates.append( d.strftime('%Y/%m/%d %H:%M:%S') )
            except: outdates.append( np.nan )
        #peakalt='{0:.2f}'.format(TRSP[3][optind]);
        sunsep='{0:.2f}'.format(sunsep); moonsep='{0:.2f}'.format(moonsep); moonperc='{0:.2f}'.format(moonperc); 
    elif 'ephem' in return_fmt.lower():
        #pyephem.Date format
        ##outdates=[ephem.Date(local_to_utc(d)) for d in np.array(TRSP)[:3,optind]] #--> this doesn't handle NaNs for dates when target always up/down
        outdates=[]
        for d in np.array(TRSP)[:3,optind]:
            try: outdates.append( ephem.Date(local_to_utc(d)) )
            except: outdates.append( np.nan )
    else: 
        #return dates as dt.datetime format
        outdates=np.array(TRSP)[:3,optind]
    
    if extra_info==True: return np.append(outdates,[TRSP[3][optind],sunsep,moonsep,moonperc]) #transit, rise, set, peak alt., sun separation, moon separation, moon percent illumination
    else: return outdates[0] #transit time (dt format)



def print_VLBA_observability_between_dates(obstarget, duration_hrs, ephemdatestart, ephemdateend, obstypestring='[]-Band', outtimefmt='%Y/%m/%d %H:%M:%S', mode='nearest'):
    """
    Prints out optimal VLBA observing start time (in UT and PT_LST) for each day between two given dates.
    
    Parameters
    ----------
    obstarget = ephem.FixedBody()
        The sky target of interest.  Initialize, e.g., as a create_ephem_target() object
    duration_hrs : float
        Duration of total observation session in hours 
    ephemdatestart : ephem.Date()
        The starting date/time to print results for
    ephemdateend : ephem.Date() 
        The end date/time to print results for
    obstypestring : str 
        Label denoting observation type, for printout.  (e.g., 'K-band')
    outtimefmt : str
        The dt.datetime.strftime format to use for the output times
    mode : str
        Specifies which transit time to return. Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    
    Example
    -------
    ngc2992=obs.create_ephem_target('NGC2992','09:45:42.05','-14:19:34.98')  \n
    obs.print_VLBA_observability_between_dates(ngc2992, 7.5, ephem.Date('2021/08/01 00:00:00'), ephem.Date('2021/08/31 23:59:59'), obstypestring='K-Band') \n
    # --> \n
    # Aug 01  --  2021/08/01 16:19:42 UTC,  2021/08/01 05:49:11 LST \n
    # Aug 02  --  2021/08/02 16:15:46 UTC,  2021/08/02 05:49:11 LST \n
    # ...
    """
    print('\nOptimal VLBA observation start times between %s and %s,\nfor %.2f hr duration %s obs. of %s:\n%s'%(ephemdatestart.datetime().strftime('%Y %B %d'),ephemdateend.datetime().strftime('%Y %B %d'),duration_hrs,obstypestring,obstarget.name,'-'*65))
    for djd in np.arange(ephemdatestart,ephemdateend):
        obsdate_i=ephem.Date(djd)
        obstime_ut=calculate_optimal_VLBAstarttime(obstarget,obsdate_i,duration_hrs,weights=None, return_fmt=outtimefmt, LST_PT=False, mode=mode)
        obstime_lst=calculate_optimal_VLBAstarttime(obstarget,obsdate_i,duration_hrs,weights=None, return_fmt=outtimefmt, LST_PT=True, mode=mode)
        print('  %s  --  %s UTC,  %s LST'%(obsdate_i.datetime().strftime('%b %d'), obstime_ut,  obstime_lst ))
        
def query_object_coords_simbad(stringname, return_fmt='dec'):
    """
    Query Simbad for coordinates of a named target.
    
    Note
    ----
    Internet connection required to complete the query.
    
    Parameters
    ----------
    stringname : str
        The common name of the desired query target. e.g., 'NGC1275'
    return_fmt : str
        'dec','decimal','sex', or 'sexagesimal'.  The desired format of the returned coordinates.
    
    Returns
    -------
    targetcoords : list
        List of the [RA,DEC] of the queried target 
    
    Example
    -------
    obs.query_object_coords_simbad('M1', return_fmt='dec') \n
    # --> -->  [83.63308333333333, 22.0145] \n
    obs.query_object_coords_simbad('NGC1275',return_fmt='sex') \n
    # -->  ['03 19 48.1597', '+41 30 42.114']
    """
    #from astroquery.simbad import Simbad
    querytable=Simbad.query_object(stringname)
    targetcoords_sex=[querytable[0]['RA'],querytable[0]['DEC']] #RA in hms, DEC in dms
    if 'dec' in return_fmt.lower(): 
        targetcoords_dec=sex2dec(*targetcoords_sex)
        return targetcoords_dec
    elif 'sex' in return_fmt.lower():
        return targetcoords_sex


def get_visible_targets_from_source_list(observer, time_start, time_end, source_list, elevation_limit=15., decbin_limits_deg=[-90,90], minimum_observability_minutes='full', coord_format='dec', skip_header=0, delimiter=',', nsteps=100, decimal_format='deg'):
    """
    For a specified observatory location, determine the observable targets from a list of sources.  Observability here means the targets will be above the specified telescope elevation limit for the full (or partial, if specified) duration.  Targets can be restricted to a certain declination range using a list of min/max limits.
    
    Note
    ----
    The format of the input source list is three columns: \n
    source name (string), RA, DEC \n
    The coordinates can be in decimal (default) or sexagesimal, to be specified by the user.\n
    Sparse lists (missing elements in a row) are not supported. \n
    Querying positions by name is not yet implemented within this function -- do this beforehand, e.g., with obs.query_object_coords_simbad()
    
    Parameters
    ----------
    observer : ephem.Observer() 
        Observatory/antenna location object to use for the calculations
    time_start : ephem.Date(), datetime.datetime, or str formatted as 'YYYY/MM/DD HH:MM:SS.s' 
        The observation start date/time, in UTC
    time_end : ephem.Date(), datetime.datetime, or str formatted as 'YYYY/MM/DD HH:MM:SS.s' 
        The observation end date/time, in UTC
    source_list : array-like sequence (list, tuple, ...)
        The list of candidate sources and coordinates, of shape Nx3: each row/element consisting of [Sourcename,RA,DEC]
    elevation_limit : float [degrees]
        The elevation limit at the observatory to use for visibility calculations, in degrees
    decbin_limits_deg : sequence (list, tuple, np.array...)
        The [minimum, maximum] declination limits for returned targets, in degrees.  
    minimum_observability_minutes : float, int,  or str 
        The minimum time (in minutes, as a float or int) that the source targets must be observable in the specified observation window. If set to 'full' (default), only returns source targets that are observable for the entire duration of the specified observation window.
    coord_format = 'dec'/'decimal' or 'sex'/'sexagesimal'.  Specifies the format of the source_list coordinates
    skip_header : int
        The number of header lines to skip when loading a source list from file, passed to np.genfromtxt
    delimiter : str
        The delimiter to use when loading a source list from file, passed to np.genfromtxt
    nsteps : int
        The number of steps over which to calculate altitudes for the visibility checks
    decimal_format : str
        'deg' or 'rad' -- The format of the source list coordinates (degrees or radians), when they are supplied as decimal (floats) instead of strings.  Passed to obs.create_ephem_target.
    
    Returns
    -------
    observable_targets : list
        A list of the sources that satisfy the minimum observability duration within the specified observing window.  
    
    Example
    -------
    obs.get_visible_targets_from_source_list([obs.vlbaMK], dt.datetime(2021,10,31,10,0,0), dt.datetime(2021,10,31,10,45,0), 'special_sources.txt', elevation_limit=15., minimum_observability_minutes=30., minimum_mutual_vis_observers='full', coord_format='sex')
    """
    
    ### Todo: implement option for local time.
    #   Add option to auto- calculate tz as in other functions
    #   First, convert times to ephem, and cast as datetime, to ensure all inputs end up as dt
    #   Next, attach the timezone to the datetime.
    #   Finally, dtaware_to_ephem() to get back in ephem.Date format, but with tz accounted for... 
    
    ephemdatestart=ephem.Date(time_start); ephemdateend=ephem.Date(time_end)

    if type(source_list) == str:
        # Load in the source list to a numpy array, using genfromtxt
        if 'dec' in coord_format.lower(): coordtype=None #coordtype=float
        else: coordtype=str
        source_list_arr=np.genfromtxt(source_list, dtype=coordtype, skip_header=skip_header, delimiter=delimiter)
    else: source_list_arr=np.array(source_list)
    
    # For each source in the source list, generate a visibility (elevation) track for
    # the specified observation duration, and include it if it satisfies the observability
    # duration.    
    observable_targets=[]
    decbin_limits_rad=[dbld*np.pi/180 for dbld in decbin_limits_deg]
    #with tqdm(total=source_list_arr.shape[0]) as pbar:
    for i in tqdm( range(source_list_arr.shape[0]), desc='Calculating visibility of target ', ascii=True ):
        target_i=create_ephem_target(*source_list_arr[i])
        if target_i.dec<decbin_limits_rad[0] or target_i.dec>decbin_limits_rad[1]: continue
        alts=compute_target_altaz(target_i,observer,ephemdatestart,ephemdateend,nsteps=nsteps)[0]
        above_min_alt=np.array(alts>elevation_limit)
        if np.nansum(above_min_alt)==0:
            #In this case, the elevation never rises above the minimum limit
            continue
        #if minimum_observability_minutes.lower()=='full':
        if np.nansum(~above_min_alt)==0:
            #All altitudes are above the limit
            observable_targets.append(source_list_arr[i])
        else: 
            #There is a mix of alts above & below the limit
            if type(minimum_observability_minutes)==str and minimum_observability_minutes.lower()=='full': 
                #Already checked above if target is always visible, and here it's not.
                continue
            
            #Generate an array of ephem times spanning the duration of the observation window
            obstimes=create_obstime_array(str(ephemdatestart),str(ephemdateend),n_steps=100)
            #The duration of visibility ('up time') for the target at this station, in dt.timedelta format
            td=obstimes[np.where(above_min_alt)[0][-1]]-obstimes[np.where(above_min_alt)[0][0]]
            
            if td.total_seconds()/60 >= minimum_observability_minutes:
                observable_targets.append(source_list_arr[i])
        
    return observable_targets
        

def get_visible_targets_from_source_list_multistation(observerlist, time_start, time_end, source_list, elevation_limit=15., decbin_limits_deg=[-90,90], minimum_observability_minutes='full', minimum_mutual_vis_observers='full', coord_format='dec', skip_header=0, delimiter=',', nsteps=100, decimal_format='deg'):
    """
    For a specified list of observatory locations, determine the observable targets from a list of sources.  Observability here means the targets will be above the specified telescope elevation limit for the full (or partial, if specified) duration, and concurrently visible across at least the specified number of stations. 
    
    Parameters
    ----------
    observer : ephem.Observer() 
        Observatory/antenna location object to use for the calculations
    time_start : ephem.Date(), datetime.datetime, or str formatted as 'YYYY/MM/DD HH:MM:SS.s' 
        The observation start date/time, in UTC
    time_end : ephem.Date(), datetime.datetime, or str formatted as 'YYYY/MM/DD HH:MM:SS.s' 
        The observation end date/time, in UTC
    source_list : array-like sequence (list, tuple, ...)
        The list of candidate sources and coordinates, of shape Nx3: each row/element consisting of [Sourcename,RA,DEC]
    elevation_limit : float [degrees]
        The elevation limit at the observatory to use for visibility calculations, in degrees
    decbin_limits_deg : sequence (list, tuple, np.array...)
        The [minimum, maximum] declination limits for returned targets, in degrees.  
    minimum_observability_minutes : float, int, or str 
        The minimum time (in minutes, as a float or int) that the source targets must be observable in the specified observation window. If set to 'full' (default), only returns source targets that are observable for the entire duration of the specified observation window.
    minimum_mutual_vis_observers : int, or str 
        The minimum number of stations across which a target must be concurrently visible.  If set to 'full' (default), only returns source targets that are observable over the entire array, for the specified duration.
    coord_format = 'dec'/'decimal' or 'sex'/'sexagesimal'.  Specifies the format of the source_list coordinates
    skip_header : int
        The number of header lines to skip when loading a source list from file, passed to np.genfromtxt
    delimiter : str
        The delimiter to use when loading a source list from file, passed to np.genfromtxt
    nsteps : int
        The number of steps over which to calculate altitudes for the visibility checks, and times for the duration checks
    decimal_format : str
        'deg' or 'rad' -- The format of the source list coordinates (degrees or radians), when they are supplied as decimal (floats) instead of strings.  Passed to obs.create_ephem_target.
    
    Returns
    -------
    observable_targets : list
        A list of the sources that satisfy the minimum observability duration within the specified observing window, with mutual visibility among at least the minimum number of stations.  
    
    Example
    -------
    obs.get_visible_targets_from_source_list_multistation([obs.vlbaMK,obs.vlbaHN], dt.datetime(2021,10,31,10,0,0), dt.datetime(2021,10,31,10,45,0), 'special_sources.txt', elevation_limit=15., minimum_observability_minutes=30., minimum_mutual_vis_observers='full', coord_format='sex')
    """
    ephemdatestart=ephem.Date(time_start); ephemdateend=ephem.Date(time_end)
    
    if type(minimum_mutual_vis_observers)==str and minimum_mutual_vis_observers.lower()=='full': 
        minimum_mutual_vis_observers=len(observerlist)
                
    if type(source_list) == str:
        # Load in the source list to a numpy array, using genfromtxt
        if 'dec' in coord_format.lower(): coordtype=None #coordtype=float
        else: coordtype=str
        source_list_arr=np.genfromtxt(source_list, dtype=coordtype, skip_header=skip_header, delimiter=delimiter)
    else: source_list_arr=np.array(source_list)
    
    # For each source in the source list, generate a visibility (elevation) track for
    # the specified observation duration -- at each observer -- and include the source
    # if it satisfies the observability duration and the minimum mutual visibility.    
    observable_targets=[]
    decbin_limits_rad=[dbld*np.pi/180 for dbld in decbin_limits_deg]
    #with tqdm(total=source_list_arr.shape[0]) as pbar:
    
    for i in tqdm( range(source_list_arr.shape[0]), desc='Calculating visibility of target ', ascii=True ):
        target_i=create_ephem_target(*source_list_arr[i])
        if target_i.dec<decbin_limits_rad[0] or target_i.dec>decbin_limits_rad[1]: continue
        above_min_alts_all=[]
        for station in observerlist:
            alts=compute_target_altaz(target_i,station,ephemdatestart,ephemdateend,nsteps=nsteps)[0]
            above_min_alt=np.array(alts>elevation_limit)
            above_min_alts_all.append(above_min_alt)
        above_min_alts_all=np.array(above_min_alts_all)
        #--> above_min_alts_all.shape is [Nobservers,Ntimesteps] 
        
        #...Now do some calcs on above_min_alts_all, to see if a minimum number of stations (i.e. having mutual visibility) satisfy the minimum duration...
        if np.nansum(~above_min_alts_all)==0:
            #All altitudes are above the limit, at all stations (and therefore also mutually visible by all)
            observable_targets.append(source_list_arr[i])
        else: 
            #There is a mix of alts above & below the limit
            if type(minimum_observability_minutes)==str and minimum_observability_minutes.lower()=='full': 
                #Already checked above if target is always visible, and here it's not.
                continue
            
            #Generate an array of ephem times spanning the duration of the observation window
            obstimes=create_obstime_array(str(ephemdatestart),str(ephemdateend),n_steps=nsteps)
            
            #The number of stations with mutual visibility at each time step
            mutual_vis_number=np.sum(above_min_alts_all,axis=0)
            satisfies_mutual_vis_number = mutual_vis_number >= minimum_mutual_vis_observers
            
            #Appling the minimum mutual vis boolean filter to the array specifying above minimum altitude
            above_min_alts_mutual=above_min_alts_all*satisfies_mutual_vis_number
            
            #The duration of visibility ('up time') at each station where the minimum number
            # of stations with mutual visibility is satisfied
            #td = np.array( [obstimes[np.where(arr)[0][-1]] - obstimes[np.where(arr)[0][0]] for arr in above_min_alts_mutual] )
            td = []
            for arr in above_min_alts_mutual:
                try:
                    td.append( obstimes[np.where(arr)[0][-1]] - obstimes[np.where(arr)[0][0]] )
                except: 
                    td.append(dt.timedelta(-1))
            
            #Determine whether the minimum duration is satisfied (with mutual visibility 
            # accounted for), for each station
            satisfies_duration = np.array([t.total_seconds() for t in td])/60 >= minimum_observability_minutes
            
            #Finally, add the target if the number of stations for which the minimum mutual
            #visibility duration is satisfied is >= the required number
            if np.sum(satisfies_duration) >= minimum_mutual_vis_observers:
                observable_targets.append(source_list_arr[i])
    
    return observable_targets


def get_visible_targets_from_source_list_VLBA(time_start, time_end, source_list, elevation_limit=15., decbin_limits_deg=[-90,90], minimum_observability_minutes='full', minimum_mutual_vis_observers='full', fast_mode=True, coord_format='dec', skip_header=0, delimiter=',', nsteps=100, decimal_format='deg'):
    """
    Determine the observable targets for the VLBA from a list of sources.  Observability here means the targets will be above the specified telescope elevation limit for the full (or partial, if specified) duration, and concurrently visible across at least the specified number of stations. 
    
    Parameters
    ----------
    time_start : ephem.Date(), datetime.datetime, or str formatted as 'YYYY/MM/DD HH:MM:SS.s' 
        The observation start date/time, in UTC
    time_end : ephem.Date(), datetime.datetime, or str formatted as 'YYYY/MM/DD HH:MM:SS.s' 
        The observation end date/time, in UTC
    source_list : array-like sequence (list, tuple, ...)
        The list of candidate sources and coordinates, of shape Nx3: each row/element consisting of [Sourcename,RA,DEC]
    elevation_limit : float [degrees]
        The elevation limit at the observatory to use for visibility calculations, in degrees
    decbin_limits_deg : sequence (list, tuple, np.array...)
        The [minimum, maximum] declination limits for returned targets, in degrees.  
    minimum_observability_minutes : float, int, or str 
        The minimum time (in minutes, as a float or int) that the source targets must be observable in the specified observation window. If set to 'full' (default), only returns source targets that are observable for the entire duration of the specified observation window.
    minimum_mutual_vis_observers : int, or str 
        The minimum number of stations across which a target must be concurrently visible.  If set to 'full' (default), only returns source targets that are observable over the entire array, for the specified duration.
    fast_mode : bool
        If set to True, will perform calculations for the four outermost N,E,S,W stations (MK,BR,HN,SC) instead of the full array of 10 stations 
    coord_format = 'dec'/'decimal' or 'sex'/'sexagesimal'.  Specifies the format of the source_list coordinates
    skip_header : int
        The number of header lines to skip when loading a source list from file, passed to np.genfromtxt
    delimiter : str
        The delimiter to use when loading a source list from file, passed to np.genfromtxt
    nsteps : int
        The number of steps over which to calculate altitudes for the visibility checks, and times for the duration checks
    decimal_format : str
        'deg' or 'rad' -- The format of the source list coordinates (degrees or radians), when they are supplied as decimal (floats) instead of strings.  Passed to obs.create_ephem_target.
    
    Returns
    -------
    observable_targets : list
        A list of the sources that satisfy the minimum observability duration within the specified observing window, with mutual visibility among at least the minimum number of stations.  
    
    Example
    -------
    obs.get_visible_targets_from_source_list_VLBA('2021/10/31 10:00:00.0', '2021/10/31 10:45:00.0', './spooky_halloween_sources.txt', elevation_limit=15., minimum_observability_minutes='full', minimum_mutual_vis_observers='full', fast_mode=True, coord_format='dec')
    """
    if fast_mode == True:
        observerlist = [vlbaMK,vlbaBR,vlbaHN,vlbaSC]
    else:
        observerlist = [vlbaBR,vlbaFD,vlbaHN,vlbaKP,vlbaLA,vlbaMK,vlbaNL,vlbaOV,vlbaPT,vlbaSC]
    observable_targets = get_visible_targets_from_source_list_multistation(observerlist, time_start, time_end, source_list, elevation_limit=elevation_limit, decbin_limits_deg=decbin_limits_deg,  minimum_observability_minutes=minimum_observability_minutes, minimum_mutual_vis_observers=minimum_mutual_vis_observers, coord_format=coord_format, skip_header=skip_header, delimiter=delimiter, nsteps=nsteps)
    
    return observable_targets



##### ------------ Utility functions for visualizations/plotting ------------ #####



def plot_year_observability(target, observer, obsyear, timezone='auto', time_of_obs='midnight', local=True, figsize=(14,8), dpi=200, savepath='', showplot=False, xaxisformatter=mdates.DateFormatter('%b')):
    """
    Plot the daily peak visibility (altitude) track for a target source over the course of a year, from a given observer site. Similar to a classic 'starobs' plot.
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The target observation source
    observer : ephem.Observer() 
        Observatory/antenna location object to use for the calculations
    obsyear : int, float, str, ephem.Date, or datetime.datetime
        The year to consider for the calculations.  Use 4-digit format 'YYYY' if given as a string.
    timezone : str or None
        Timezone (standard Olson database names) for local time on upper x-axis.  Options: \n
            'none' or None : Do not add local time info to the upper axis
            string of a standard Olson database name [e.g., 'US/Mountain' or 'America/Chicago'] :  use this directly for dt calculations \n
            'auto' or 'calculate' : compute the timezone from the observer lat/lon using module tzwhere.  Takes a while (but is convenient!) 
    time_of_obs : str  ['noon', 'night', 'midnight', 'peak' or 'HH:MM:SS'] 
        Denotes how to generate the times to use for each day. Options: \n
            'middark' = calculate the altitudes at each night's middle-dark time \n
            'noon' = calculate the values at 12:00 local time each day \n
            'night' or 'midnight' = use 23:59:59 local time each day (for night obs.  00:00:00 would be for the previous night.) \n
            'peak' = return values during daily peak altitude \n
    local : bool 
        True returns times in the local observer site time. False returns in UTC time.
    figsize : tuple
        matplotlib Figure size
    dpi : int
        matplotlib image resolution (dpi = dots per inch)
    savepath : str
        The path for saving the figure, including the filename and extension. 
    showplot : bool
        Set to True to invoke plt.show() and display the interactive plot
    xaxisformatter: mdates.DateFormatter
        Formatter object for the xaxis dates.  mdates.DateFormatter('%b') will print abbreviated month names, '%B' will print full month names.
    
    Example
    -------
    crab=obs.create_ephem_target('Crab Nebula','05:34:31.94','22:00:52.2') #'M1' \n
    obs.plot_year_observability(crab, obs.vlbaMK, '2021', showplot=True)
    """    
    #Calculate the daily transit times, rise times, set times, and peak altitudes for the specified year
    if time_of_obs=='middark': TRSP=compute_yearly_target_data(target, observer, obsyear, timezone=timezone, time_of_obs='23:59:59', peak_alt=True, local=local)
    else: TRSP=compute_yearly_target_data(target, observer, obsyear, timezone=timezone, time_of_obs=time_of_obs, peak_alt=True, local=local)
    #Caculate the daily alt values
    if time_of_obs == 'peak':  pass
    elif time_of_obs == 'middark':
        suntimes=np.array([calculate_twilight_times(observer,t.replace(hour=12)) for t in TRSP[0]])[:,0].astype(ephem.Date)
        middarks=np.squeeze([ suntimes[:-1,1] + (suntimes[:-1,1]-suntimes[1:,0]) ])
        middarks=np.array([ephem.Date(m).datetime() for m in middarks])
        alts=np.array([compute_target_altaz_single(target,observer,m) for m in middarks])[:,0]
    else:
        if type(obsyear) in [str,float,int]: yearint=int(obsyear)
        elif type(obsyear) is dt.datetime: yearint=obsyear.year
        elif type(obsyear) is ephem.Date: yearint=obsyear.datetime().year
        else: raise Exception('Input obsyear="%s" not understood. Must be of type int, float, str, dt.datetime, or ephem.Date'%(obsyear))
        if 'night' in time_of_obs.lower(): 
            obsstart=dt.datetime(yearint,1,1,23,59,59); obsend=dt.datetime(yearint,12,31,23,59,59); 
        elif 'noon' in time_of_obs.lower(): 
            obsstart=dt.datetime(yearint,1,1,12,0,0); obsend=dt.datetime(yearint,12,31,12,0,0); 
        else: 
            obsstart=dt.datetime(yearint,1,1,int(time_of_obs[0:2]), int(time_of_obs[3:5]), int(time_of_obs[6:8]));
            obsend=dt.datetime(yearint,12,31,int(time_of_obs[0:2]), int(time_of_obs[3:5]), int(time_of_obs[6:8]));
        alts=compute_target_altaz(target,observer,ephem.Date(obsstart),ephem.Date(obsend),nsteps=365)[0]
    
    fig1=plt.figure(1,figsize=figsize)
    axin=fig1.add_subplot(111)#,adjustable='box-forced',aspect=5e-3)
    
    #sun_alts=compute_sun_tracks(observer,ephem.Date(TRSP[0][0]),ephem.Date(TRSP[0][-1]),nsteps=len(TRSP[0]))[0]
    #moon_alts=compute_moon_tracks(observer,ephem.Date(TRSP[0][0]),ephem.Date(TRSP[0][-1]),nsteps=len(TRSP[0]))[0]
    sunseps=sunsep_timearray(target,observer,TRSP[0])
    moonseps=moonsep_timearray(target,observer,TRSP[0])
    moonpercs=[compute_moonphase(ephem.Date(d),return_fmt='percent') for d in TRSP[0]][0]
    
    if time_of_obs=='peak': axin.plot(TRSP[0],TRSP[3],color='0.2',label='Peak altitude',zorder=5)
    elif time_of_obs=='middark': 
        axin.plot(middarks,alts,color='0.2',label='Middle-dark altitudes',zorder=5)
    else: 
        ts=create_obstime_array(obsstart, obsend, timezone_string='UTC', n_steps=365)
        axin.plot(ts,alts,color='0.2',label='Altitudes at '+time_of_obs,zorder=5)
    #axin.plot(TRSP[0],moon_alts,color='0.6',lw=0.8,ls='--',label='Moon',alpha=0.5,zorder=1)
    #axin.plot(TRSP[0],sun_alts,color='orange',lw=0.8,ls='--',label='Sun',alpha=0.5,zorder=1)
    axin.plot(TRSP[0],moonseps,color='0.6',lw=0.8,ls='--',label='Moon Separation (deg)',alpha=0.5,zorder=1)
    axin.plot(TRSP[0],sunseps,color='orange',lw=0.8,ls='--',label='Sun Separation (deg)',alpha=0.5,zorder=1)
    
    ### Vertical line denoting the optimal observing date
    if time_of_obs=='peak':
        optind=np.argmax(TRSP[3]); optdate=TRSP[0][optind]
    elif time_of_obs=='middark': 
        optind=np.argmax(alts); optdate=ephem.Date(middarks[optind]).datetime()
    else: 
        optind=np.argmax(alts); optdate=ts[optind]
    axin.axvline(optdate, ls='--',color='0.5',zorder=-1)
    axin.annotate('Optimal observing date: %s (%s)'%(optdate.strftime('%Y/%m/%d %H:%M:%S'),['local time' if local==True else 'UTC'][0]),xy=[optdate,85], xytext=[optdate+dt.timedelta(days=-7),85], va='top', rotation=90,color='0.5')
    
    axin.legend(loc='best')#loc='upper right')
    if time_of_obs=='peak': axin.set_xlim(TRSP[0][0],TRSP[0][-1]); 
    elif time_of_obs=='middark': axin.set_xlim(middarks[0],middarks[-1]);
    else: axin.set_xlim(obsstart,obsend);
    axin.set_ylim(0,90)
    axin.xaxis.set_major_locator(mdates.MonthLocator())
    axin.xaxis.set_major_formatter(ticker.NullFormatter())
    axin.xaxis.set_minor_formatter(xaxisformatter)
    axin.xaxis.set_minor_locator(mdates.MonthLocator(bymonthday=16))
    #axin.xaxis.set_major_formatter(xaxisformatter)
    #axin.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=16))
    for tick in axin.xaxis.get_minor_ticks():
        tick.tick1line.set_markersize(0); tick.tick2line.set_markersize(0); tick.label1.set_horizontalalignment('center')
    #fig1.autofmt_xdate(rotation=0,ha='left')
    axin.set_xlabel('Month in %i'%(TRSP[0][0].year)); 
    axin.set_ylabel('Altitude [deg]')
    plt.grid(True,color='0.92')#,'both')
    
    axin2=axin.twinx(); #axin2.set_ylim(0,90)
    axin2.set_ylim(axin.get_ylim()) #--> This actually returns slightly different lims!!  Not [0,90]
    axin2.set_yticks([10,30,50,70,90]); #Just compute and set them manually on the ALTITUDE frame
    axin2.set_yticklabels(np.round(alt2airmass(np.array([10,30,50,70,90])),decimals=2)) 
    axin2.set_ylabel('Airmass [sec(z)]',rotation=-90,va='bottom')
    
    ### Possibly add a second subplot, with N hours on y-axis, plotting number 
    ### of hours target is up during night, and the number of sunless hours -- 
    ### like in starobs plots.
    #axin3.plot(TRSP[0],moonseps,color='0.6',lw=0.8,ls='--',label='Moon (%i%%)'%(moonseps),alpha=0.5,zorder=1)
    #axin3.plot(TRSP[0],sunseps,color='r',lw=0.8,ls='--',label='Sun',alpha=0.5,zorder=1)
    
    if savepath is not '': plt.savefig(savepath,bbox_inches='tight',dpi=dpi)
    if showplot==True: plt.show()
    plt.clf(); plt.close('all')

def plot_year_RST(target, observer, obsyear, timezone='auto', time_of_obs='night', local=True, showsun=True, showmoon=False, figsize=(14,8), dpi=200, savepath='', showplot=False, xaxisformatter=mdates.DateFormatter('%b'), yaxisformatter=mdates.DateFormatter('%H:%M')):
    """
    For a specified target source and observer site, plot the time of day the target rises, sets, and transits over the year.
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The target observation source
    observer : ephem.Observer() 
        Observatory/antenna location object to use for the calculations
    obsyear : int, float, str, ephem.Date, or datetime.datetime
        The year to consider for the calculations.  Use 4-digit format 'YYYY' if given as a string.
    timezone : str or None
        Timezone (standard Olson names) for local time on upper x-axis.  Options: \n
            'none' or None : Do not add local time info to the upper axis
            string of a standard Olson database name [e.g., 'US/Mountain' or 'America/Chicago'] :  use this directly for dt calculations \n
            'auto' or 'calculate' : compute the timezone from the observer lat/lon using module tzwhere.  Takes a while (but is convenient!)
    time_of_obs : str  ['noon', 'night', 'midnight', 'peak' or 'HH:MM:SS'] 
        Denotes how to generate the times to use for each day. Options: \n
            'noon' = calculate the values at 12:00 local time each day \n
            'night' or 'midnight' = use 23:59:59 local time each day (for night obs.  00:00:00 would be for the previous night.) \n
            'peak' = return values during daily peak altitude \n
    local : bool 
        True returns times in the local observer site time. False returns in UTC time.
    showsun : bool
        Set to True to include a track on the plot for the Sun
    showmoon : bool
        Set to True to include a track on the plot for the Moon
    figsize : tuple
        matplotlib Figure size
    dpi : int
        matplotlib image resolution (dpi = dots per inch)
    savepath : str
        The path for saving the figure, including the filename and extension. 
    showplot : bool
        Set to True to invoke plt.show() and display the interactive plot
    xaxisformatter: mdates.DateFormatter
        Formatter object for the xaxis dates.  Default mdates.DateFormatter('%b') will print abbreviated month names, '%B' will print full month names.
    yaxisformatter: mdates.DateFormatter
        Formatter object for the yaxis times.  Default mdates.DateFormatter('%H:%M') will print the hour and minute.
    
    Example
    -------
    wht = ephem.create_ephem_observer('WHT', '-17 52 53.8', '28 45 37.7', 2344); \n
    crab = obs.create_ephem_target('Crab Nebula','05:34:31.94','22:00:52.2') #'M1' \n
    obs.plot_year_RST(crab, wht, '2021', showplot=True)
    """
    #if True in [i in timezone.lower() for i in ['auto','calc']]:
    #    #Use tzwhere to compute the timezone based on the observer lat/lon (input in degrees)
    #    try: timezone=tzwhere.tzwhere().tzNameAt(observer.lat*180/np.pi, wrap_pm180(observer.lon*180/np.pi))
    #    except: timezone=tzwhere.tzwhere(forceTZ=True).tzNameAt(observer.lat*180/np.pi, wrap_pm180(observer.lon*180/np.pi), forceTZ=True)
    timezone = tz_from_observer(observer)
    
    #Calculate the daily transit times, rise times, set times, and peak altitudes for the specified year
    TRSP=compute_yearly_target_data(target, observer, obsyear, timezone=timezone, time_of_obs=time_of_obs, peak_alt=True, local=local)
    
    fig1=plt.figure(1,figsize=figsize)
    axin=fig1.add_subplot(111)#,adjustable='box-forced',aspect=5e-3)
    
    dt2hrfrac=lambda dt_in: dt_in.hour+dt_in.minute/60.+dt_in.second/3600.
    hourfracs_rise=np.array([dt2hrfrac(t) for t in TRSP[1]])
    hourfracs_transit=np.array([dt2hrfrac(t) for t in TRSP[0]])
    hourfracs_set=np.array([dt2hrfrac(t) for t in TRSP[2]])
    
    #sun_alts=compute_sun_tracks(observer,ephem.Date(TRSP[0][0]),ephem.Date(TRSP[0][-1]),nsteps=len(TRSP[0]))[0]
    #moon_alts=compute_moon_tracks(observer,ephem.Date(TRSP[0][0]),ephem.Date(TRSP[0][-1]),nsteps=len(TRSP[0]))[0]
    #sunseps=sunsep_timearray(target,observer,TRSP[0])
    #moonseps=moonsep_timearray(target,observer,TRSP[0])
    #moonpercs=[compute_moonphase(ephem.Date(d),return_fmt='percent') for d in TRSP[0]][0]
    if showmoon==True: 
        moontimes=np.array([calculate_moon_times(observer,t,outtype='dt') for t in TRSP[0]]) #dt format, [Ndays,2]
        hourfracs_moonrise=np.array([dt2hrfrac(utc_to_local(t,timezone)) for t in moontimes[:,0]])
        hourfracs_moonset=np.array([dt2hrfrac(utc_to_local(t,timezone)) for t in moontimes[:,1]])
    if showsun==True:
        suntimes=np.array([calculate_twilight_times(observer,t.replace(hour=12)) for t in TRSP[0]])[:,0].astype(ephem.Date) #ephem.Date format, [sunset,sunrise]
        for i in range(suntimes.shape[0]): 
            suntimes[i,0]=ephem.Date(suntimes[i,0]).datetime(); suntimes[i,1]=ephem.Date(suntimes[i,1]).datetime();
        hourfracs_sunrise=np.array([dt2hrfrac(utc_to_local(t,timezone)) for t in suntimes[:,1]])
        hourfracs_sunset=np.array([dt2hrfrac(utc_to_local(t,timezone)) for t in suntimes[:,0]])
    
    axin.scatter(TRSP[1],hourfracs_rise,color='#3B5B92',label='Rise Time',zorder=5)
    axin.scatter(TRSP[0],hourfracs_transit,color='#4B006E',label='Transit Time',zorder=5)
    axin.scatter(TRSP[2],hourfracs_set,color='#A83C09',label='Set Time',zorder=5)
    
    #Now fill in the times between rise/transit and transit/set
    def fill_RTStimes(times_earlier,times_later,facecolor='k',alpha=0.3,lw=0.0,zorder=0):
        splitmask=(times_earlier<times_later)
        #switchinds=np.where(np.diff(splitmask)) #indices where switches happen (True to False or False to True)
        #switchinds=(np.diff(splitmask) !=0)[0] 
        switchinds=np.array([0,]+list(np.unique(np.where(np.diff(splitmask))[0]+1))+[len(times_earlier),])
        chunkinds=[[switchinds[i],switchinds[i+1]] for i in list(range(len(switchinds)-1))]
        for chunk in chunkinds:
            chunkmask=np.array([False,]*len(splitmask)); chunkmask[chunk[0]:chunk[1]]=True
            if splitmask[chunk[0]]==True:
                axin.fill_between(TRSP[1][chunkmask],times_earlier[chunkmask],y2=times_later[chunkmask], facecolor=facecolor, alpha=alpha, lw=lw, zorder=zorder)
            else: 
                axin.fill_between(TRSP[1][chunkmask],times_later[chunkmask],y2=np.zeros(np.sum(chunkmask)), facecolor=facecolor, alpha=alpha, lw=lw, zorder=zorder)
                axin.fill_between(TRSP[1][chunkmask],times_earlier[chunkmask],y2=np.zeros(np.sum(chunkmask))+24, facecolor=facecolor, alpha=alpha, lw=lw, zorder=zorder)
    
    fill_RTStimes(hourfracs_rise,hourfracs_transit,facecolor='#3B5B92',alpha=0.3,lw=0.0,zorder=3)
    fill_RTStimes(hourfracs_transit,hourfracs_set,facecolor='#A83C09',alpha=0.3,lw=0.0,zorder=3)
    
    if showmoon==True: 
        axin.scatter(TRSP[1],hourfracs_moonrise,color='0.2',s=3,label='Moon',zorder=4)
        axin.scatter(TRSP[1],hourfracs_moonset,color='0.2',s=3,zorder=4)
        fill_RTStimes(hourfracs_moonrise,hourfracs_moonset,facecolor='0.2',alpha=0.1,lw=0.0,zorder=1)
    if showsun==True: 
        axin.scatter(TRSP[1],hourfracs_sunrise,color='#FFFD74',s=3,label='Sun',zorder=4)
        axin.scatter(TRSP[1],hourfracs_sunset,color='#FFFD74',s=3,zorder=4)
        fill_RTStimes(hourfracs_sunrise,hourfracs_sunset,facecolor='#FFFD74',alpha=0.3,lw=0.0,zorder=2)
    
    axin.legend(loc='best')#loc='upper right')
    axin.set_xlim(TRSP[0][0],TRSP[0][-1]); 
    axin.xaxis.set_major_locator(mdates.MonthLocator())
    axin.xaxis.set_major_formatter(ticker.NullFormatter())
    axin.xaxis.set_minor_formatter(xaxisformatter)
    axin.xaxis.set_minor_locator(mdates.MonthLocator(bymonthday=16))
    for tick in axin.xaxis.get_minor_ticks():
        tick.tick1line.set_markersize(0); tick.tick2line.set_markersize(0); tick.label1.set_horizontalalignment('center')
    #fig1.autofmt_xdate(rotation=0,ha='left')
    #axin.set_ylim(dt.datetime.strptime(obsyear+'/01/01 00:00:00','%Y/%m/%d %H:%M:%S'),dt.datetime.strptime(obsyear+'/12/31 23:59:59','%Y/%m/%d %H:%M:%S'))
    axin.set_ylim(0,24)
    axin.set_yticks(np.arange(0,24.1,3)); axin.set_yticklabels(['{0:.0f}:00'.format(h) for h in np.arange(0,24.1,3)]); 
    #axin.yaxis.set_major_formatter(yaxisformatter) #This works, but I don't know how to properly set ylims in hours...
    #axin.yaxis.set_major_locator(mdates.HourLocator(interval=3)) 
    axin.set_xlabel('Month in %i'%(TRSP[0][0].year)); 
    axin.set_ylabel('Time of day [hour, %s]'%(['local time' if local==True else 'UTC'][0]))
    plt.grid(True,color='0.92')#,'both')
    
    if savepath is not '': plt.savefig(savepath,bbox_inches='tight',dpi=dpi)
    if showplot==True: plt.show()
    plt.clf(); plt.close('all')

def plot_year_darktime(target, observer, obsyear, timezone='auto', time_of_obs='night', local=True, showmoon=True, figsize=(14,8), dpi=200, savepath='', showplot=False, xaxisformatter=mdates.DateFormatter('%b'), yaxisformatter=mdates.DateFormatter('%H:%M')):
    """
    Plot the dark time over the course of a year, from a specified observing site, and plot the target altitude track. 
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The target observation source
    observer : ephem.Observer() 
        Observatory/antenna location object to use for the calculations
    obsyear : int, float, str, ephem.Date, or datetime.datetime
        The year to consider for the calculations.  Use 4-digit format 'YYYY' if given as a string.
    timezone : str or None
        Timezone (standard Olson names) for local time on upper x-axis.  Options: \n
            'none' or None : Do not add local time info to the upper axis
            string of a standard Olson database name [e.g., 'US/Mountain' or 'America/Chicago'] :  use this directly for dt calculations \n
            'auto' or 'calculate' : compute the timezone from the observer lat/lon using module tzwhere.  Takes a while (but is convenient!) 
    time_of_obs : str  ['noon', 'night', 'midnight', 'peak' or 'HH:MM:SS'] 
        Denotes how to generate the times to use for each day. Options: \n
            'noon' = calculate the values at 12:00 local time each day \n
            'night' or 'midnight' = use 23:59:59 local time each day (for night obs.  00:00:00 would be for the previous night.) \n
            'peak' = return values during daily peak altitude \n
    local : bool 
        True returns times in the local observer site time. False returns in UTC time.
    showmoon : bool
        Set to True to include a track on the plot for the Moon
    figsize : tuple
        matplotlib Figure size
    dpi : int
        matplotlib image resolution (dpi = dots per inch)
    savepath : str
        The path for saving the figure, including the filename and extension. 
    showplot : bool
        Set to True to invoke plt.show() and display the interactive plot
    xaxisformatter: mdates.DateFormatter
        Formatter object for the xaxis dates.  Default mdates.DateFormatter('%b') will print abbreviated month names, '%B' will print full month names.
    yaxisformatter: mdates.DateFormatter
        Formatter object for the yaxis times.  Default mdates.DateFormatter('%H:%M') will print the hour and minute.
    
    Example
    -------
    wht = ephem.Observer(); wht.name='WHT' \n
    wht.lon='-17 52 53.8'; wht.lat='28 45 37.7'; wht.elevation=2344; \n
    crab=obs.create_ephem_target('Crab Nebula','05:34:31.94','22:00:52.2') #'M1' \n
    obs.plot_year_darktime(crab, wht, '2021', showplot=True)
    """
    #if True in [i in timezone.lower() for i in ['auto','calc']]:
    #    #Use tzwhere to compute the timezone based on the observer lat/lon (input in degrees)
    #    try: timezone=tzwhere.tzwhere().tzNameAt(observer.lat*180/np.pi, wrap_pm180(observer.lon*180/np.pi))
    #    except: timezone=tzwhere.tzwhere(forceTZ=True).tzNameAt(observer.lat*180/np.pi, wrap_pm180(observer.lon*180/np.pi), forceTZ=True)
    timezone = tz_from_observer(observer)
    
    #Calculate the daily transit times, rise times, set times, and peak altitudes for the specified year
    TRSP=compute_yearly_target_data(target, observer, obsyear, timezone=timezone, time_of_obs=time_of_obs, peak_alt=True, local=local)
    
    if showmoon==True: 
        moontimes=np.array([calculate_moon_times(observer,t,outtype='dt') for t in TRSP[0]]) #dt format, [Ndays,2]
        nomoonhours=np.array([wrap_24hr(t.total_seconds()/3600.) for t in moontimes[:,1]-moontimes[:,0]])
    
    suntimes=np.array([calculate_twilight_times(observer,t.replace(hour=12)) for t in TRSP[0]])[:,0].astype(ephem.Date) #ephem.Date format, [sunset,sunrise]
    for i in range(suntimes.shape[0]): 
        suntimes[i,0]=ephem.Date(suntimes[i,0]).datetime(); suntimes[i,1]=ephem.Date(suntimes[i,1]).datetime(); #in UTC still
    nosunhours=np.array([wrap_24hr(t.total_seconds()/3600.) for t in suntimes[:,1]-suntimes[:,0]]) 
    
    #Calculate the number of dark hours for the target each day (considering rise/set time, but not Moon)
    midnights = create_obstime_array( dt.datetime(int(obsyear),1,1,23,59,59), dt.datetime(int(obsyear),12,31,23,59,59), timezone_string=timezone, output_as_utc=False, n_steps=365)
    
    targetdarkhours=np.zeros(len(TRSP[0]))
    for i in list(range(len(TRSP[0]))):
        targetdarkhours[i]=calculate_target_darktime_singlenight(target, observer, midnights[i])
    
    fig1=plt.figure(1,figsize=figsize)
    axin=fig1.add_subplot(111)#,adjustable='box-forced',aspect=5e-3)
    
    axin.plot(TRSP[1],targetdarkhours,color='#A83C09',zorder=5,label='Target dark time')
    
    axin.plot(TRSP[0],nosunhours,color='#FFFD74',zorder=4,label='Sun')
    axin.fill_between(TRSP[1],nosunhours,y2=np.zeros(len(nosunhours)),facecolor='0.1',zorder=1)
    axin.fill_between(TRSP[1],nosunhours,y2=np.zeros(len(nosunhours))+24.,facecolor='#FFFD74',alpha=0.3)
    
    if showmoon==True:
        #for each time step: fill_between with alpha based on fraction of night that moon is up   (***and brightness??)
        for i in list(range(len(moontimes)-1)):
            axin.fill_between(TRSP[1][i:i+2],nosunhours[i:i+2],y2=[0,0],facecolor='#D9D5C5', alpha=(24.-nomoonhours[i])/24., lw=0., zorder=3)
    
    axin.legend(loc='best')#loc='upper right')
    axin.set_xlim(TRSP[0][0],TRSP[0][-1]); 
    axin.xaxis.set_major_locator(mdates.MonthLocator())
    axin.xaxis.set_major_formatter(ticker.NullFormatter())
    axin.xaxis.set_minor_formatter(xaxisformatter)
    axin.xaxis.set_minor_locator(mdates.MonthLocator(bymonthday=16))
    for tick in axin.xaxis.get_minor_ticks():
        tick.tick1line.set_markersize(0); tick.tick2line.set_markersize(0); tick.label1.set_horizontalalignment('center')
    axin.set_ylim(0,24)
    axin.set_yticks(np.arange(0,24.1,3)); axin.set_yticklabels(['{0:.0f}:00'.format(h) for h in np.arange(0,24.1,3)]); 
    axin.set_xlabel('Month in %i'%(TRSP[0][0].year)); 
    axin.set_ylabel('Hours of Darkness')
    plt.grid(True,color='0.92')#,'both')
    
    if savepath is not '': plt.savefig(savepath,bbox_inches='tight',dpi=dpi)
    if showplot==True: plt.show()
    plt.clf(); plt.close('all')

    

def fill_twilights(axin,obsframe,startdate,offsetdatetime=0.,timetype='offset',bgcolor='k'):
    """
    Fills a matplotlib axis with shading for twilight.
    
    Parameters
    ----------
    axin : matplotlib.axes._subplots.AxesSubplot
        The matplotlib axis object to plot to, created with plt.subplot(), plt.axis() etc.
    obsframe : ephem.Observer()
        pyephem Observer() frame
    startdate : ephem.Date or str, formatted as e.g., 'YYYY/MM/DD hh:mm:ss'
        The starting date/time (the evening at the beginning observations, for classical night sessions) 
    offsetdatetime : datetime.datetime
        Reference time to use for plotting offset times
    timetype : str, 'offset' or 'abs'
        String denoting the style to use for displaying the time: either 'abs' for absolute time, or 'offset' for offset/delta time from speccified offsetdatetime.
    bgcolor : str, or valid matplotlib color
        Base color to use for plot background fill shading.  
    """
    sunset,t_civil,t_naut,t_astro=calculate_twilight_times(obsframe,startdate)
    if timetype=='abs':
        s_m=[convert_ephem_datetime(ephem.Date(sunset[0])),convert_ephem_datetime(ephem.Date(sunset[1]))]
        tc_m=[convert_ephem_datetime(ephem.Date(t_civil[0])),convert_ephem_datetime(ephem.Date(t_civil[1]))]
        tn_m=[convert_ephem_datetime(ephem.Date(t_naut[0])),convert_ephem_datetime(ephem.Date(t_naut[1]))]
        ta_m=[convert_ephem_datetime(ephem.Date(t_astro[0])),convert_ephem_datetime(ephem.Date(t_astro[1]))]
    else: 
        s_m=(sunset-offsetdatetime)/ephem.hour; tc_m=(t_civil-offsetdatetime)/ephem.hour;
        tn_m=(t_naut-offsetdatetime)/ephem.hour; ta_m=(t_astro-offsetdatetime)/ephem.hour;
    axin.axvspan(*s_m,color=bgcolor,alpha=0.1,zorder=-10,ec=None)
    axin.axvspan(*tc_m,color=bgcolor,alpha=0.3,zorder=-9,ec=None)
    axin.axvspan(*tn_m,color=bgcolor,alpha=0.5,zorder=-8,ec=None)
    axin.axvspan(*ta_m,color=bgcolor,alpha=0.7,zorder=-7,ec=None)
    #axin.annotate('Sunset %s'%(s_m[0]), xy=(s_m[0],80), xytext=(s_m[0]-10./60./24,80), arrowprops=None)
    if timetype=='abs':
        #axin.text(convert_ephem_datetime(ephem.Date(sunset[0]-10./60./24)),80,'Sunset %s'%('{}:{}'.format(*deg2dms(wrap_24hr(24+s_m[0])))),color='k',fontsize=5)
        axin.text(convert_ephem_datetime(ephem.Date(sunset[0]-7./60./24)),80.5,'Sunset %s'%(s_m[0].time().strftime('%H:%M')),color='k',fontsize=6,rotation=90,va='bottom') 
        axin.text(convert_ephem_datetime(ephem.Date(t_astro[0]-7./60./24)),80.5,'Ast.Twi. %s'%(ta_m[0].time().strftime('%H:%M')),color='w',fontsize=6,rotation=90,va='bottom') 
        axin.text(convert_ephem_datetime(ephem.Date(t_astro[1]-7./60./24)),80.5,'Ast.Twi. %s'%(ta_m[1].time().strftime('%H:%M')),color='w',fontsize=6,rotation=90,va='bottom') 
        axin.text(convert_ephem_datetime(ephem.Date(sunset[1]-7./60./24)),80.5,'Sunrise %s'%(s_m[1].time().strftime('%H:%M')),color='k',fontsize=6,rotation=90,va='bottom') 

def fill_twilights_light(axin, obsframe, startdate, plottimerange, offsetdatetime=0., bgcolor='#95D0FC', timetype='offset'):
    """
    Fills a matplotlib axis with shading for twilight -- alternate variant for light colors.
    
    Parameters
    ----------
    axin : matplotlib.axes._subplots.AxesSubplot
        The matplotlib axis object to plot to, created with plt.subplot(), plt.axis() etc.
    obsframe : ephem.Observer()
        pyephem Observer() frame
    startdate : ephem.Date or str, formatted as e.g., 'YYYY/MM/DD hh:mm:ss'
        The starting date/time (the evening at the beginning observations, for classical night sessions) 
    plottimerange : array-like [list, tuple, ...]
        The beginning and end times (dt.datetime format) use for the plot x-range span
    offsetdatetime : datetime.datetime
        Reference time to use for plotting offset times
    timetype : str, 'offset' or 'abs'
        String denoting the style to use for displaying the time: either 'abs' for absolute time, or 'offset' for offset/delta time from speccified offsetdatetime.
    bgcolor : str, or valid matplotlib color
        Base color to use for plot background fill shading.  
    """
    midtime=(ephem.Date(ephem.Date(plottimerange[0]))+ephem.Date(plottimerange[-1]))/2.
    sunset,t_civil,t_naut,t_astro=calculate_twilight_times(obsframe,midtime)
    
    if timetype=='abs':
        s_m=[convert_ephem_datetime(ephem.Date(sunset[0])),convert_ephem_datetime(ephem.Date(sunset[1]))]
        tc_m=[convert_ephem_datetime(ephem.Date(t_civil[0])),convert_ephem_datetime(ephem.Date(t_civil[1]))]
        tn_m=[convert_ephem_datetime(ephem.Date(t_naut[0])),convert_ephem_datetime(ephem.Date(t_naut[1]))]
        ta_m=[convert_ephem_datetime(ephem.Date(t_astro[0])),convert_ephem_datetime(ephem.Date(t_astro[1]))]
    else: 
        s_m=(sunset-offsetdatetime)/ephem.hour; tc_m=(t_civil-offsetdatetime)/ephem.hour;
        tn_m=(t_naut-offsetdatetime)/ephem.hour; ta_m=(t_astro-offsetdatetime)/ephem.hour;
    axin.axvspan(plottimerange[0],s_m[0],color=bgcolor,zorder=-10,alpha=1.0,ec=None)
    axin.axvspan(s_m[0],tc_m[0],color=bgcolor,zorder=-10,alpha=.7,ec=None)
    axin.axvspan(tc_m[0],tn_m[0],color=bgcolor,zorder=-10,alpha=.5,ec=None)
    axin.axvspan(tn_m[0],ta_m[0],color=bgcolor,zorder=-10,alpha=.2,ec=None)
    axin.axvspan(ta_m[1],tn_m[1],color=bgcolor,zorder=-10,alpha=.2,ec=None)
    axin.axvspan(tn_m[1],tc_m[1],color=bgcolor,zorder=-10,alpha=.5,ec=None)
    axin.axvspan(tc_m[1],s_m[1],color=bgcolor,zorder=-10,alpha=.7,ec=None)
    axin.axvspan(s_m[1],plottimerange[-1],color=bgcolor,zorder=-10,alpha=1.0,ec=None)
    if timetype=='abs':
        #axin.text(convert_ephem_datetime(ephem.Date(sunset[0]-10./60./24)),80,'Sunset %s'%('{}:{}'.format(*deg2dms(wrap_24hr(24+s_m[0])))),color='k',fontsize=5)
        if ephem.Date(s_m[0])>ephem.Date(plottimerange[0]):
            axin.text(convert_ephem_datetime(ephem.Date(sunset[0]-7./60./24)),80.5,'Sunset %s'%(s_m[0].time().strftime('%H:%M')),color='k',fontsize=6,rotation=90,va='bottom') 
        if ephem.Date(t_astro[0])>ephem.Date(plottimerange[0]):
            axin.text(convert_ephem_datetime(ephem.Date(t_astro[0]-7./60./24)),80.5,'Ast.Twi. %s'%(ta_m[0].time().strftime('%H:%M')),color='k',fontsize=6,rotation=90,va='bottom') 
        if ephem.Date(t_astro[1])<ephem.Date(plottimerange[-1]):
            axin.text(convert_ephem_datetime(ephem.Date(t_astro[1]-7./60./24)),80.5,'Ast.Twi. %s'%(ta_m[1].time().strftime('%H:%M')),color='k',fontsize=6,rotation=90,va='bottom') 
        if ephem.Date(s_m[1])<ephem.Date(plottimerange[-1]):
            axin.text(convert_ephem_datetime(ephem.Date(sunset[1]-7./60./24)),80.5,'Sunrise %s'%(s_m[1].time().strftime('%H:%M')),color='k',fontsize=6,rotation=90,va='bottom') 

def fill_twilights_split(axin,obsframe,startdate,enddate, bgcolor='k'):
    """
    Fill a matplotlib axis with shading for twilight; for daytime observations, so centered in the daylight hours with twilight shading split into two parts.
    
    Parameters
    ----------
    axin : matplotlib.axes._subplots.AxesSubplot
        The matplotlib axis object to plot to, created with plt.subplot(), plt.axis() etc.
    obsframe : ephem.Observer()
        pyephem Observer() frame
    startdate : ephem.Date or str, formatted as e.g., 'YYYY/MM/DD hh:mm:ss'
        The starting date/time (the morning at the beginning of observations, for daytime sessions) 
    enddate : ephem.Date or str, formatted as e.g., 'YYYY/MM/DD hh:mm:ss'
        The ending date/time (the evening at the end of observations, for daytime sessions) 
    bgcolor : str, or valid matplotlib color
        Base color to use for plot background fill shading.  
    
    Notes
    -----
    Civil, Nautical, and Astronomical twilight defined as when Sun is > [6,12,18] deg below horizon.
    """
    start_time_dt=dt.datetime.strptime(startdate,'%Y/%m/%d %H:%M:%S')
    end_time_dt=dt.datetime.strptime(enddate,'%Y/%m/%d %H:%M:%S')
    start_time=ephem.Date(startdate); end_time=ephem.Date(enddate)
    
    tmpobsframe=obsframe.copy()
    tmpobsframe.date=startdate
    
    #Determine if start date/time is in night or day
    if tmpobsframe.next_rising(ephem.Sun(), use_center=True) < tmpobsframe.next_setting(ephem.Sun(), use_center=True):
        #In this case, the start date/time occurs before sunrise
        sunrise=tmpobsframe.next_rising(ephem.Sun(),use_center=True) 
        sunset=tmpobsframe.next_setting(ephem.Sun(),use_center=True)
        ##Civil Twilight
        #tmpobsframe.horizon='+6'
        #t_civil1=tmpobsframe.next_setting(ephem.Sun(),use_center=True), tmpobsframe.next_rising(ephem.Sun(),use_center=True)]
        #
        ##Nautical Twilight
        #tmpobsframe.horizon='+12'
        #t_naut=[tmpobsframe.previous_setting(ephem.Sun(),use_center=True), tmpobsframe.next_rising(ephem.Sun(),use_center=True)]
        #
        ##Astronomical Twilight
        #tmpobsframe.horizon='+18'
        #t_astro=[tmpobsframe.previous_setting(ephem.Sun(),use_center=True), tmpobsframe.next_rising(ephem.Sun(),use_center=True)]
        # 
        #sunset1,t_civil1,t_naut1,t_astro1=calculate_twilight_times(obsframe,startdate)
        #sunset2,t_civil2,t_naut2,t_astro2=calculate_twilight_times(obsframe,enddate)
        
        #### Eventually need to figure out more clever way of doing this...
        
    else:
        #In this case, the start date/time occurs between sunrise and sunset.  Proceed as in standard nighttime case. 
        sunset,t_civil,t_naut,t_astro=calculate_twilight_times(obsframe,startdate)
        s_m=[convert_ephem_datetime(ephem.Date(sunset[0])),convert_ephem_datetime(ephem.Date(sunset[1]))]
        tc_m=[convert_ephem_datetime(ephem.Date(t_civil[0])),convert_ephem_datetime(ephem.Date(t_civil[1]))]
        tn_m=[convert_ephem_datetime(ephem.Date(t_naut[0])),convert_ephem_datetime(ephem.Date(t_naut[1]))]
        ta_m=[convert_ephem_datetime(ephem.Date(t_astro[0])),convert_ephem_datetime(ephem.Date(t_astro[1]))]
        
        axin.axvspan(*s_m,color=bgcolor,alpha=0.1,zorder=-10,ec=None)
        axin.axvspan(*tc_m,color=bgcolor,alpha=0.3,zorder=-9,ec=None)
        axin.axvspan(*tn_m,color=bgcolor,alpha=0.5,zorder=-8,ec=None)
        axin.axvspan(*ta_m,color=bgcolor,alpha=0.7,zorder=-7,ec=None)

        axin.text(convert_ephem_datetime(ephem.Date(sunset[0]-7./60./24)),70.5,'Sunset %s'%(s_m[0].time().strftime('%H:%M')),color=bgcolor,fontsize=6,rotation=90,va='bottom') 
        axin.text(convert_ephem_datetime(ephem.Date(t_astro[0]-7./60./24)),70.5,'Ast.Twi. %s'%(ta_m[0].time().strftime('%H:%M')),color=bgcolor,fontsize=6,rotation=90,va='bottom') 
        axin.text(convert_ephem_datetime(ephem.Date(t_astro[1]-7./60./24)),80.5,'Ast.Twi. %s'%(ta_m[1].time().strftime('%H:%M')),color=bgcolor,fontsize=6,rotation=90,va='bottom') 
        axin.text(convert_ephem_datetime(ephem.Date(sunset[1]-7./60./24)),80.5,'Sunrise %s'%(s_m[1].time().strftime('%H:%M')),color=bgcolor,fontsize=6,rotation=90,va='bottom') 


def plot_night_observing_tracks(target_list, observer, obsstart, obsend, weights=None, mode='nearest', plotmeantransit=False, toptime='local', timezone='auto', n_steps=1000, simpletracks=False, azcmap='rainbow',light_fill=False, bgcolor='k', xaxisformatter=mdates.DateFormatter('%H:%M'), figsize=(14,8), dpi=200, savepath='',showplot=False):
    """
    Plot visibility tracks of an astronomical target on the sky, over the course of a night - like a classic 'staralt' plot. 
    
    Parameters
    ----------
    target_list : ephem.FixedBody(), or array-like [list, tuple, np.array...] 
        List of pyephem FixedBody sky sources.  Can also input a single target source.
    observer : ephem.Observer() 
        Observatory/antenna/location object to use for the calculations
    obsstart : emphem.Date() or str [formatted as 'YYYY/MM/DD HH:MM:SS.s']
        Observation start time
    obsend : emphem.Date() or str [formatted as 'YYYY/MM/DD HH:MM:SS.s']
        Observation end time
    weights : np.array, list, or other array-like
        Optional (float) weights to give to certain stations.  (e.g. to prioritize certain baselines)   Larger values give higher importance. Passed to calculate_targets_mean_transit_time()
    mode : str
        Specifies which transit time to return. Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    plotmeantransit : bool 
        Set to True to plot a vertical line denoting the mean of the target transit times
    toptime : str ['local', 'LMST', 'none'] 
        Type of time to plot on the top axis: local time, LMST, or nothing
    timezone : str or None
        Timezone (standard Olson names) for local time on upper x-axis.  Options: \n
            'none' or None : Do not add local time info to the upper axis
            string of a standard Olson database name [e.g., 'US/Mountain' or 'America/Chicago'] :  use this directly for dt calculations \n
            'auto' or 'calculate' : compute the timezone from the observer lat/lon using module tzwhere.  Takes a while (but is convenient!)
    n_steps : int 
        The number of positions to calculate in the altitude tracks
    simpletracks : bool 
        Whether to plot the visibility tracks as simple lines (True) or with a colormap (False). \n
        colormap option only works for a single target...
    azcmap : str, or matplotlib.colors.LinearSegmentedColormap
        The colormap to use for azimuth, passed to plt.scatter().  Any valid matplotlib colormap should work (including strings for builtin cmaps such as 'inferno').
    light_fill : bool, or 'none'
        Whether to fill the background with light (True) color, dark (False, default) color, or nothing ('none')
    bgcolor : str, or valid matplotlib color
        Base color to use for plot background fill shading.  
    xaxisformatter: mdates.DateFormatter
        Formatter object for the xaxis dates.  Default mdates.DateFormatter('%b') will print abbreviated month names, '%B' will print full month names.    
    figsize : tuple
        matplotlib Figure size
    dpi : int
        matplotlib image resolution (dpi = dots per inch)
    savepath : str
        The path for saving the figure, including the filename and extension. 
    showplot : bool
        Set to True to invoke plt.show() and display the interactive plot
    
    Example
    -------
    ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751') \n
    mrk348=obs.create_ephem_target('Mrk348','00:48:47.14','+31:57:25.1') \n
    obsstart=obs.dtaware_to_ephem(obs.construct_datetime('2021/09/15 16:00:00','dt',timezone='US/Pacific')) \n
    obsend=obs.dtaware_to_ephem(obs.construct_datetime('2021/09/16 10:00:00','dt',timezone='US/Pacific')) \n
    obs.plot_night_observing_tracks([ngc1052,mrk348],obs.vlbaBR,obsstart,obsend, plotmeantransit=False, simpletracks=True, toptime='local', timezone='calculate', n_steps=1000, azcmap='rainbow', light_fill=True, savepath='test.jpg', showplot=False)
    """
    if obsend<obsstart: raise Exception('plot_night_observing_tracks(): obsend is earlier than obsstart!')
    #if type(obsstart) is str: obsstart=ephem.Date(obsstart)
    #if type(obsend) is str: obsend=ephem.Date(obsend)
    obsstart=ephem.Date(obsstart)
    obsend=ephem.Date(obsend)
    #If the input target_list is just a single object, put it in a list, so it plays nice with the script.
    #   --> checking if np.ndim(target_list)==0 to check if it's array_like, as only lists etc. will have ndim>0.
    #       Per the official docs, using np.ndim(object)==0 is preferable to np.isscalar(object).
    if np.ndim(target_list)==0: target_list=[target_list,] #Putting the single input target into a list
    
    fig1=plt.figure(1,figsize=figsize)
    axin=fig1.add_subplot(111)#,adjustable='box-forced',aspect=5e-3)
    
    times_utc=create_obstime_array(obsstart.datetime().strftime('%Y/%m/%d %H:%M:%S'), obsend.datetime().strftime('%Y/%m/%d %H:%M:%S'), timezone_string='UTC', n_steps=n_steps)
    #mid_night=ephem.Date(dt.datetime.strptime( obsstart.datetime().strftime('%Y/%m/%d')[:10]+' 23:59:59','%Y/%m/%d %H:%M:%S' ) )
    meanobstime=np.mean([obsstart,obsend])
    obsday=obsstart.datetime().strftime('%Y/%m/%d')[:10]
    
    sun_alts=compute_sun_tracks(observer,obsstart,obsend,nsteps=n_steps)[0]
    moon_alts=compute_moon_tracks(observer,obsstart,obsend,nsteps=n_steps)[0]
    moonphase_perc=compute_moonphase(obsstart,return_fmt='percent')
    moontimes=calculate_moon_times(observer,obsstart,outtype='dt') #[moonrise,moonset]
    
    if light_fill is True: 
        #R,S=calculate_rise_set_times_single(target,observer,obsstart,mode='nearest',return_fmt='ephem')
        fill_twilights_light(axin, observer, meanobstime, times_utc, bgcolor='#FEFFCA', timetype='abs')
    elif light_fill is False: fill_twilights(axin, observer, meanobstime, timetype='abs', bgcolor=bgcolor)
    else: pass #Don't do any background fill
    
    transit_times=np.zeros(len(target_list)).astype(str)
    for target,i in zip(target_list,list(range(len(target_list)))):
        alts,azs=compute_target_altaz(target,observer,obsstart,obsend,nsteps=n_steps)
        transit_times[i]=calculate_transit_time_single(target,observer,meanobstime, return_fmt='str', mode=mode) 
        #moonseps=moonsep_timearray(target,observer,times_utc)
        #sunseps=sunsep_timearray(target,observer,times_utc)
        moonsep_transit=moonsep_single(target,observer,transit_times[i])
        sunsep_transit=sunsep_single(target,observer,transit_times[i])
        if target.name is '': target.name='Target_%i'%(i)
        if simpletracks==True: axin.plot(times_utc,alts,lw=1.5,label=target.name+'\n (Moon sep. = %i deg)\n (Sun sep. = %i deg)'%(moonsep_transit,sunsep_transit), zorder=5)
        else: colortrack=axin.scatter(times_utc,alts,c=azs,label=target.name+'\n (Moon sep. = %i deg)'%(moonsep_transit), lw=0, s=10, cmap=azcmap, zorder=5)
        
        axin.axvline(dt.datetime.strptime(transit_times[i],'%Y/%m/%d %H:%M:%S'), ls=':',color='0.5',zorder=-1)
        axin.annotate('Transit time of %s =  %s (UTC)'%(target.name,transit_times[i][-8:]),xy=[dt.datetime.strptime(transit_times[i],'%Y/%m/%d %H:%M:%S'),50.],xytext=[dt.datetime.strptime(transit_times[i],'%Y/%m/%d %H:%M:%S')+dt.timedelta(minutes=-25),50.], va='center', rotation=90, color='0.5')
    
    if simpletracks is False: fig1.colorbar(colortrack,ax=axin,pad=.12).set_label('Azimuth [deg]')
    
    axin.plot(times_utc,moon_alts,color='0.6',lw=0.8,ls='--',label='Moon (%i%%)'%(moonphase_perc),alpha=0.5,zorder=1)
    axin.plot(times_utc,sun_alts,color='r',lw=0.8,ls='--',label='Sun',alpha=0.5,zorder=1)
    #fill_twilights_split(axin,vlbaBR,'2021/04/15 00:00:00','2021/04/15 23:59:00')
    axin.legend(loc='best')#loc='upper right')
    
    axin.set_ylim(0,90)
    axin.xaxis.set_major_formatter(xaxisformatter)#hourformat)
    axin.set_xlabel('UTC Time, for Local Starting Night %s'%(str(obsstart.datetime().date()))); plt.ylabel('Altitude [deg]')
    plt.grid(True,color='0.92')#,'both')
    
    if plotmeantransit==True: 
        meantransit=calculate_targets_mean_transit_time(target_list,observer,meanobstime,weights=weights)
        axin.axvline(dt.datetime.strptime(meantransit,'%Y/%m/%d %H:%M:%S'), ls='--',color='0.5',zorder=-1)
        axin.annotate('%sMean of target transit times =  %s (UTC)'%(['' if weights is None else 'Weighted '][0],meantransit[-8:]),xy=[dt.datetime.strptime(meantransit,'%Y/%m/%d %H:%M:%S'),85.],xytext=[dt.datetime.strptime(meantransit,'%Y/%m/%d %H:%M:%S')+dt.timedelta(minutes=-25),85.], va='top', rotation=90, color='0.5')
    
    axin.set_xlim(times_utc[0],times_utc[-1]); 
    
    axin2=axin.twinx(); #axin2.set_ylim(0,90)
    axin2.set_ylim(axin.get_ylim()) #--> This actually returns slightly different lims!!  Not [0,90]
    axin2.set_yticks([10,30,50,70,90]); #Just compute and set them manually on the ALTITUDE frame
    axin2.set_yticklabels(np.round(alt2airmass(np.array([10,30,50,70,90])),decimals=2)) 
    axin2.set_ylabel('Airmass [sec(z)]',rotation=-90,va='bottom')
    
    axin3=axin.twiny(); 
    if 'loc' in toptime.lower():
        #if True in [i in timezone.lower() for i in ['auto','calc']]:
        #    #Use tzwhere to compute the timezone based on the observer lat/lon (input in degrees)
        #    try: timezone=tzwhere.tzwhere().tzNameAt(observer.lat*180/np.pi, wrap_pm180(observer.lon*180/np.pi))
        #    except: timezone=tzwhere.tzwhere(forceTZ=True).tzNameAt(observer.lat*180/np.pi, wrap_pm180(observer.lon*180/np.pi), forceTZ=True)
        timezone = tz_from_observer(observer)
        utcoffset=calculate_dtnaive_utcoffset(ephem.Date(meanobstime).datetime(),local_timezone=timezone)
        times_local=np.zeros_like(times_utc)
        for i in range(len(times_utc)): times_local[i]=times_utc[i].astimezone(pytz.timezone(timezone))
        axin3.set_xlim(times_local[0],times_local[-1],auto=True) #Local time... still need to set tz set below...
        axin3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M',tz=times_local[0].tzinfo)) #NEED to set tz for axis update!
        axin3.set_xlabel('Local Time $\\rightarrow$', x=-.08, va='top') 
        axin3.text(-.04, 0.95, '(%s)\nUTC offset = %i'%(times_local[0].tzinfo.zone,utcoffset),va='bottom', ha='right', transform=axin3.transAxes, color='k', fontsize=8)
    
    elif toptime.upper()=='LMST':
        ### LMST, calculate each major ticklabel
        axin3.set_xlim(times_utc[0],times_utc[-1]); axin3.xaxis.set_major_formatter(xaxisformatter)#hourformat)
        fig1.canvas.draw()#For some reason, forcing it to draw is necessary to make the ticklabels grabbable.
        utctime0=compute_sidereal_time(observer,times_utc[0],as_type='datetime')
        #utc_ticks=axin3.xaxis.get_major_ticklabels() #Use ax3 instead of ax1, because ax1 ticklabels have already been altered
        #utc_ticks=axin3.xaxis.get_major_ticks() #Use ax3 instead of ax1, because ax1 ticklabels have already been altered
        utc_ticks=axin3.get_xticklabels(); xtl=[]
        for t in utc_ticks:
            t_utc=dt.datetime.strptime(obsday+' '+t.get_text(),'%Y/%m/%d %H:%M').replace(tzinfo=pytz.utc) #need the Y/M/D!
            if t_utc.time()<utctime0.time(): t_utc=t_utc.replace(day=t_utc.day+1)
            #t_lst=LMST_from_UTC(t_utc,observer.lon*180./np.pi)
            t_lst=LST_from_local(t_utc,observer.lon*180./np.pi)
            #t_lst=compute_sidereal_time(observer,t_utc,as_type='dms')
            xtl.append('%i:%i'%(t_lst[0],t_lst[1]))
        axin3.xaxis.set_ticklabels(xtl)
        axin3.set_xlabel('LMST $\\rightarrow$',x=-.05,va='top') 
    #for t in utc_ticks:
    #    t_utc=dt.datetime.strptime(t.get_text(),'%H:%M')
    #    #t_utc=dt.datetime.strptime(t.label2.get_text(),'%H:%M')
    #    t_lst=LMST_from_UTC(t_utc,observer.lon*180./np.pi)
    #    #t.label2.set_text('%i:%i'%(t_lst[0],t_lst[1]))
    #    t.set_text(t_lst.strftime('%H:%M'))
    #utc_ticks=axin.xaxis.get_major_ticks(); lst_ticks=axin3.xaxis.get_major_ticks();
    #for t1,t3 in zip(utc_ticks,lst_ticks):
    #    t_utc=dt.datetime.strptime(t1.label2.get_text(),'%H:%M')
    #    t_lst=LMST_from_UTC(t_utc,observer.lon*180./np.pi)
    #    t3.label2.set_text('%i:%i'%(t_lst[0],t_lst[1]))
    #xtl = ['(%s)'%l.get_text() for l in axin.xaxis.get_ticklabels()]
    #for i in range(len(xtl)): xtl[i]=LMST_from_UTC(dt.datetime.strptime(xtl[i],'(%H:%M)'),observer.lon*180./np.pi) 
    #axin3.xaxis.set_ticklabels(xtl)
    
    axin3.text(-.05, 0.8, 'Moonrise %s\nMoonset %s'%(moontimes[0].strftime('%H:%M'),moontimes[1].strftime('%H:%M')),va='bottom', ha='right', transform=axin.transAxes, color='k', fontsize=10)

    #plt.suptitle('Moon: Separation = %.2f deg, Illumination = %.2f %%'%(moonsep,moon.moon_phase*100),y=.999)
    
    #plt.subplots_adjust(left=0.2,right=0.98)
    plt.savefig(savepath,bbox_inches='tight',dpi=dpi)
    if showplot==True: plt.show()
    plt.clf(); plt.close('all')


def plot_visibility_tracks_toaxis(target_list, observer, obsstart, obsend, axin, weights=None, mode='nearest', duration_hours=0, plotmeantransit=False, timezone='auto', xaxisformatter=mdates.DateFormatter('%H:%M')):
    """
    Plot visibility tracks for the given targets, for a given ephem observer (telescope location) and obstimes. Includes the mean of the transit times from each target.  Plots onto a specified existing matplotlib/pyplot axis -- useful for creating the plot & allowing further manual editing.  The Moon separation and Sun separation at the target transit times are included in the legend entries.
    
    Parameters
    ----------
    target_list : ephem.FixedBody(), or array-like [list, tuple, np.array...] 
        List of pyephem FixedBody sky sources.  Can also input a single target source.
    observer : ephem.Observer() 
        Observatory/antenna/location object to use for the calculations
    obsstart : emphem.Date() or str [formatted as 'YYYY/MM/DD HH:MM:SS.s']
        Observation start time
    obsend : emphem.Date() or str [formatted as 'YYYY/MM/DD HH:MM:SS.s']
        Observation end time
    axin : matplotlib/pyplot axis object
        The axis on which to plot
    weights : np.array, list, or other array-like
        Optional (float) weights to give to certain stations.  (e.g. to prioritize certain baselines)   Larger values give higher importance. Passed to calculate_targets_mean_transit_time()
    mode : str
        Specifies which transit time to return. Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    duration_hours : float
        The desired duration of observations, for calculating the optimal start time.  If greater than zero, calculate start time and plot a vertical line denoting it. 
    plotmeantransit : bool 
        Set to True to plot a vertical line denoting the mean of the target transit times
    timezone : str or None
        Timezone (standard Olson names) for local time on upper x-axis.  Options: \n
            'none' or None : Do not add local time info to the upper axis
            string of a standard Olson database name [e.g., 'US/Mountain' or 'America/Chicago'] :  use this directly for dt calculations \n
            'auto' or 'calculate' : compute the timezone from the observer lat/lon using module tzwhere.  Takes a while (but is convenient!)
    xaxisformatter: mdates.DateFormatter
        Formatter object for the xaxis dates.  Default mdates.DateFormatter('%b') will print abbreviated month names, '%B' will print full month names.    
    
    Example
    -------
    fig1=plt.figure(1,figsize=(14,8)) \n
    ax1=fig1.add_subplot(111) \n
    ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751') \n
    mrk348=obs.create_ephem_target('Mrk348','00:48:47.14','+31:57:25.1') \n
    obs.plot_visibility_tracks_toaxis([ngc1052,mrk348],vlbaBR, ephem.Date('2021/04/15 00:00:00'),ephem.Date('2021/04/15 23:59:59'), ax1, timezone='US/Pacific') \n
    plt.savefig('ngc1052_fullVLBA_april15.jpg'); plt.clf(); plt.close('all')
    """
    #If the input target_list is just a single object, put it in a list, so it plays nice with the script.
    #   --> checking if np.ndim(target_list)==0 to check if it's array_like, as only lists etc. will have ndim>0.
    #       Per the official docs, using np.ndim(object)==0 is preferable to np.isscalar(object).
    if np.ndim(target_list)==0: target_list=[target_list,] #Putting the single input target into a list
    
    if type(obsstart) is str: obsstart=ephem.Date(obsstart)
    if type(obsend) is str: obsend=ephem.Date(obsend)
    meanobstime=np.mean([obsstart,obsend])
    times_utc=create_obstime_array(obsstart.datetime().strftime('%Y/%m/%d %H:%M:%S'), obsend.datetime().strftime('%Y/%m/%d %H:%M:%S'), timezone_string='UTC', n_steps=100)
    
    transit_times=np.zeros(len(target_list)).astype(str)
    for target,i in zip(target_list,list(range(len(target_list)))):
        alts,azs=compute_target_altaz(target,observer,obsstart,obsend,nsteps=100)
        transit_times[i]=calculate_transit_time_single(target,observer,meanobstime, return_fmt='str', mode=mode)
        #moonseps=moonsep_timearray(target,observer,times_utc)
        #sunseps=sunsep_timearray(target,observer,times_utc)
        moonsep_transit=moonsep_single(target,observer,transit_times[i])
        sunsep_transit=sunsep_single(target,observer,transit_times[i])
        if target.name is '': target.name='Target_%i'%(i)
        axin.plot(times_utc,alts,lw=1.5,label=target.name+'\n (Moon sep. = %i deg)\n (Sun sep. = %i deg)'%(moonsep_transit,sunsep_transit), zorder=5)
        
        axin.axvline(dt.datetime.strptime(transit_times[i],'%Y/%m/%d %H:%M:%S'), ls=':',color='0.5',zorder=-1)
        axin.annotate('Transit time of %s =  %s (UTC)'%(target.name,transit_times[i][-8:]),xy=[dt.datetime.strptime(transit_times[i],'%Y/%m/%d %H:%M:%S'),50.],xytext=[dt.datetime.strptime(transit_times[i],'%Y/%m/%d %H:%M:%S')+dt.timedelta(minutes=-25),50.], va='center', rotation=90,color='0.5')
    
    sun_alts=compute_sun_tracks(observer,obsstart,obsend,nsteps=100)[0]
    moon_alts=compute_moon_tracks(observer,obsstart,obsend,nsteps=100)[0]
    moonphase_perc=compute_moonphase(obsstart,return_fmt='percent')
    
    axin.plot(times_utc,moon_alts,color='0.6',lw=0.8,ls='--',label='Moon (%i%%)'%(moonphase_perc),alpha=0.5,zorder=1)
    axin.plot(times_utc,sun_alts,color='r',lw=0.8,ls='--',label='Sun',alpha=0.5,zorder=1)
    #fill_twilights_split(axin,vlbaBR,'2021/04/15 00:00:00','2021/04/15 23:59:00')
    axin.legend(loc='best')#loc='upper right')
    #plt.colorbar(pad=.12).set_label('Azimuth [deg]')
    axin.set_xlim(times_utc[0],times_utc[-1]); 
    axin.set_ylim(0,90)
    axin.xaxis.set_major_formatter(xaxisformatter)
    axin.set_xlabel('UTC Time, Starting %s'%(str(obsstart.datetime().date()))); plt.ylabel('Altitude [deg]')
    plt.grid(True,color='0.92')#,'both')
    
    if plotmeantransit==True: 
        meantransit=calculate_targets_mean_transit_time(target_list,observer,meanobstime,weights=weights)
        axin.axvline(dt.datetime.strptime(meantransit,'%Y/%m/%d %H:%M:%S'), ls='--',color='0.5',zorder=-1)
        axin.annotate('%sMean of target transit times =  %s (UTC)'%(['' if weights is None else 'Weighted '][0],meantransit[-8:]),xy=[dt.datetime.strptime(meantransit,'%Y/%m/%d %H:%M:%S'),85.],xytext=[dt.datetime.strptime(meantransit,'%Y/%m/%d %H:%M:%S')+dt.timedelta(minutes=-25),85.], va='top', rotation=90, color='0.5')
    
    axin2=axin.twinx(); #axin2.set_ylim(0,90)
    axin2.set_ylim(axin.get_ylim()) #--> This actually returns slightly different lims!!  Not [0,90]
    axin2.set_yticks([10,30,50,70,90]); #Just compute and set them manually on the ALTITUDE frame
    axin2.set_yticklabels(np.round(alt2airmass(np.array([10,30,50,70,90])),decimals=2)) 
    axin2.set_ylabel('Airmass [sec(z)]',rotation=-90,va='bottom')
    
    ### Local Time on top axis
    axin3=axin.twiny(); 
    #axin3.set_xlim(axin.get_xlim()) #--> This actually returns slightly different lims!!  Not [0,90]
    if timezone is not None and 'none' not in timezone.lower():
        #if True in [i in timezone.lower() for i in ['auto','calc']]:
        #    #Use tzwhere to compute the timezone based on the observer lat/lon (input in degrees)
        #    try: timezone=tzwhere.tzwhere().tzNameAt(observer.lat*180/np.pi, wrap_pm180(observer.lon*180/np.pi))
        #    except: timezone=tzwhere.tzwhere(forceTZ=True).tzNameAt(observer.lat*180/np.pi, wrap_pm180(observer.lon*180/np.pi), forceTZ=True)
        timezone = tz_from_observer(observer)
        utcoffset=calculate_dtnaive_utcoffset(ephem.Date(meanobstime).datetime(),local_timezone=timezone)
        times_local=np.zeros_like(times_utc)
        for i in range(len(times_utc)): times_local[i]=times_utc[i].astimezone(pytz.timezone(timezone))
        axin3.set_xlim(times_local[0],times_local[-1],auto=True) #Local time... still need to set tz set below...
        axin3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M',tz=times_local[0].tzinfo)) #NEED to set tz for axis update!
        axin3.set_xlabel('Local Time $\\rightarrow$', x=-.08, va='top') 
        axin3.text(-.04, 0.95, '(%s)\nUTC offset = %i'%(times_local[0].tzinfo.zone,utcoffset),va='bottom', ha='right', transform=axin3.transAxes, color='k', fontsize=8)

    ###LMST, Full axis calculation
    #axin3.set_xlim(times_lmst[0],times_lmst[-1])# L(M)ST
    #axin3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    ### LMST, calculate each major ticklabel
    ##axin3.set_xlim(times_abs[0],times_abs[-1]); axin3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M')); 
    #fig1.canvas.draw()#For some reason, forcing it to draw is necessary to make the ticklabels grabbable.
    #utctime0=compute_sidereal_time(wht,times_abs[0],as_type='datetime')
    #utc_ticks=axin3.get_xticklabels(); xtl=[]
    #for t in utc_ticks:
    #    t_utc=datetime.datetime.strptime(obsday+' '+t.get_text(),'%Y/%m/%d %H:%M').replace(tzinfo=pytz.utc) #need the Y/M/D!
    #    if t_utc.time()<utctime0.time(): t_utc=t_utc.replace(day=t_utc.day+1)
    #    #t_lst=LMST_from_UTC(t_utc,wht.lon*180./np.pi)
    #    t_lst=LST_from_local(t_utc,wht.lon*180./np.pi)
    #    #t_lst=compute_sidereal_time(wht,t_utc,as_type='dms')
    #    xtl.append('%i:%i'%(t_lst[0],t_lst[1]))
    #axin3.xaxis.set_ticklabels(xtl)
    #axin3.set_xlabel('LMST $\\rightarrow$',x=-.05,va='top') 
    
def plot_visibility_tracks(target_list,observer,obsstart,obsend, weights=None, mode='nearest', duration_hours=0, plotmeantransit=False, timezone='auto', xaxisformatter=mdates.DateFormatter('%H:%M'), figsize=(14,8), dpi=200,savepath='',showplot=False):
    """
    Plot visibility tracks for targets observed from a single observer station, for a given start and end time.  Creates the figure and axis, calls plot_visibility_tracks_toaxis() to do the plotting, then saves/shows.
    
    Parameters
    ----------
    target_list : ephem.FixedBody(), or array-like [list, tuple, np.array...] 
        List of pyephem FixedBody sky sources.  Can also input a single target source.
    observer : ephem.Observer() 
        Observatory/antenna/location object to use for the calculations
    obsstart : emphem.Date() or str [formatted as 'YYYY/MM/DD HH:MM:SS.s']
        Observation start time
    obsend : emphem.Date() or str [formatted as 'YYYY/MM/DD HH:MM:SS.s']
        Observation end time
    weights : np.array, list, or other array-like
        Optional (float) weights to give to certain stations.  (e.g. to prioritize certain baselines)   Larger values give higher importance. 
    mode : str
        Specifies which transit time to return. Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    duration_hours : float
        The desired duration of observations, for calculating the optimal start time.  If greater than zero, calculate start time and plot a vertical line denoting it. 
    plotmeantransit : bool 
        Set to True to plot a vertical line denoting the mean of the target transit times
    timezone : str or None
        Timezone (standard Olson names) for local time on upper x-axis.  Options: \n
            'none' or None : Do not add local time info to the upper axis
            string of a standard Olson database name [e.g., 'US/Mountain' or 'America/Chicago'] :  use this directly for dt calculations \n
            'auto' or 'calculate' : compute the timezone from the observer lat/lon using module tzwhere.  Takes a while (but is convenient!)
    xaxisformatter: mdates.DateFormatter
        Formatter object for the xaxis dates.  Default mdates.DateFormatter('%b') will print abbreviated month names, '%B' will print full month names.    
    figsize : tuple
        matplotlib Figure size
    dpi : int
        matplotlib image resolution (dpi = dots per inch)
    savepath : str
        The path for saving the figure, including the filename and extension. 
    showplot : bool
        Set to True to invoke plt.show() and display the interactive plot
    """
    fig1=plt.figure(1,figsize=figsize)
    ax1=fig1.add_subplot(111)#,adjustable='box-forced',aspect=5e-3)
    plot_visibility_tracks_toaxis(target_list,observer,obsstart,obsend, ax1, weights=weights, mode=mode, duration_hours=duration_hours, plotmeantransit=plotmeantransit, timezone=timezone, xaxisformatter=xaxisformatter)
    if savepath is not '': plt.savefig(savepath,bbox_inches='tight',dpi=dpi)
    if showplot==True: plt.show()
    plt.clf(); plt.close('all')


def plot_VLBA_visibility_tracks_toaxis(target,obsstart,obsend, axin, weights=None, mode='nearest', duration_hours=0, xaxisformatter=hourformat):
    """
    Plot visibility tracks for the VLBA stations, for a given target and obsdate. 
    Includes the mean of the transit times from each station in the full VLBA array. 
    Plots onto a user-supplied matplotlib/pyplot axis (useful for creating the plot & allowing further manual editing). 
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The target observation source
    obsstart : emphem.Date() or str [formatted as 'YYYY/MM/DD HH:MM:SS.s']
        Observation start time
    obsend : emphem.Date() or str [formatted as 'YYYY/MM/DD HH:MM:SS.s']
        Observation end time
    axin : matplotlib/pyplot axis object
        The axis on which to plot
    weights : np.array, list, or other array-like
        Optional (float) weights to give to certain stations.  (e.g. to prioritize certain baselines)   Larger values give higher importance. Passed to calculate_VLBA_mean_transit_time().  The order of stations is alphabetical: [BR,FD,HN,KP,LA,MK,NL,OV,PT,SC].
    mode : str
        Specifies which transit time to return. (passed to calculate_VLBA_mean_transit()) Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    duration_hours : float
        The desired duration of observations, for calculating the optimal start time.  If greater than zero, calculate start time and plot a vertical line denoting it. 
    xaxisformatter: mdates.DateFormatter
        Formatter object for the xaxis times.  Default obs.hourformat [mdates.DateFormatter('%H:%M')] prints the hour and minute.
    
    Example
    -------
    fig1=plt.figure(1,figsize=(14,8)) \n
    ax1=fig1.add_subplot(111) \n
    ngc1052=obs.create_ephem_target('NGC1052','02:41:04.7985','-08:15:20.751') \n
    obs.plot_VLBA_visibility_tracks_toaxis(ngc1052,ephem.Date('2021/04/15 00:00:00'),ephem.Date('2021/04/15 23:59:59'),ax1) \n
    plt.savefig('ngc1052_fullVLBA_april15.jpg'); plt.clf(); plt.close('all')
    """
    if type(obsstart) is str: obsstart=ephem.Date(obsstart)
    if type(obsend) is str: obsend=ephem.Date(obsend)
    times_utc=create_obstime_array(obsstart.datetime().strftime('%Y/%m/%d %H:%M:%S'),obsend.datetime().strftime('%Y/%m/%d %H:%M:%S'),timezone_string='UTC',n_steps=100)
    
    target_meanVLBAtransit=calculate_VLBA_mean_transit_time(target, obsstart, weights=weights, mode=mode)
    
    sunsepminmax,sunsepminmaxnames=sunsep_VLBAminmax(target,obsstart,return_names=True)
    
    #fig1=plt.figure(1,figsize=figsize)
    #ax1=fig1.add_subplot(111)#,adjustable='box-forced',aspect=5e-3)
    
    ### Plot the visibility tracks for each station position (calculate on the fly)
    for station,lab in zip([vlbaBR,vlbaFD,vlbaHN,vlbaKP,vlbaLA,vlbaMK,vlbaNL,vlbaOV,vlbaPT,vlbaSC],['BR','FD','HN','KP','LA','MK','NL','OV','PT','SC']):
        #target_altitudes_BR,target_azimuths_BR=compute_target_altaz(target,vlbaBR,obsstart,obsend,nsteps=100)
        axin.plot(times_utc,compute_target_altaz(target,station,obsstart,obsend,nsteps=100)[0],label=lab)
    
    ### Vertical line denoting the mean transit time among all the antennae
    axin.axvline(dt.datetime.strptime(target_meanVLBAtransit,'%Y/%m/%d %H:%M:%S'), ls='--',color='0.5',zorder=-1)
    axin.annotate('Mean transit time: %s (UTC)'%(target_meanVLBAtransit[-8:]),xy=[dt.datetime.strptime(target_meanVLBAtransit,'%Y/%m/%d %H:%M:%S'),85], xytext=[dt.datetime.strptime(target_meanVLBAtransit,'%Y/%m/%d %H:%M:%S')+dt.timedelta(minutes=-25),85], va='top', rotation=90,color='0.5')
    
    ### If duration_hours > 0, calculate the optimal start time and draw a vertical line denoting it
    if duration_hours>0: 
        starttime=calculate_optimal_VLBAstarttime(target,obsstart,duration_hours,return_fmt='dt', mode=mode)
        axin.axvline(starttime, ls=':',color='0.5',zorder=-1)
        axin.annotate('Optimal start time for %s hours duration on %s: %s (UTC)'%(duration_hours,dt.datetime.strftime(obsstart.datetime(),'%Y/%m/%d'),dt.datetime.strftime(starttime,'%H:%M:%S')), xy=[starttime,85], xytext=[starttime+dt.timedelta(minutes=-25),85], va='top', rotation=90, color='0.5')
    
    ### Fill a hatched area denoting the slew limit for the antennas.  According to the sched user manual Sec.3.3,
    #   This info can be found in $SCHED/catalogs/stations.dat with each station's elevation limit given by AXLIM2
    #   VLBA antennae have elevation limits [2.5,90].  VLA antennae are [8,90] ...
    #   
    axin.fill_between([times_utc[0],times_utc[-1]],[2.5,2.5],facecolor='none',hatch='XXXX',edgecolor='0.9',lw=0.0,zorder=-2)
    plt.grid(True,color='0.92')#,'both')
    
    axin.legend(loc='best')
    axin.set_title('%s elevations from VLBA stations on %s (UTC). Closest separation from Sun = %i deg [%s]'%(target.name,dt.datetime.strftime(obsstart.datetime(),'%Y/%m/%d'),sunsepminmax[0],sunsepminmaxnames[0]))
    axin.set_xlim(times_utc[0],times_utc[-1])
    axin.set_ylim(0,90)
    #axin.xaxis.set_major_formatter(mdates.DateFormatter('%H'))#,mdates.DateFormatter('%H:%M'))
    axin.xaxis.set_major_formatter(xaxisformatter)
    axin.set_xlabel('UTC Time'); plt.ylabel('Altitude [deg]')
    plt.tight_layout()


def plot_VLBA_visibility_tracks(target, obsstart, obsend, weights=None, mode='nearest', duration_hours=0, figsize=(14,8), dpi=200, savepath='', showplot=False):
    """
    Plot visibility tracks for the VLBA stations, for a given target and obsdate.
    Creates the figure and axis, calls plot_VLBA_visibility_tracks_toaxis() to do the plotting, then saves/shows.
    
    Parameters
    ----------
    target : ephem.FixedBody()
        The target observation source
    obsstart : emphem.Date() or str [formatted as 'YYYY/MM/DD HH:MM:SS.s']
        Observation start time
    obsend : emphem.Date() or str [formatted as 'YYYY/MM/DD HH:MM:SS.s']
        Observation end time
    weights : np.array, list, or other array-like
        Optional (float) weights to give to certain stations.  (e.g. to prioritize certain baselines)   Larger values give higher importance. Passed to calculate_VLBA_mean_transit_time().  The order of stations is alphabetical: [BR,FD,HN,KP,LA,MK,NL,OV,PT,SC].
    mode : str
        Specifies which transit time to return. Options are \n
        'previous'/'before' = calculate the last transit before the input time \n
        'next'/'after' = calculate the next transit after the input time \n
        'nearest' = calculate the nearest transit to the input time
    duration_hours : float
        The desired duration of observations, for calculating the optimal start time.  If greater than zero, calculate start time and plot a vertical line denoting it. 
    figsize : tuple
        matplotlib Figure size
    dpi : int
        matplotlib image resolution (dpi = dots per inch)
    savepath : str
        The path for saving the figure, including the filename and extension. 
    showplot : bool
        Set to True to invoke plt.show() and display the interactive plot
    """
    fig1=plt.figure(1,figsize=figsize)
    ax1=fig1.add_subplot(111)#,adjustable='box-forced',aspect=5e-3)
    plot_VLBA_visibility_tracks_toaxis(target,obsstart,obsend, ax1, weights=weights, mode=mode, duration_hours=duration_hours)
    if savepath is not '': plt.savefig(savepath,bbox_inches='tight',dpi=dpi)
    if showplot==True: plt.show()
    plt.clf(); plt.close('all')


##### ------------ Finder plots and survey cutouts ------------ #####

### General ###

def download_cutout(coords, boxwidth, savepath, survey='DSS2 Red', boxwidth_units='pix', frame='icrs', overwrite=False, search_name=False):
    """
    Downloads and saves a single fits image. If multiple surveys given, only first will be saved.
    SkyView.get_images() requires size in pixels.
    
    Parameters
    ----------
    coords : list or str
        The list of target [RA,DEC] coordinates to use in astroquery.skyview.SkyView.get_images, in decimal (not sexagesimal).  \n
        If search_name is True, this is a string name for the target source, to use with astroquery.
    boxwidth : int or float
        The width of the box sides of the desired image.  Units (e.g., pixels or arcminutes) specified with boxwidth_units
    savepath : str
        The path for saving the downloaded file, including the filename and desired filetype (.fits, .FIT, etc.)
    survey : str
        The survey to download from.  E.g., 'SDSSr','2MASS-J','WISE 3.4' ...  \n
        List all available surveys with: \n
        from astroquery.skyview import SkyView; SkyView.list_surveys()
    boxwidth_units: str
        The units of the specified box width.  Options: 'pixels'/'pix', 'degree','deg', 'arcmin'/'amin', 'arcsec'/'asec'
    frame : str
        The coordinate reference frame.  Defaults to 'icrs'. Other options include 'galactic','fk5', ...
    overwrite : bool
        Set to True to overwrite the existing file at the specified savepath \n
        If you have deleted a fits file and the old one keeps reappearing, try deleting the contents of your astroquery download cache, ~/.astopy/cache/download/url/ 
    search_name : bool
        Set to True to query SkyView by a common name, rather than by coordinates.  In this case, coords should be the name (string)
    
    Examples
    --------
    obs.download_cutout([83.63725, 22.0145], 10., 'Crab_2massJ_10arcmin.fits', survey='2MASS-J', boxwidth_units='arcmin') #Crab Nebula, M1 \n
    obs.download_cutout( obs.sex2dec('05:34:32.94','22:00:52.2'), 10., 'Crab_2massJ_10arcmin.fits', survey='2MASS-J', boxwidth_units='arcmin') \n
    obs.download_cutout('NGC1275', 90., './NGC1275_SDSSr_90asec.fits', survey='SDSSr', search_name=True, boxwidth_units='asec')
    """
    if overwrite == False and os.path.exists(savepath)==True: return
    
    if 'pix' in boxwidth_units.lower(): boxwidth_pix=int(np.ceil(boxwidth))
    else: 
        asec_per_pix={'SDSSg':.396, 'SDSSi':.396, 'SDSSr':.396, 'SDSSu':.396, 'SDSSz':.396, '2MASS-J':1., '2MASS-H':1., '2MASS-K':1., 'UKIDSS-J':0.4, 'UKIDSS-H':0.4, 'UKIDSS-K':0.4, 'DSS2 Blue':1., 'DSS2 Red':1., 'DSS2 IR':1., 'GALEX Near UV':1.5, 'GALEX Far UV':1.5, 'WISE 3.4':1.37484, 'WISE 4.6':1.37484, 'WISE 12':1.37484, 'WISE 22':1.37484}
        asec_per_unit={'arcsec':1., 'arcmin':60., 'degree':3600., 'asec':1., 'amin':60., 'deg':3600.}
        boxwidth_pix=int(np.ceil(boxwidth*asec_per_unit[boxwidth_units]/asec_per_pix[survey]))
    if search_name==True: 
        cutout_im=SkyView.get_images(position=coords,survey=survey,pixels=boxwidth_pix)[0][0]
    else: 
        coords_formatted = coordinates.SkyCoord(coords[0], coords[1], unit=('deg','deg'), frame=frame)
        cutout_im = SkyView.get_images(position=coords_formatted, survey=survey, pixels=boxwidth_pix)[0][0] 
    cutout_im.writeto(savepath, overwrite=overwrite)

def download_cutouts(coords, boxwidth, savepathstem, surveybands=['DSS2 Blue','DSS2 Red','DSS2 IR'], boxwidth_units='pix', frame='icrs', overwrite=False, search_name=False):
    """
    Downloads and saves multiple cutout fits images using astroquery.skyview.  Useful for, e.g., collecting images to make 3-color plots, or to quickly assemble plots in different bands.
    
    Parameters
    ----------
    coords : list or str
        The list of target [RA,DEC] coordinates to use in astroquery.skyview.SkyView.get_images, in decimal (not sexagesimal).  \n
        If search_name is True, this is a string name for the target source, to use with astroquery.
    boxwidth : int or float
        The width of the box sides of the desired image.  Units (e.g., pixels or arcminutes) specified with boxwidth_units
    savepathstem : str
        The path string stem to use for saving the downloaded files.  Filenames will be appended with '_<surveyband>.fits'
    boxwidth_units: str
        The units of the specified box width.  Options: 'pixels'/'pix', 'degree','deg', 'arcmin'/'amin', 'arcsec'/'asec'
    frame : str
        The coordinate reference frame.  Defaults to 'icrs'. Other options include 'galactic','fk5', ...
    overwrite : bool
        Set to True to overwrite the existing file at the specified savepath
    search_name : bool
        Set to True to query SkyView by a common name (e.g., 'M33').  Specify this name with the coords parameter instead of a list of [RA,DEC] coordinates.
    
    Example
    -------
    obs.download_cutouts('M104', 90., './M104_90asec', search_name=True, boxwidth_units='arcsec')
    """
    for surveyband in surveybands:
        savestring=savepathstem+'_%s.fits'%(surveyband).replace(' ','_')
        download_cutout(coords, boxwidth, savestring, survey=surveyband,  boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)




### SDSS ###

def download_sdss_cutouts(coords, boxwidth, savepathstem, boxwidth_units='pix', frame='icrs', overwrite=False, search_name=False, SDSSbands=['SDSSg','SDSSi','SDSSr']):
    """
    Convenience function to download multiple SDSS cutout images.  \n
    Calls download_cutouts() with SDSS bands being passed to surveybands.
    
    Parameters
    ----------
    coords : list or str
        The list of target [RA,DEC] coordinates to use in astroquery.skyview.SkyView.get_images, in decimal (not sexagesimal).  \n
        If search_name is True, this is a string name for the target source, to use with astroquery.
    boxwidth : int or float
        The width of the box sides of the desired image.  Units (e.g., pixels or arcminutes) specified with boxwidth_units
    savepathstem : str
        The path string stem to use for saving the downloaded files.  Filenames will be appended with '_SDSS<band>.fits'
    boxwidth_units: str
        The units of the specified box width.  Options: 'pixels'/'pix', 'degree','deg', 'arcmin'/'amin', 'arcsec'/'asec'
    frame : str
        The coordinate reference frame.  Defaults to 'icrs'. Other options include 'galactic','fk5', ...
    overwrite : bool
        Set to True to overwrite the existing file at the specified savepath
    search_name : bool
        Set to True to query SkyView by a common name (e.g., 'M33').  Specify this name with the coords parameter instead of a list of [RA,DEC] coordinates.
    SDSSbands : list
        The list of SDSS bands to download.  'SDSSg','SDSSi','SDSSr','SDSSu','SDSSz', 'SDSSdr7g','SDSSdr7i',... 
    
    Example
    -------
    obs.download_sdss_cutouts('M104', 90., './M104_90asec', search_name=True, boxwidth_units='arcsec')
    """
    download_cutouts(coords, boxwidth, savepathstem, surveybands=SDSSbands, boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)


### DSS ###

def download_dss_cutouts(coords, boxwidth, savepathstem, boxwidth_units='pix', frame='icrs', overwrite=False, search_name=False):
    """
    Convenience function to download cutouts of the DSS2 images: B, R, IR.  \n
    Calls download_cutouts() with DSS2 bands being passed to surveybands.
    
    Parameters
    ----------
    coords : list or str
        The list of target [RA,DEC] coordinates to use in astroquery.skyview.SkyView.get_images, in decimal (not sexagesimal).  \n
        If search_name is True, this is a string name for the target source, to use with astroquery.
    boxwidth : int or float
        The width of the box sides of the desired image.  Units (e.g., pixels or arcminutes) specified with boxwidth_units
    savepathstem : str
        The path string stem to use for saving the downloaded files.  Filenames will be appended with '_DSS2<band>.fits'
    boxwidth_units: str
        The units of the specified box width.  Options: 'pixels'/'pix', 'degree','deg', 'arcmin'/'amin', 'arcsec'/'asec'
    frame : str
        The coordinate reference frame.  Defaults to 'icrs'. Other options include 'galactic','fk5', ...
    overwrite : bool
        Set to True to overwrite the existing file at the specified savepath
    search_name : bool
        Set to True to query SkyView by a common name (e.g., 'M33').  Specify this name with the coords parameter instead of a list of [RA,DEC] coordinates.
    
    Example
    -------
    obs.download_dss_cutouts('M104', 90., './M104_90asec', search_name=True, boxwidth_units='arcsec')
    """
    download_cutouts(coords, boxwidth, savepathstem, surveybands=['DSS2 Blue','DSS2 Red','DSS2 IR'], boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)


### GALEX ###

def download_galex_cutouts(coords, boxwidth, savepathstem, boxwidth_units='pix', overwrite=False, search_name=False):
    """
    Convenience function to download cutouts of the GALEX NUV and FUV files.  \n
    Calls download_cutouts() with GALEX bands being passed to surveybands.
    
    Parameters
    ----------
    coords : list or str
        The list of target [RA,DEC] coordinates to use in astroquery.skyview.SkyView.get_images, in decimal (not sexagesimal).  \n
        If search_name is True, this is a string name for the target source, to use with astroquery.
    boxwidth : int or float
        The width of the box sides of the desired image.  Units (e.g., pixels or arcminutes) specified with boxwidth_units
    savepathstem : str
        The path string stem to use for saving the downloaded files.  Filenames will be appended with '_GALEX<band>.fits'
    boxwidth_units: str
        The units of the specified box width.  Options: 'pixels'/'pix', 'degree','deg', 'arcmin'/'amin', 'arcsec'/'asec'
    frame : str
        The coordinate reference frame.  Defaults to 'icrs'. Other options include 'galactic','fk5', ...
    overwrite : bool
        Set to True to overwrite the existing file at the specified savepath
    search_name : bool
        Set to True to query SkyView by a common name (e.g., 'M33').  Specify this name with the coords parameter instead of a list of [RA,DEC] coordinates.
    
    Example
    -------
    obs.download_galex_cutouts('M104', 90., './M104_90asec', search_name=True, boxwidth_units='arcsec')
    """
    download_cutout(coords, boxwidth, savepathstem+'_galexNUV.fits', survey='GALEX Near UV', boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)
    download_cutout(coords, boxwidth, savepathstem+'_galexFUV.fits', survey='GALEX Far UV', boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)


### WISE ###

def download_wise_cutouts(coords, boxwidth, savepathstem, boxwidth_units='pix', frame='icrs', overwrite=False, search_name=False):
    """
    Convenience function to download cutouts of the four WISE band images.  \n
    Calls download_cutouts() with GALEX bands being passed to surveybands.
    
    Parameters
    ----------
    coords : list or str
        The list of target [RA,DEC] coordinates to use in astroquery.skyview.SkyView.get_images, in decimal (not sexagesimal).  \n
        If search_name is True, this is a string name for the target source, to use with astroquery.
    boxwidth : int or float
        The width of the box sides of the desired image.  Units (e.g., pixels or arcminutes) specified with boxwidth_units
    savepathstem : str
        The path string stem to use for saving the downloaded files.  Filenames will be appended with '_WISE<band>.fits'
    boxwidth_units: str
        The units of the specified box width.  Options: 'pixels'/'pix', 'degree','deg', 'arcmin'/'amin', 'arcsec'/'asec'
    frame : str
        The coordinate reference frame.  Defaults to 'icrs'. Other options include 'galactic','fk5', ...
    overwrite : bool
        Set to True to overwrite the existing file at the specified savepath
    search_name : bool
        Set to True to query SkyView by a common name (e.g., 'M33').  Specify this name with the coords parameter instead of a list of [RA,DEC] coordinates.
    
    Example
    -------
    obs.download_wise_cutouts('M104', 90., './M104_90asec', search_name=True, boxwidth_units='arcsec')
    """
    download_cutouts(coords, boxwidth, savepathstem, surveybands=['WISE 3.4','WISE 4.6','WISE 12','WISE 22'], boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)


### J,H,K ###

def download_2mass_cutouts(coords, boxwidth, savepathstem, boxwidth_units='pix', frame='icrs', overwrite=False, search_name=False):
    """
    Convenience function to download cutouts of the three 2MASS images, in J,H,K bands.  \n
    Calls download_cutouts() with GALEX bands being passed to surveybands.
    
    Parameters
    ----------
    coords : list or str
        The list of target [RA,DEC] coordinates to use in astroquery.skyview.SkyView.get_images, in decimal (not sexagesimal).  \n
        If search_name is True, this is a string name for the target source, to use with astroquery.
    boxwidth : int or float
        The width of the box sides of the desired image.  Units (e.g., pixels or arcminutes) specified with boxwidth_units
    savepathstem : str
        The path string stem to use for saving the downloaded files.  Filenames will be appended with '_2MASS<band>.fits'
    boxwidth_units: str
        The units of the specified box width.  Options: 'pixels'/'pix', 'degree','deg', 'arcmin'/'amin', 'arcsec'/'asec'
    frame : str
        The coordinate reference frame.  Defaults to 'icrs'. Other options include 'galactic','fk5', ...
    overwrite : bool
        Set to True to overwrite the existing file at the specified savepath
    search_name : bool
        Set to True to query SkyView by a common name (e.g., 'M33').  Specify this name with the coords parameter instead of a list of [RA,DEC] coordinates.
    
    Example
    -------
    obs.download_2mass_cutouts('M104', 90., './M104_90asec', search_name=True, boxwidth_units='arcsec')
    """
    download_cutouts(coords, boxwidth, savepathstem, surveybands=['2MASS-J','2MASS-H','2MASS-K'], boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)

def download_ukidss_cutouts(coords, boxwidth, savepathstem, boxwidth_units='pix', frame='icrs', overwrite=False, search_name=False):
    """
    Downloads cutouts of the three UKIDSS bands: J,H,K
    
    Parameters
    ----------
    coords : list or str
        The list of target [RA,DEC] coordinates to use in astroquery.skyview.SkyView.get_images, in decimal (not sexagesimal).  \n
        If search_name is True, this is a string name for the target source, to use with astroquery.
    boxwidth : int or float
        The width of the box sides of the desired image.  Units (e.g., pixels or arcminutes) specified with boxwidth_units
    savepathstem : str
        The path string stem to use for saving the downloaded files.  Filenames will be appended with '_UKIDSS<band>.fits'
    boxwidth_units: str
        The units of the specified box width.  Options: 'pixels'/'pix', 'degree','deg', 'arcmin'/'amin', 'arcsec'/'asec'
    frame : str
        The coordinate reference frame.  Defaults to 'icrs'. Other options include 'galactic','fk5', ...
    overwrite : bool
        Set to True to overwrite the existing file at the specified savepath
    search_name : bool
        Set to True to query SkyView by a common name (e.g., 'M33').  Specify this name with the coords parameter instead of a list of [RA,DEC] coordinates.
    
    Example
    -------
    obs.download_ukidss_cutouts('M104', 90., './M104_90asec', search_name=True, boxwidth_units='arcsec')
    """
    download_cutouts(coords, boxwidth, savepathstem, surveybands=['UKIDSS-J','UKIDSS-H','UKIDSS-K'], boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)


def download_ukidss_cutouts2(coords, boxwidth_asec, savepathstem, overwrite=False):
    """
    Alternative method to download UKIDSS images using astroquery.ukidss instead of astroquery.skyview
    """
    from astroquery.ukidss import Ukidss
    import astropy.units as u
    coords_formatted = coordinates.SkyCoord(coords[0], coords[1], unit=('deg','deg'), frame='icrs')
    cutout_imJ = Ukidss.get_images(coords_formatted, waveband='J', radius=boxwidth_asec*u.arcmin)[0][0] 
    cutout_imJ.writeto(savepathstem+'_ukidssJ.fits', overwrite=overwrite)
    cutout_imH = Ukidss.get_images(coords_formatted, waveband='H', radius=boxwidth_asec*u.arcmin)[0][0] 
    cutout_imH.writeto(savepathstem+'_ukidssH.fits', overwrite=overwrite)
    cutout_imK = Ukidss.get_images(coords_formatted, waveband='K', radius=boxwidth_asec*u.arcmin)[0][0] 
    cutout_imK.writeto(savepathstem+'_ukidssK.fits', overwrite=overwrite)


### Producing the finder plots

def add_sizebar(axin, length_pixels, label, loc=4, sep=5, borderpad=0.8, frameon=False, path_effects=[PathEffects.withStroke(foreground='k', linewidth=1.75)], color=None, **kwargs):
    """
    Add a sizebar to the specified matplotlib axis
    """
    ### See matplotlib anchoredartist tutorial:  http://matplotlib.org/examples/pylab_examples/anchored_artists.html
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse, AnchoredSizeBar
    try: from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText  #Matplotlib <2.1
    except: from matplotlib.offsetbox import AnchoredText                   #Matplotlib >=2.1
    #from mpl_toolkits.axes_grid.anchored_artists import AnchoredDrawingArea
    #asb =  AnchoredSizeBar(axin.transData,length_pixel,label,loc=loc,borderpad=borderpad, sep=sep,frameon=frameon,**kwargs)
    asb =  AnchoredSizeBar(axin.transData, length_pixels, label, loc=loc, borderpad=borderpad, sep=sep, frameon=frameon, **kwargs)
    if color is not None:
        scalebar = asb.size_bar.get_children()[0]
        scalebar.set_ec(color)
        asb.txt_label._text.set_color(color)
    if path_effects!=[PathEffects.withStroke(foreground='k', linewidth=1.75)]:
        for a in [asb.size_bar._children[0],asb.txt_label._text]: a.set_path_effects(path_effects)
    axin.add_artist(asb)
    return asb

def bgpars2pix(arrin,headerin,precise=True):
    """
    ### Gets the regions in terms of image coordinates: 
    takes array [RA (dec), DEC (dec), Rmaj (deg), Rmin (deg), PA (deg., N from W)] \n
    Converts these angular (degree) RA,DEC,Rmaj,Rmin values to pixel values using the supplied fits header.
    """
    cpix=mcf.convsky2pix(headerin,arrin[0],arrin[1],precise=precise); 
    try: cd2=mcf.getcdelts(headerin)[1]
    except: raise(Exception('Invalid CDELTs in WCS'))
    return [cpix[0],cpix[1],arrin[2]/cd2,arrin[3]/cd2,arrin[4]]

def make_finder_plot_singleband(targetname, coords, boxwidth, survey='DSS2 Red', boxwidth_units='pix', frame='icrs', overwrite=False, search_name=False, refregs=None, cmap='gist_yarg', dpi=200, tickcolor='0.2', mfc='r', mec='w', bs_amin=0., surveytitle=False, savepathstem='./', filetype='jpg'):
    """
    Downloads a cutout image (with astroquery.skyview, if not already on disk) of the specified survey band, and creates a finder plot.  Optional refregs can be specified to plot reference positions on the plot (notable guide stars, etc).
    
    Parameters
    ----------
    targetname : str
        The name of the target, to use in filenames and the plot title
    coords : list or str
        The list of target [RA,DEC] coordinates to use in astroquery.skyview.SkyView.get_images, in decimal (not sexagesimal).  \n
        If search_name is True, this is a string name for the target source, to use with astroquery.
    boxwidth : int or float
        The width of the box sides of the desired image.  Units (e.g., pixels or arcminutes) specified with boxwidth_units
    boxwidth_units: str
        The units of the specified box width.  Options: 'pixels'/'pix', 'degree','deg', 'arcmin'/'amin', 'arcsec'/'asec'
    survey : str
        The survey to download from.  E.g., 'SDSSr', '2MASS-J', 'WISE 3.4' ...  \n
        List all available surveys with: \n
        from astroquery.skyview import SkyView; SkyView.list_surveys()
    frame : str
        The coordinate reference frame.  Defaults to 'icrs'. Other options include 'galactic','fk5', ...
    overwrite : bool
        Set to True to overwrite the existing file at the specified savepath \n
        If you have deleted a fits file and the old one keeps reappearing, try deleting the contents of your astroquery download cache, ~/.astopy/cache/download/url/ \n
        Or, astropy.utils.data.clear_download_cache()
    search_name : bool
        Set to True to query SkyView by a common name (e.g., 'M33').  Specify this name with the coords parameter instead of a list of [RA,DEC] coordinates.
    refregs : list
        refregs = list of reference position points in format [RA,DEC,'Name']
    cmap : str or matplotlib LinearSegmentedColormap
        The matplotlib colormap (or cmap specifier string) to apply to the image.  
    dpi : int
        The raster resolution (dots per inch) to apply to the saved image
    tickcolor : str or matplotlib color specifier 
        The color to use for the tick marks
    mfc : str or matplotlib color specifier
        Marker face color
    mec : str or matplotlib color specifier
        Marker edge color
    bs_amin : float
        The (angular) bar size to plot, in arcmin 
    surveytitle : bool
        Set to True to include the survey used for the background image
    savepathstem : str
        The path string stem to use for saving the downloaded file.  
    filetype : str
        The filetype to use for the saved image.  'jpg','png','pdf','svg','eps'...
    
    Example
    -------
    refstars=[['5:34:42.3','22:10:34.5','HD 244988'], ['5:34:36.6','21:37:19.9','HD 36707'],]
    obs.make_finder_plot_singleband('Crab Nebula', 'M1', 50., boxwidth_units='arcmin', survey='DSS2 Red', search_name=True, refregs=refstars, cmap='gist_yarg', dpi=200, tickcolor='0.2', mfc='r', mec='w', bs_amin=10., )
    """
    
    rcparams_initial = {'xtick.major.size':rcParams['xtick.major.size'], 'ytick.major.size':rcParams['ytick.major.size'], 'xtick.minor.size':rcParams['xtick.minor.size'], 'ytick.minor.size':rcParams['ytick.minor.size'], 'grid.linestyle':rcParams['grid.linestyle'], 'grid.linewidth':rcParams['grid.linewidth'], 'xtick.direction':rcParams['xtick.direction'], 'ytick.direction':rcParams['ytick.direction'], 'xtick.minor.visible':rcParams['xtick.minor.visible'], 'ytick.minor.visible':rcParams['ytick.minor.visible']}
    #{'grid.alpha':rcParams['grid.alpha'], 'grid.color':rcParams['grid.color'], 'grid.linestyle':rcParams['grid.linestyle'], 'grid.linewidth':rcParams['grid.linewidth']}
    rcparams_update={'xtick.major.size':7, 'ytick.major.size':7, 'xtick.minor.size':3, 'ytick.minor.size':3, 'grid.linestyle':':', 'grid.linewidth':0.3, 'xtick.direction':'in', 'ytick.direction':'in', 'xtick.minor.visible':True, 'ytick.minor.visible':True}
    rcParams.update(rcparams_update)
    
    fitssavestring=savepathstem+targetname+'_'+survey+'.fits'; 
    fitssavestring=fitssavestring.replace(' ','_')
    download_cutout(coords, boxwidth, fitssavestring, survey=survey, boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)
    
    cutout_dat,cutout_hdr = pyfits.getdata(fitssavestring, header=True);
    
    fig1=plt.figure(1)#,figsize=(9,3.5))
    if surveytitle is True: fig1.suptitle(targetname+'  (%s)'%survey, fontsize=14, x=.5, y=.98) # single row
    else: fig1.suptitle(targetname, fontsize=14, x=.5, y=.98) 
    
    cutout_wcs = WCS(cutout_hdr)

    ax1=fig1.add_subplot(111, projection=cutout_wcs)
    ax1.imshow(cutout_dat,origin='lower',interpolation='nearest',cmap=cmap)
    rapars = ax1.coords[0]
    decpars = ax1.coords[1]
    rapars.set_ticks(color=tickcolor); decpars.set_ticks(color=tickcolor)
    rapars.set_ticks(number=6); decpars.set_ticks(number=6); 
    rapars.set_ticklabel(size=8); decpars.set_ticklabel(size=8); 
    #rapars.set_ticks(spacing=10*u.arcmin, color='white', exclude_overlapping=True)
    #decpars.set_ticks(spacing=5*u.arcmin, color='white', exclude_overlapping=True)
    rapars.display_minor_ticks(True); #rapars.set_minor_frequency(10)
    decpars.display_minor_ticks(True)
    rapars.set_major_formatter('hh:mm:ss.s')
    decpars.set_major_formatter('dd:mm')
    rapars.set_separator(('$^\mathrm{H}$', "'", '"'))
    #decpars.ticklabels.set_rotation(45) #Rotate ticklabels
    #decpars.ticklabels.set_color(xkcdrust) #Ticklabel Color
    
    ax1.set_xlabel('RA (J2000)'); ax1.set_ylabel('DEC (J2000)');
    
    if type(coords) is str: 
        tpx=[cutout_dat.shape[1]/2,cutout_dat.shape[0]/2]
    else:
        tpx=mcf.convsky2pix(cutout_hdr,coords[0],coords[1],precise=True)
    ax1.plot(tpx[0],tpx[1],marker='+',mec=mfc,ms=8,zorder=10)
    
    if refregs!=None: 
        for reg in refregs:
            if type(reg[0])==str: regcoords=sex2dec(reg[0],reg[1])
            else: regcoords=reg[:2]
            rpx=bgpars2pix([regcoords[0],regcoords[1],(boxwidth/80.)/60.,(boxwidth/80.)/60.,0.],cutout_hdr)
            c1=patches.Ellipse((rpx[0],rpx[1]),2*rpx[2],2*rpx[3],angle=rpx[4], linewidth=.8, ec=mfc, fill=False, zorder=5, path_effects=[PathEffects.withStroke(linewidth=1.2,foreground=mec)])
            #plt.setp(c1,path_effects=[PathEffects.withStroke(linewidth=1,foreground=mec)])
            ax1.add_patch(c1)
            ax1.text(rpx[0],rpx[1]+rpx[2]*1.1,reg[-1],va='bottom',ha='center', color=mfc, fontsize=8, path_effects=[PathEffects.withStroke(linewidth=.3,foreground=mec)],zorder=6)
    
    ax1.grid(b=True,which='major',color='0.2',linestyle='-',linewidth=0.4)#,alpha=0.8)
    
    ### Size Bar
    if bs_amin>0:
        #bs_amin = Length of inset size bar, in arcmin
        barsizepix=bs_amin/(60.*mcf.getdegperpix(cutout_hdr))
        dbar=add_sizebar(ax1, barsizepix, '%i $^{\prime}$'%(bs_amin), loc=3, sep=5, borderpad=0.3, frameon=False, path_effects=[PathEffects.withStroke(foreground=mec, linewidth=1.2)], color=mfc)
    
    figsavestring=savepathstem+targetname+'_'+survey.replace(' ','')+'.'+filetype
    plt.savefig(figsavestring.replace(' ','_'),dpi=dpi,bbox_inches='tight')
    plt.clf(); plt.close()
    
    rcParams.update(rcparams_initial) #Revert back so as not to disrupt previous settings


def make_finder_plot_simpleRGB(targetname, coords, boxwidth, surveyR, surveyG, surveyB, boxwidth_units='pix', frame='icrs', overwrite=False, search_name=False, refregs=None, Rscalefunc='sqrt', Gscalefunc='sqrt', Bscalefunc='sqrt', dpi=200, tickcolor='0.2', mfc='r', mec='w', bs_amin=0., surveytitle=False, savepathstem='./', filetype='jpg'):
    """
    Downloads a cutout image (with astroquery.skyview, if not already on disk) of the specified survey band, and creates a simple R,G,B stacked image finder plot.  Optional refregs can be specified to plot reference positions on the plot (notable guide stars, etc).
    
    Parameters
    ----------
    targetname : str
        The name of the target, to use in filenames and the plot title
    coords : list or str
        The list of target [RA,DEC] coordinates to use in astroquery.skyview.SkyView.get_images, in decimal (not sexagesimal).  \n
        If search_name is True, this is a string name for the target source, to use with astroquery.
    boxwidth : int or float
        The width of the box sides of the desired image.  Units (e.g., pixels or arcminutes) specified with boxwidth_units
    surveyR : str
        The Red-channel survey to download from.  E.g., 'SDSSr', '2MASS-J', 'WISE 3.4' ...  \n
        List all available surveys with: \n
        from astroquery.skyview import SkyView; SkyView.list_surveys()
    surveyG : str 
        The Green-channel survey to download from.  E.g., 'SDSSr', '2MASS-J',' WISE 3.4' ...
    surveyB : str
        The Blue-channel survey to download from.  E.g., 'SDSSr', '2MASS-J', 'WISE 3.4' ...
    boxwidth_units: str
        The units of the specified box width.  Options: 'pixels'/'pix', 'degree','deg', 'arcmin'/'amin', 'arcsec'/'asec'
    frame : str
        The coordinate reference frame.  Defaults to 'icrs'. Other options include 'galactic','fk5', ...
    overwrite : bool
        Set to True to overwrite the existing file at the specified savepath \n
        If you have deleted a fits file and the old one keeps reappearing, try deleting the contents of your astroquery download cache, ~/.astopy/cache/download/url/ 
    search_name : bool
        Set to True to query SkyView by a common name (e.g., 'M33').  Specify this name with the coords parameter instead of a list of [RA,DEC] coordinates.
    refregs : list
        refregs = list of reference position points in format [RA,DEC,'Name']
    Rscalefunc : str
        The scaling function from multicolorfits to use for scaling in the Red channel: 'linear', 'sqrt', 'log', 'asinh' 
    Gscalefunc : str
        The scaling function from multicolorfits to use for scaling in the Green channel: 'linear', 'sqrt', 'log', 'asinh'
    Bscalefunc : str 
        The scaling function from multicolorfits to use for scaling in the Blue channel: 'linear', 'sqrt', 'log', 'asinh'
    dpi : int
        The raster resolution (dots per inch) to apply to the saved image
    tickcolor : str or matplotlib color specifier 
        The color to use for the tick marks
    mfc : str or matplotlib color specifier
        Marker face color
    mec : str or matplotlib color specifier
        Marker edge color
    bs_amin : float
        The (angular) bar size to plot, in arcmin 
    surveytitle : bool
        Set to True to include the surveys used for the background image
    savepathstem : str
        The path string stem to use for saving the downloaded file.  
    filetype : str
        The filetype to use for the saved image.  'jpg','png','pdf','svg','eps'...
    Example
    -------
    refstars=[['5:34:42.3','22:10:34.5','HD 244988'], ['5:34:36.6','21:37:19.9','HD 36707'],]
    obs.make_finder_plot_simpleRGB('M1', obs.sex2dec('05:34:31.94','22:00:52.2'), 50., 'DSS2 IR','DSS2 Red','DSS2 Blue', boxwidth_units='arcmin', refregs=refstars, Rscalefunc='linear', Gscalefunc='linear', Bscalefunc='linear', dpi=200, tickcolor='0.8', mfc='#C11B17',mec='w',bs_amin=10., filetype='jpg')
    """
    rcparams_initial = {'xtick.major.size':rcParams['xtick.major.size'], 'ytick.major.size':rcParams['ytick.major.size'], 'xtick.minor.size':rcParams['xtick.minor.size'], 'ytick.minor.size':rcParams['ytick.minor.size'], 'grid.linestyle':rcParams['grid.linestyle'], 'grid.linewidth':rcParams['grid.linewidth'], 'xtick.direction':rcParams['xtick.direction'], 'ytick.direction':rcParams['ytick.direction'], 'xtick.minor.visible':rcParams['xtick.minor.visible'], 'ytick.minor.visible':rcParams['ytick.minor.visible']}
    rcparams_update={'xtick.major.size':7, 'ytick.major.size':7, 'xtick.minor.size':3, 'ytick.minor.size':3, 'grid.linestyle':':', 'grid.linewidth':0.3, 'xtick.direction':'in', 'ytick.direction':'in', 'xtick.minor.visible':True, 'ytick.minor.visible':True}
    rcParams.update(rcparams_update)
    
    fitssavestringR=savepathstem+targetname+'_'+surveyR+'.fits'
    fitssavestringR=fitssavestringR.replace(' ','_')
    fitssavestringG=savepathstem+targetname+'_'+surveyG+'.fits'
    fitssavestringG=fitssavestringG.replace(' ','_')
    fitssavestringB=savepathstem+targetname+'_'+surveyB+'.fits'
    fitssavestringB=fitssavestringB.replace(' ','_')
    
    download_cutout(coords, boxwidth, fitssavestringR, survey=surveyR, boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)
    download_cutout(coords, boxwidth, fitssavestringG, survey=surveyG, boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)
    download_cutout(coords, boxwidth, fitssavestringB, survey=surveyB, boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)
    
    cutoutR_dat,cutoutR_hdr=pyfits.getdata(fitssavestringR, header=True);
    cutoutG_dat,cutoutG_hdr=pyfits.getdata(fitssavestringG, header=True);
    cutoutB_dat,cutoutB_hdr=pyfits.getdata(fitssavestringB, header=True);
    
    #All dl files should have the same spatial grid.  Otherwise, mcf.reproject_to()...
    
    R_scale=(mcf.scaling_fns[Rscalefunc]() + mcf.ManualInterval(vmin=np.nanmin(cutoutR_dat), vmax=np.nanmax(cutoutR_dat)))(cutoutR_dat)
    G_scale=(mcf.scaling_fns[Gscalefunc]() + mcf.ManualInterval(vmin=np.nanmin(cutoutG_dat), vmax=np.nanmax(cutoutG_dat)))(cutoutG_dat)
    B_scale=(mcf.scaling_fns[Bscalefunc]() + mcf.ManualInterval(vmin=np.nanmin(cutoutB_dat), vmax=np.nanmax(cutoutB_dat)))(cutoutB_dat)
    
    RGBimage = np.dstack([ R_scale, G_scale, B_scale])
    plothdr=cutoutR_hdr.copy()
    
    fig1=plt.figure(1)#,figsize=(9,3.5))
    if surveytitle is True: 
        fig1.suptitle(targetname+'  (%s, %s, %s)'%(surveyR,surveyG,surveyB), fontsize=12, x=.5, y=.98) 
    else: fig1.suptitle(targetname, fontsize=12, x=.5, y=.98) 
    cutout_wcs = WCS(plothdr)

    ax1=fig1.add_subplot(111, projection=cutout_wcs)
    ax1.imshow(RGBimage,origin='lower',interpolation='nearest')
    rapars = ax1.coords[0]
    decpars = ax1.coords[1]
    rapars.set_ticks(color=tickcolor); decpars.set_ticks(color=tickcolor)
    rapars.set_ticks(number=6); decpars.set_ticks(number=6); 
    rapars.set_ticklabel(size=8); decpars.set_ticklabel(size=8); 
    #rapars.set_ticks(spacing=10*u.arcmin, color='white', exclude_overlapping=True)
    #decpars.set_ticks(spacing=5*u.arcmin, color='white', exclude_overlapping=True)
    rapars.display_minor_ticks(True); #rapars.set_minor_frequency(10)
    decpars.display_minor_ticks(True)
    rapars.set_major_formatter('hh:mm:ss.s')
    decpars.set_major_formatter('dd:mm')
    rapars.set_separator(('$^\mathrm{H}$', "'", '"'))
    #decpars.ticklabels.set_rotation(45) #Rotate ticklabels
    #decpars.ticklabels.set_color(xkcdrust) #Ticklabel Color
    
    #ax1.set_xlabel("$\Delta$ RA",size=8)
    #ax1.set_ylabel("$\Delta$ DEC",size=8)
    ax1.set_xlabel('RA (J2000)'); ax1.set_ylabel('DEC (J2000)');
    
    if type(coords) is str: 
        tpx=[cutoutR_dat.shape[1]/2,cutoutR_dat.shape[0]/2]
    else:
        tpx=mcf.convsky2pix(cutoutR_hdr,coords[0],coords[1],precise=True)
    ax1.plot(tpx[0],tpx[1],marker='+',mec=mfc,ms=8,zorder=10)
    
    if refregs!=None: 
        for reg in refregs:
            if type(reg[0])==str: regcoords=sex2dec(reg[0],reg[1])
            else: regcoords=reg[:2]
            rpx = bgpars2pix( [regcoords[0], regcoords[1], (boxwidth/80.)/60., (boxwidth/80.)/60., 0.], plothdr)
            c1 = patches.Ellipse( (rpx[0],rpx[1]), 2*rpx[2], 2*rpx[3], angle=rpx[4], linewidth=.8, ec=mfc, fill=False, zorder=5, path_effects=[PathEffects.withStroke(linewidth=1.2,foreground=mec)] )
            #plt.setp( c1, path_effects=[PathEffects.withStroke(linewidth=1,foreground=mec)] )
            ax1.add_patch(c1)
            ax1.text(rpx[0], rpx[1]+rpx[2]*1.1, reg[-1], va='bottom', ha='center', color=mfc, fontsize=8, path_effects=[PathEffects.withStroke(linewidth=.3,foreground=mec)], zorder=6 )
    
    ax1.grid(b=True, which='major', color='0.8', linestyle='-', linewidth=0.4)#,alpha=0.8)
    
    ### Size Bar
    if bs_amin>0:
        #bs_amin = Length of inset size bar, in arcmin
        barsizepix=bs_amin/(60.*mcf.getdegperpix(plothdr))
        dbar=add_sizebar(ax1, barsizepix, '%i $^{\prime}$'%(bs_amin), loc=3, sep=5, borderpad=0.3, frameon=False, path_effects=[PathEffects.withStroke(foreground=mec, linewidth=1.2)], color=mfc)
    
    figsavestring=savepathstem+targetname+'_RGB'.replace(' ','')+'.'+filetype
    plt.savefig(figsavestring.replace(' ','_'),dpi=dpi,bbox_inches='tight')
    plt.clf(); plt.close()
    
    rcParams.update(rcparams_initial) #Revert back so as not to disrupt previous settings
    

#SDSSrlim=[-.001,2]; SDSSglim=[-.001,2]; SDSSulim=[-.001,2] #Suggested min/max values for SDSS R,G,U bands
#JHKjlim=[None,None]; JHKhlim=[None,None]; JHKklim=[None,None]
'''
def plotsimpleRGB(Rdata,Rhdr, Gdata,Ghdr, Bdata,Bhdr, Rminmax=None, Gminmax=None, Bminmax=None, projecttoband='R', Rscalefunc='sqrt', Gscalefunc='sqrt', Bscalefunc='sqrt', alphachan=1, ax=None, overlayhdr=None):
    """
    Rdata, Gdata, Bdata = 2D fits data files loaded with astropy.io.fits  -- trim to 2D slices beforehand
    Rhdr, Ghdr, Bhdr = fits headers
    Rminmax,Gminmax,Bminmax = the [min,max] list of values for each band, for scaling.
    projecttoband = the band to reproject to: 'R', 'G', or 'B'.  Default is R
    Rscalefunc, Gscalefunc, Bscalefunc = the scaling function from multicolorfits to use for scaling luminance: 'linear', 'sqrt', 'log', 'asinh'
    ax = the axis object to plot to (if None, will plot in new 1,1,1 subplot) 
    overlayhdr = the header used in pywcsgrid2 for plotting to grid (should be Rhdr, Ghdr, or Bhdr -- not the grid_helper header)
    alphachan = fraction 0 to 1 of overall alpha value
    
    Example
    -------
    Rdata,Rhdr=pyfits.getdata('M101_90asec_SDSSr.fits',header=True)
    Gdata,Ghdr=pyfits.getdata('M101_90asec_SDSSi.fits',header=True)
    Bdata,Bhdr=pyfits.getdata('M101_90asec_SDSSg.fits',header=True)
    obs.plotsimpleRGB(Rdata,Rhdr, Gdata,Ghdr, Bdata,Bhdr, Rminmax=None, Gminmax=None, Bminmax=None, projecttoband='R', Rscalefunc='log', Gscalefunc='log', Bscalefunc='log', alphachan=1, ax=None, overlayhdr=None)
    plt.show(); plt.clf()
    """
    import multicolorfits as mcf
    if projecttoband=='R': 
        R_reproj=Rdata.copy()
        G_reproj=mcf.reproject2D(Gdata,Ghdr,Rhdr)
        B_reproj=mcf.reproject2D(Bdata,Bhdr,Rhdr)
        if overlayhdr==None: overlayhdr=Rhdr
    elif projecttoband=='G':
        R_reproj=mcf.reproject2D(Rdata,Rhdr,Ghdr)
        G_reproj=Gdata.copy()
        B_reproj=mcf.reproject2D(Bdata,Bhdr,Ghdr)
        if overlayhdr==None: overlayhdr=Ghdr
    elif projecttoband=='B': 
        R_reproj=mcf.reproject2D(Rdata,Rhdr,Bhdr)
        G_reproj=mcf.reproject2D(Gdata,Ghdr,Bhdr)
        B_reproj=Bdata.copy()
        if overlayhdr==None: overlayhdr=Bhdr
    else: 
        R_reproj=Rdata.copy(); G_reproj=Gdata.copy(); B_reproj=Bdata.copy()
        if overlayhdr==None: overlayhdr=Bhdr
        
    if Rminmax==None: Rminmax=[np.nanmin(Rdata),np.nanmax(Rdata)]
    if Gminmax==None: Gminmax=[np.nanmin(Gdata),np.nanmax(Gdata)]
    if Bminmax==None: Bminmax=[np.nanmin(Bdata),np.nanmax(Bdata)]
    
    img = np.zeros((R_reproj.shape[0], R_reproj.shape[1], 4), dtype=float)
    img[:,:,3]+=alphachan #Alpha channel
    
    #R_scale=scalefuncs[Rscalefunc](R_reproj,scale_min=Rminmax[0],scale_max=Rminmax[1])
    #G_scale=scalefuncs[Gscalefunc](G_reproj,scale_min=Gminmax[0],scale_max=Gminmax[1])
    #B_scale=scalefuncs[Bscalefunc](B_reproj,scale_min=Bminmax[0],scale_max=Bminmax[1])
    R_scale=(mcf.scaling_fns[Rscalefunc]() + mcf.ManualInterval(vmin=Rminmax[0],vmax=Rminmax[1]))(R_reproj)
    G_scale=(mcf.scaling_fns[Gscalefunc]() + mcf.ManualInterval(vmin=Gminmax[0],vmax=Gminmax[1]))(G_reproj)
    B_scale=(mcf.scaling_fns[Bscalefunc]() + mcf.ManualInterval(vmin=Bminmax[0],vmax=Bminmax[1]))(B_reproj)
    
    img[:,:,0] = R_scale  #R   
    img[:,:,1] = G_scale  #G
    img[:,:,2] = B_scale  #B
    
    if ax!=None: ax[overlayhdr].imshow_affine(img,origin='lower') #pcolormesh can't natively do the RGB like imshow.
    else: plt.imshow(img,origin='lower',interpolation='nearest')
'''

def make_finder_plot_multicolor(targetname, coords, boxwidth, surveyscolors=[['DSS2 IR','#C11B17'], ['DSS2 Red','#FFFD74'], ['DSS2 Blue','#1569C7'],], boxwidth_units='pix', frame='icrs', overwrite=False, search_name=False, refregs=None, scalefuncs=['sqrt','sqrt','sqrt'], dpi=200, tickcolor='0.2', mfc='r', mec='w', bs_amin=0., surveytitle=False, savepathstem='./', filetype='jpg'):
    """
    Downloads cutout images (with astroquery.skyview, if not already on disk) of the specified survey bands, and creates a multicolor image finder plot.  Optional refregs can be specified to plot reference positions on the plot (notable guide stars, etc).
    
    Parameters
    ----------
    targetname : str
        The name of the target, to use in filenames and the plot title
    coords : list or str
        The list of target [RA,DEC] coordinates to use in astroquery.skyview.SkyView.get_images, 
        in decimal (not sexagesimal).  \n
        If search_name is True, this is a string name for the target source, to use with astroquery.
    boxwidth : int or float
        The width of the box sides of the desired image.  Units (e.g., pixels or arcminutes) 
        specified with boxwidth_units
    surveyscolors : list
        Each element in this list is a pair of a ['surveyname','#hexcolor'], a string 
        astroquery.skyview survey name, and a hex/html color to apply to that image.  
        (E.g., 'SDSSr', '2MASS-J', 'WISE 3.4' ...)  This list can include as many or as 
        few images as desired -- not limited to exactly 3 colors.
    boxwidth_units: str
        The units of the specified box width.  
        Options: 'pixels'/'pix', 'degree','deg', 'arcmin'/'amin', 'arcsec'/'asec'
    frame : str
        The coordinate reference frame.  Defaults to 'icrs'. Other options include 'galactic','fk5', ...
    overwrite : bool
        Set to True to overwrite the existing file at the specified savepath \n
        If you have deleted a fits file and the old one keeps reappearing, try deleting the 
        contents of your astroquery download cache, ~/.astopy/cache/download/url/ 
    search_name : bool
        Set to True to query SkyView by a common name (e.g., 'M33').  Specify this name with 
        the coords parameter instead of a list of [RA,DEC] coordinates.
    refregs : list
        refregs = list of reference position points in format [RA,DEC,'Name']
    scalefuncs : list
        List of the scaling functions from multicolorfits to use for scaling in the various images: 
        'linear', 'sqrt', 'log', 'asinh' 
    dpi : int
        The raster resolution (dots per inch) to apply to the saved image
    tickcolor : str or matplotlib color specifier 
        The color to use for the tick marks
    mfc : str or matplotlib color specifier
        Marker face color
    mec : str or matplotlib color specifier
        Marker edge color
    bs_amin : float
        The (angular) bar size to plot, in arcmin 
    surveytitle : bool
        Set to True to include the surveys used for the background image
    savepathstem : str
        The path string stem to use for saving the downloaded file.  
    filetype : str
        The filetype to use for the saved image.  'jpg','png','pdf','svg','eps'...
    
    Example
    -------
    refstars=[['5:34:42.3','22:10:34.5','HD 244988'], ['5:34:36.6','21:37:19.9','HD 36707'],]
    obs.make_finder_plot_multicolor('M1', obs.sex2dec('05:34:31.94','22:00:52.2'), 50., [['DSS2 Red','#A83C09'], ['DSS2 Blue','#336699']], boxwidth_units='arcmin', refregs=refstars, scalefuncs=['linear','linear'], dpi=200, tickcolor='0.8', mfc='#C11B17',mec='w',bs_amin=10., filetype='jpg')
    """
    rcparams_initial = {'xtick.major.size':rcParams['xtick.major.size'], 'ytick.major.size':rcParams['ytick.major.size'], 'xtick.minor.size':rcParams['xtick.minor.size'], 'ytick.minor.size':rcParams['ytick.minor.size'], 'grid.linestyle':rcParams['grid.linestyle'], 'grid.linewidth':rcParams['grid.linewidth'], 'xtick.direction':rcParams['xtick.direction'], 'ytick.direction':rcParams['ytick.direction'], 'xtick.minor.visible':rcParams['xtick.minor.visible'], 'ytick.minor.visible':rcParams['ytick.minor.visible']}
    rcparams_update={'xtick.major.size':7, 'ytick.major.size':7, 'xtick.minor.size':3, 'ytick.minor.size':3, 'grid.linestyle':':', 'grid.linewidth':0.3, 'xtick.direction':'in', 'ytick.direction':'in', 'xtick.minor.visible':True, 'ytick.minor.visible':True}
    rcParams.update(rcparams_update)
    
    fitssavestrings=[]
    #cutout_dats=[]; cutout_hdrs=[]
    images_colorized=[]
    
    for i in range(len(surveyscolors)):
        fitssavestring=savepathstem+targetname+'_'+surveyscolors[i][0]+'.fits'
        fitssavestrings.append(fitssavestring.replace(' ','_'))
        download_cutout(coords, boxwidth, fitssavestrings[i], survey=surveyscolors[i][0], boxwidth_units=boxwidth_units, frame=frame, overwrite=overwrite, search_name=search_name)
        cutout_dat,cutout_hdr=pyfits.getdata(fitssavestrings[i], header=True);
        
        if i==0: plothdr=cutout_hdr.copy()
        
        img_greyRGB=mcf.greyRGBize_image(cutout_dat,rescalefn=scalefuncs[i], gamma=2.2, checkscale=False)
        img_color=mcf.colorize_image(img_greyRGB, surveyscolors[i][1], colorintype='hex')
        images_colorized.append(img_color)
    
    img_multicolor=mcf.combine_multicolor(images_colorized,gamma=2.2)
    
    surveyliststring=str([sc[0] for sc in surveyscolors])
    
    fig1=plt.figure(1)#,figsize=(9,3.5))
    if surveytitle is True: 
        fig1.suptitle(targetname+'  %s'%(surveyliststring), fontsize=12, x=.5, y=.98) # single row
    else: fig1.suptitle(targetname, fontsize=12, x=.5, y=.98)
    cutout_wcs = WCS(plothdr)

    ax1=fig1.add_subplot(111, projection=cutout_wcs)
    ax1.imshow(img_multicolor,origin='lower',interpolation='nearest')
    rapars = ax1.coords[0]
    decpars = ax1.coords[1]
    rapars.set_ticks(color=tickcolor); decpars.set_ticks(color=tickcolor)
    rapars.set_ticks(number=6); decpars.set_ticks(number=6); 
    rapars.set_ticklabel(size=8); decpars.set_ticklabel(size=8); 
    #rapars.set_ticks(spacing=10*u.arcmin, color='white', exclude_overlapping=True)
    #decpars.set_ticks(spacing=5*u.arcmin, color='white', exclude_overlapping=True)
    rapars.display_minor_ticks(True); #rapars.set_minor_frequency(10)
    decpars.display_minor_ticks(True)
    rapars.set_major_formatter('hh:mm:ss.s')
    decpars.set_major_formatter('dd:mm')
    rapars.set_separator(('$^\mathrm{H}$', "'", '"'))
    #decpars.ticklabels.set_rotation(45) #Rotate ticklabels
    #decpars.ticklabels.set_color(xkcdrust) #Ticklabel Color
    
    #ax1.set_xlabel("$\Delta$ RA",size=8)
    #ax1.set_ylabel("$\Delta$ DEC",size=8)
    ax1.set_xlabel('RA (J2000)'); ax1.set_ylabel('DEC (J2000)');
    
    if type(coords) is str: tpx=[img_multicolor.shape[1]/2,cutoutR_dat.shape[0]/2]
    else: tpx=mcf.convsky2pix(plothdr,coords[0],coords[1],precise=True)
    ax1.plot(tpx[0],tpx[1],marker='+',mec=mfc,ms=8,zorder=10)
    
    if refregs!=None: 
        for reg in refregs:
            if type(reg[0])==str: regcoords=sex2dec(reg[0],reg[1])
            else: regcoords=reg[:2]
            rpx = bgpars2pix( [regcoords[0], regcoords[1], (boxwidth/80.)/60., (boxwidth/80.)/60., 0.], plothdr)
            c1 = patches.Ellipse( (rpx[0],rpx[1]), 2*rpx[2], 2*rpx[3], angle=rpx[4], linewidth=.8, ec=mfc, fill=False, zorder=5, path_effects=[PathEffects.withStroke(linewidth=1.2,foreground=mec)] )
            #plt.setp( c1, path_effects=[PathEffects.withStroke(linewidth=1,foreground=mec)] )
            ax1.add_patch(c1)
            ax1.text(rpx[0], rpx[1]+rpx[2]*1.1, reg[-1], va='bottom', ha='center', color=mfc, fontsize=8, path_effects=[PathEffects.withStroke(linewidth=.3,foreground=mec)], zorder=6 )
    
    ax1.grid(b=True, which='major', color='0.8', linestyle='-', linewidth=0.4)#,alpha=0.8)
    
    ### Size Bar
    if bs_amin>0:
        #bs_amin = Length of inset size bar, in arcmin
        barsizepix=bs_amin/(60.*mcf.getdegperpix(plothdr))
        dbar=add_sizebar(ax1, barsizepix, '%i $^{\prime}$'%(bs_amin), loc=3, sep=5, borderpad=0.3, frameon=False, path_effects=[PathEffects.withStroke(foreground=mec, linewidth=1.2)], color=mfc)
    
    figsavestring=savepathstem+targetname+'_multicolor'.replace(' ','')+'.'+filetype
    plt.savefig(figsavestring.replace(' ','_'),dpi=dpi,bbox_inches='tight')
    plt.clf(); plt.close()
    
    rcParams.update(rcparams_initial) #Revert back so as not to disrupt previous settings
    


###### ------------ Sensitivity / Integration Time Equations ------------ #######

def vlbi_image_sensitivity_homogeneous(SEFD_Jy,eta,Nstations,BW_total_Hz,t_int_sec,pol='dual'):
    """
    Calculate radio image sensitivity for arrays with homogeneous antennae, following: \n
    https://science.nrao.edu/facilities/vlba/docs/manuals/oss/imag-sens
    
    Parameters
    ----------
    SEFD_Jy : float
        System equivalent flux density, in units of Jansky
    eta : float 
        System inefficiency, in range [0,1]
    Nstations : int
        Number of homogeneous observing stations in the array
    BW_total_Hz : float
        Bandwidth in Hz in a single polarization
    t_int_sec : float
        Total integration time, in seconds
    pol : str
        Polarization descriptor.  ['single','r','l','x','y','1'] for single polarization, ['dual','full','rl','lr','xy','yx','2'] for dual polarization.
    
    Returns
    -------
    Sensitivity_total : float
        The estimated image sensitivity, in Jy/beam
    
    Notes
    -----
    For a single polarization: \n
    Sensitivity (Jy/beam) = SEFD / [ eta*(N*(N−1)*BW*t_int)^1/2 ] \n
    
    SEFD or "system equivalent flux density" is the system noise expressed in Janskys. 
        Values: https://science.nrao.edu/facilities/vlba/docs/manuals/oss/bands-perf#Table%205.1 \n
    eta <= 1  accounts for the VLBI system inefficiency (primarily quantization in the data recording).
        Walker (1995) and Kogan (1995b) provide the combination of scaling factors and 
        inefficiencies appropriate for VLBA visibility data. 
        The sensitivity values presented in table 5.1 were calculated using eta=0.8, while the EVN 
        Sensitivity Calculator uses eta=0.7. \n
    The bandwidth in Hz is BW (for a single polarization). 
        For a continuum target, use the channel bandwidth or the full recorded bandwidth, depending on the 
        fringe-fitting mode; for a line target, use the channel bandwidth divided by the number of spectral 
        points that span the channel. \n
    N is the number of observing stations  \n
    t_int is the total integration time on source in seconds.\n
      \n
    For simultaneous dual polarization, the image sensitivity for Stokes I, Q, U, or V is:\n
    Sensitivity_total = Sensitivity_singlePol/sqrt(2)
    
    
    Returns
    -------
    Sensitivity_total : float
        Estimated image sensitivity, in Jy/beam
    
    Examples
    --------
    #VLBA vc_b setting (6cm, SEFD=210 Jy) for 1 hour, dual pol. \n
    #(512MHz each pol), with NRAO eta = 0.8 \n
    #Sensitivity = 210. / (0.8*(10*(10-1)*512e6*3600)**.5) / np.sqrt(2) \n
    #            = 1.44e-5 = 14.4 uJy/bm \n
    obs.vlbi_image_sensitivity_homogeneous(210.,0.8,10,512e6,60*60,pol='dual') # 1.44e-5
    """
    
    Sensitivity_single = SEFD_Jy / ( eta*np.sqrt((Nstations*(Nstations-1)*BW_total_Hz*t_int_sec)) )
    
    if pol.lower() in ['single','r','l','x','y','1']:
        Sensitivity_total=Sensitivity_single
    elif pol.lower() in ['dual','full','rl','lr','xy','yx','2']:
        Sensitivity_total=Sensitivity_single/np.sqrt(2)
    else: 
        raise Exception('vlbi_image_sensitivity_homogeneous(): pol="%s" is invalid'%(pol))
    
    return Sensitivity_total #Image sensitivity in Jy/bm


def theoretical_vlbi_integration_time_homogeneous(Sensitivity_Jybm, SEFD_Jy, eta, Nstations, BW_total_Hz, pol='dual'):
    """
    Calculates the integration time required to achieve a given image sensitivity, for homogeneous antennae.
    Essentially inverting the equation in vlbi_image_sensitivity_homogeneous(): \n
    t_int_sec = (SEFD_Jy/(sqrt(#polarizations)*Sensitivity_Jybm*eta))^2 * 1./(N*(N-1)*BW_total_Hz)
    
    Parameters
    ----------
    Sensitivity_Jybm : float
        Image sensitivity, in Jy/beam
    SEFD_Jy : float
        System equivalent flux density, in units of Jansky
    eta : float 
        System inefficiency, in range [0,1]
    Nstations : int
        Number of homogeneous observing stations in the array
    BW_total_Hz : float
        Bandwidth in Hz in a single polarization
    pol : str
        Polarization descriptor.  ['single','r','l','x','y','1'] for single polarization, ['dual','full','rl','lr','xy','yx','2'] for dual polarization.
    
    Returns
    -------
    t_int_sec : float 
        The total integration time on source in seconds
    
    Example
    -------
    #To achieve 25 uJy/bm for VLBA obs at 5.0GHz (C-band,SEFD=210Jy), dual pol, 512MHz/pol, for eta=0.8: \n
    #t_int = (210/(np.sqrt(2)*25e-6*0.8))**2 / (10*(10-1)*512e6)  = 1196sec, or 0.33 hrs or 20 min \n
    obs.theoretical_vlbi_integration_time_homogeneous(25e-6,210,0.8,10,512e6,pol='dual') #--> 1196.29
    """
    if pol.lower() in ['single','r','l','x','y']: nr_pols=1
    elif pol.lower() in ['dual','full','rl','lr','xy','yx']: nr_pols=2
    else: raise Exception('theoretical_vlbi_integration_time_homogeneous(): pol="%s" is invalid'%(pol))
    t_int_sec = (SEFD_Jy/(Sensitivity_Jybm*eta))**2. * 1./(Nstations*(Nstations-1)*BW_total_Hz*nr_pols)
    return t_int_sec


def get_VLBA_sensitivity_table_info(query_freq_GHz):
    """
    Takes a nominal observation frequency in GHz and returns the relevant VLBA
    sensitivity calculation values from: \n
    https://science.nrao.edu/facilities/vlba/docs/manuals/oss/bands-perf \n
    (LBA Frequency Bands & Performance.  Table 5.1. )
    
    Parameters
    ----------
    query_freq_GHz : float
        The frequency of interest, in GHz.
    
    Returns
    -------
    VLBA_sensitivity_info : list
        [receiver_name, zenith_SEFD(Jy), peak_gain(K/Jy), baseline_sensitivity(60min), image_sensitivity(8hrs)]
    
    Example:
        obs.get_VLBA_sensitivity_table_info(5.)
        #Returns --> ['6 cm', 210, 0.119, 2.1, 5]
    """
    receiver=['90 cm', '50 cm', '21 cm', '18 cm ', '13 cm', '13 cm', '6 cm', '7 ghz', '4 cm', '4 cm', '2 cm', '1 cm', '24 ghz', '7 mm', '3 mm'] #receiver band designation
    freqrange_min=np.array([0.312, 0.596, 1.35, 1.35, 2.2, 2.2, 3.9, 3.9, 8.0, 8.0, 12.0, 21.7, 21.7, 41.0, 80.0]) #in GHz, Nominal frequency range edge minima
    freqrange_max=np.array([0.342, 0.626, 1.75, 1.75, 2.4, 2.4, 7.9, 7.9, 8.8, 8.8, 15.4, 24.1, 24.1, 45.0, 90.0]) #in GHz, Nominal frequency range edge maxima
    zenith_SEFD=np.array([2742, 2744, 289, 314, 347, 359, 210, 278, 327, 439, 543, 640, 534, 1181, 4236]) #in Jy, Typical Zenith SEFD
    SEFD_centerfreq=np.array([0.326, 0.611, 1.438, 1.658, 2.269, 2.269, 4.993, 6.66, 8.419, 8.419, 15.363, 22.236, 23.801, 43.124, 86.2]) #in GHz, Central Frequency for SEFD
    peak_gain=np.array([0.077, 0.078, 0.11, 0.112, 0.087, 0.085, 0.119, 0.103, 0.118, 0.105, 0.111, 0.11, 0.118, 0.09, 0.033]) #in K/Jy, Typical Peak Gain
    baseline_sensitivity=np.array([111, 443, 2.9, 3.2, 3.5, 3.6, 2.1, 2.8, 3.3, 4.4, 5.5, 6.5, 5.4, 12, 60]) #in mJy [512Mbps, 1 minute integration]
    image_sensitivity=np.array([266, 753, 10, 11, 12, 12, 5, 7, 8, 11, 13, 16, 13, 29, 184]) #in microJy/bm [4Gbps, 8 hour integration]
    
    try: ind=np.argwhere((query_freq_GHz>=freqrange_min)&(query_freq_GHz<freqrange_max))[0][0]
    except: raise Exception('VLBA_sensitivity_table_info(): %s GHz Not within table frequency range definitions'%(query_freq_GHz))
    
    return [receiver[ind], zenith_SEFD[ind], peak_gain[ind], baseline_sensitivity[ind], image_sensitivity[ind]]

def theoretical_VLBA_image_sensitivity(freq_GHz,BW_total_GHz,t_int_hr,eta=0.8,Nstations=10,pol='dual'):
    """
    Calculate the theoretical image sensitivity for VLBA observations of duration t_int_hr hours. 
    Calls vlbi_image_sensitivity_homogeneous() using values relevant to VLBA.
    
    Parameters
    ----------
    freq_GHz : float 
        The nominal frequency of the observations, in GHz. SEFD estimated from this value.
    BW_total_GHz : float 
        The total bandwidth of the receiver setup, in GHz, for a single polarization. 
        Default is 0.512 (= 512 MHz per polarization)
    t_int_hr : float  
        Integration time, in hours
    eta : float  
        The system (in)efficiency value, eta_s.  NRAO uses 0.8, the EVN calculator uses 0.7. 
    Nstations : int 
        Number of stations to use for the calculation. [default = 10, for the full array]
    pol : str
        Polarization descriptor.  ['single','r','l','x','y','1'] for single polarization, ['dual','full','rl','lr','xy','yx','2'] for dual polarization.
    
    Returns
    -------
    Sensitivity_Jybm : float
        Image sensitivity, in Jy/beam
    
    Example
    -------
    # 5.0 GHz obs in dual polarization, 512MHz BW per pol., 1 hr obs, eta=0.8, full array \n
    obs.theoretical_VLBA_image_sensitivity(5.0,0.512,1.0,eta=0.8,Nstations=10,pol='dual') \n
    # --> 1.4411421628371519e-05  or 14.41 microJy/beam
    """
    
    VLBA_SEFD_Jy=get_VLBA_sensitivity_table_info(freq_GHz)[1]
    #[receiver_name, zenith_SEFD(Jy), peak_gain(K/Jy), baseline_sensitivity(60min), image_sensitivity(8hrs)]
    
    Sensitivity_Jybm = vlbi_image_sensitivity_homogeneous(VLBA_SEFD_Jy, eta, Nstations, BW_total_GHz*1e9, t_int_hr*3600, pol=pol)
    
    return Sensitivity_Jybm

def theoretical_VLBA_image_integration_time(freq_GHz,BW_total_GHz,Sensitivity_uJybm, eta=0.8, Nstations=10, pol='dual', return_fmt='hours'):
    """
    Calculate the the theoretical integration time to achieve a certain sensitivity level with the VLBA.
    Calls theoretical_vlbi_integration_time_homogeneous()  using values relevant to the VLBA. 
    
    Parameters
    ----------
    freq_GHz : float 
        The nominal frequency of the observations, in GHz. SEFD estimated from this value.
    BW_total_GHz : float 
        The total bandwidth of the receiver setup, in GHz, for a single polarization. 
        Default is 0.512 (= 512 MHz per polarization)
    Sensitivity_uJybm : float 
        Desired image sensitivity level, in microJy/beam
    eta : float  
        The system (in)efficiency value, eta_s.  NRAO uses 0.8, the EVN calculator uses 0.7. 
    Nstations : int 
        Number of stations to use for the calculation. [default = 10, for the full array]
    pol : str
        Polarization descriptor.  ['single','r','l','x','y','1'] for single polarization, ['dual','full','rl','lr','xy','yx','2'] for dual polarization.
    return_fmt : str  ['sec', 'min', or 'hours']
        Format for the returned integration time
    
    Example
    -------
    #Goal = 25uJy/bm at 5.0 GHz in dual polarization, 512MHz BW per pol., eta=0.8, full array \n
    obs.theoretical_VLBA_image_integration_time(5.0, 0.512, 25.0, eta=0.8, Nstations=10, pol='dual', return_fmt='hours') \n
    # --> 0.3323025173611111  [hours]
    """
    
    VLBA_SEFD_Jy=get_VLBA_sensitivity_table_info(freq_GHz)[1]
    #[receiver_name, zenith_SEFD(Jy), peak_gain(K/Jy), baseline_sensitivity(60min), image_sensitivity(8hrs)]
    
    t_int_sec=theoretical_vlbi_integration_time_homogeneous(Sensitivity_uJybm*1e-6, VLBA_SEFD_Jy, eta, Nstations, BW_total_GHz*1e9, pol=pol)
    
    if return_fmt.lower() in ['hour','hours','hr','hrs']: t_int=t_int_sec/3600.
    elif return_fmt.lower() in ['minute','minutes','min','mins']: t_int=t_int_sec/60.
    elif return_fmt.lower() in ['second','seconds','sec','secs']: t_int=t_int_sec
    else: raise Exception('theoretical_VLBA_image_integration_time(): return_fmt="%s" is invalid. Use one of ["sec","min","hour"]'%(return_fmt))
    
    return t_int

def VLBA_freq_minimum_sun_separation(frequency_GHz):
    """
    Calculates the recommended minimum source separation (in degrees) from the   
    Sun for pointings in VLBA observations at a specified radio frequency in GHz.  \n
    The minimum separations vs frequency trend very cleanly follows a line in 
    log space, so arbitrary frequencies (beyond the exact defined values) are 
    calculated by interpolating from the defined NRAO list accordingly.
    
    Parameters
    ----------
    frequency_GHz : float
        The radio frequency of interest, in GHz
    
    Returns
    -------
    minsep : float
        The minimum separation from the Sun recommended for VLBA observations at the given frequency
        
    Note
    ----
    Following the definitions in Table 17.1:  \n
    https://science.nrao.edu/facilities/vlba/docs/manuals/oss/referencemanual-all-pages 
    
    Example
    -------
    obs.VLBA_freq_minimum_sun_separation(6.5)  #--> 19.74 (degrees minimum here at C band)
    """
    obsfreqs=[0.33,0.61,1.6,2.3,5.0,8.4,15,22,43] #GHz
    minsolseps=[117,81,45,36,23,17,12,9,6] #deg
    return 10**np.interp(np.log10(frequency_GHz),np.log10(obsfreqs),np.log10(minsolseps))




###### Misc Radio Tools #######

def band_from_freq(query_freq_GHz,print_table=False):
    """
    Takes an input frequency in GHz, and returns the radio band name string it corresponds to.
    
    Parameters
    ----------
    query_freq_GHz : float
        The radio frequency to query, in GigaHertz. Current definitions span [0.054,50.0] GHz (Band '4' through Q-band).
    print_table : bool
        Set to True to also print the entire table of radio band definitions (band name, frequency range, wavelength range) to screen
    
    Returns
    -------
    resultband : str
        The common name of the radio band corresponding to the input frequency
    
    Example
    -------
    obs.band_from_freq(7.3)  \n
    #-->  'C'
    """
    ### VLA Frequency Bands & Tunability.  Table 3.3.1.  
    # https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/vla-frequency-bands-and-tunability
    #Band names
    bands=['4','P', 'L', 'S', 'C', 'X', 'Ku', 'K', 'Ka', 'Q'] 
    #Nominal band edges (GHz for freqs, cm for wavs)
    freqrange_min=np.array([0.054, 0.20, 1.0, 2.0, 4.0, 8.0, 12.0, 18.0, 26.5, 40.0]) #in GHz
    freqrange_max=np.array([0.0842, 0.503, 2.04, 4.0, 8.0, 12.0, 18.0, 26.5, 40.0, 50.0]) #in GHz
    wavrange_min=299792458./(freqrange_max*1e9) *1e2 #in cm
    wavrange_max=299792458./(freqrange_min*1e9) *1e2 #in cm
    if print_table is True:
        print('Band Name   Frequency Range (GHz)   Wavelength Range (cm)\n%s'%('-'*60))
        for b in range(len(bands)):
            print('    %-2s         %7.3f - %-7.3f      %7.2f - %-7.2f'%(bands[b], freqrange_min[b],freqrange_max[b], wavrange_min[b],wavrange_max[b],))

    try: resultband=bands[np.argwhere((query_freq_GHz>=freqrange_min)&(query_freq_GHz<freqrange_max))[0][0]]
    except: resultband='%s GHz Not within NRAO band definitions'%(query_freq_GHz)

    return resultband

def info_from_idifits(idifits_path,print_style='short'):
    """
    Prints and returns info about the contents of an idifits file, in either long or short form.
    
    Parameters
    ----------
    idifits_path : str
        The path to the idifits file (relative & absolute paths OK)
    print_style : str  ('short' or 'long')
        Specifies the level of detail to print to screen/terminal
    
    Returns
    -------
    OBSCODE, DATE_OBS, SOURCE, sourcecoordlist, REF_FREQ, totalBW, NR_BANDS, BAND_FREQS : str, str, list, list, float, float, int, np.array
    
    Examples
    --------
    obs.info_from_idifits('myVLBAquasar_Kband.idifits',print_style='long'); \n
    obs.info_from_idifits('/home/observerdude/data/NGC3079_Kband.idifits',print_style='short'); \n
        \n
    for f in sorted(glob('/some/directory/with/idifits/')): \n
        obs.info_from_idifits(f,print_style='short');  \n
    #Note the semicolon at the end, to suppress the return output in ipython...
    """
    #pyfits.info(idifits_path)
    primhdr=pyfits.getheader(idifits_path,'PRIMARY') #Primary header
    srchdr=pyfits.getheader(idifits_path,'SOURCE') #Header with info about sources
    freqhdr=pyfits.getheader(idifits_path,'FREQUENCY') #Header with info about frequencies
    
    srctable=pyfits.getdata(idifits_path,'SOURCE') #Table of observed source info
    freqtable=pyfits.getdata(idifits_path,'FREQUENCY') #Table of frequency info
    
    totalbw=freqhdr['NO_BAND']*freqhdr['NO_CHAN']*freqhdr['CHAN_BW']/1e9
    sourcecoordlist = [ [ra,dec] for ra,dec in zip(srctable['RAOBS'],srctable['DECOBS']) ]
    
    if print_style=='long':
        ## Print overall info about the observing session & correlation
        print('\nOBSCODE  %s,  Observation Date  %s,  Correlation Date = %s\n'%(srchdr['OBSCODE'], primhdr['DATE-OBS'], primhdr['DATE-MAP']))
        ## Print info about each target source
        print('%sName  ID   RA_pointing          DEC_pointing'%(' '*11))
        try:
            for i in range(srchdr['NAXIS']):
                print('%15s  %2i  %18.13f  %18.13f'%(srctable['SOURCE'][i], srctable['SOURCE_ID'][i], srctable['RAOBS'][i], srctable['DECOBS'][i],))
        except: 
            print('%15s  %2i  %18.13f  %18.13f'%(srctable['SOURCE'], srctable['SOURCE_ID'], srctable['RAOBS'], srctable['DECOBS'],))
        #
        ## Print info about the frequency settings in the observations
        print('\n  Reference Frequency = %.3f GHz,  %i Bands,  %i Channels,  Channel Bandwidth = %.3f kHz,  Total BW = %.3f GHz'%(freqhdr['REF_FREQ']/1e9, freqhdr['NO_BAND'], freqhdr['NO_CHAN'], freqhdr['CHAN_BW']/1e3, freqhdr['NO_BAND']*freqhdr['NO_CHAN']*freqhdr['CHAN_BW']/1e9))
        for i in range(len(freqtable['FREQID'])):
            print('  FreqID  = %i'%(freqtable['FREQID'][i]))
            print('  Offset (MHz)    Channel Width (kHz)    Band Bandwidth (MHz)   Sideband   Baseband Channel #')
            for freq in range(len(freqtable['BANDFREQ'][i])):
                print('  %7.1f  %s%.1f  %s%.1f  %s%i %s%i'%(freqtable['BANDFREQ'][i][freq]/1e6, ' '*15,freqtable['CH_WIDTH'][i][freq]/1e3, ' '*15,freqtable['TOTAL_BANDWIDTH'][i][freq]/1e6, ' '*12,freqtable['SIDEBAND'][i][freq], ' '*12,freqtable['BB_CHAN'][i][freq]) )
        print('')
    elif print_style=='short':
        bandname=band_from_freq(freqhdr['REF_FREQ']/1e9,print_table=False)
        try: filename=idifits_path.split('/')[-1]
        except: filename=idifits_path
        print('Obscode = %s,  Targets = %s,  Band = %s,  Ref.Freq = %.3f GHz,  filename = %s'%(srchdr['OBSCODE'], str(srctable['SOURCE']), bandname, freqhdr['REF_FREQ']/1e9, filename))
    else: raise Exception('info_from_idifits(): print_style = "%s" not valid.  Must be "short" or "long"'%(print_style))
    
    ### Return header information:
    #   Obscode, Obsdate, Source list, List of Coords (RA/DEC), Ref.Freq.(GHz), TotalBandwidth (GHz), # Bands, Band Offsets (MHz)
    return srchdr['OBSCODE'], primhdr['DATE-OBS'], srctable['SOURCE'], sourcecoordlist, freqhdr['REF_FREQ']/1e9, totalbw, freqhdr['NO_BAND'], freqtable['BANDFREQ']



###### VLBA stations #######

### Create the VLBA telescope ephem.Observer objects -- e.g. Brewster, Pietown, ... 
#'BR','FD','HN','KP','LA','MK','NL','OV','PT','SC'

#vlbanames=np.array(['BR','FD','HN','KP','LA','MK','NL','OV','PT','SC'], dtype='<U2')
#eastlongs_dec=np.array([240.31672222,256.05518333,288.01341944,248.38757778, 253.75440278,204.54449722,268.42586667,241.72295278, 251.88081667,295.41636944])
#northlats_dec=np.array([48.13122778,30.63503056,42.93360833,31.95630556,35.775125,19.80138056,41.771425,37.23165278,34.30100278, 17.75657778])
#elevations_m=np.array([250,1606,296,1902,1962,3763,222,1196,2365,-15])
#ut_offsets=np.array([-7,-6,-4,-7,-7,-10,-6,-7,-7,-4])

### To calculate precise station positions from the TRF solutions positions (WGS84/ECEF XYZ Cartesian coords):
#   Brewster XYZ coords: X=-2112065168.35, Y=-3705356499.15,  Z=4726813692.71
#       --> values here from usn2020d global solution TRF .sta , and these values are in millimeters
#   Convert using the Vincenty formulae for the WGS-84 ellipsoidal earth model
#       e.g., https://mrjean1.github.io/PyGeodesy/docs/pygeodesy.ellipsoidalVincenty-module.html
#       pygeodesy.ellipsoidalVincenty.Cartesian expects coords in meters, so need to divide TRF mm coords by 1e3
#   
#   from pygeodesy.ellipsoidalVincenty import LatLon, Cartesian
#   latlon_N_W_BR=Cartesian(-2112065168.35/1e3,-3705356499.15/1e3,4726813692.71/1e3).toLatLon() #lat/lon in +N,+W
#   print(latlon_N_W_BR)     # 48.131224°N, 119.68328°W, +250.48m  --> May print as +W, but actually stored as +E
#   print(latlon_N_W_BR.lon, wrap_360(latlon_N_W_BR.lon)) # -119.6832798093127, 240.3167201906873
#   #pyephem takes lat/lon coords either in radians (floats), or as strings (D:M:S, delimited by colons)  
#   #To format into strings (also useful for printing), can do this: 
#   lon_lat_BR=dec2sex(wrap_360(latlon_N_W_BR.lon),latlon_N_W_BR.lat,as_string=True,decimal_places=10,RA=False)
#   height_BR=str(latlon_N_W_BR.height)
#   print(lon_lat_BR+[height_BR,]) #['240:19:00.1926864743', '48:7:52.4051932358','250.4771973239258']
#   #Note that input has 13 decimal places --> in principle could print string out to 1e-13/3600 ~ 1e-17 or 17 decimal places...

#vlbaBR=ephem.Observer() #Must be initialized with blank parentheses
##vlbaBR.lat=48.1312*np.pi/180; vlbaBR.lon=240.3167*np.pi/180; vlbaBR.elevation=250.  
##===> NOTE that pyephem expects decimal coords to be in RADIANS!!
## pyephem longitudes are +E, or (increasing) East coords (pygeodesy may print coords as +West, if the +E are negative)
#vlbaBR.lon='240:19:00.20'; vlbaBR.lat='48:07:52.42'; vlbaBR.elevation=250; vlbaBR.name='vlbaBR'

#vlbaFD=ephem.Observer(); vlbaFD.lon='256:03:18.66'; vlbaFD.lat='30:38:06.11'; vlbaFD.elevation=1606; vlbaFD.name='vlbaFD'
#vlbaHN=ephem.Observer(); vlbaHN.lon='288:00:48.31'; vlbaHN.lat='42:56:00.99'; vlbaHN.elevation=296;  vlbaHN.name='vlbaHN'
#vlbaKP=ephem.Observer(); vlbaKP.lon='248:23:15.28'; vlbaKP.lat='31:57:22.70'; vlbaKP.elevation=1902; vlbaKP.name='vlbaKP'
#vlbaLA=ephem.Observer(); vlbaLA.lon='253:45:15.85'; vlbaLA.lat='35:46:30.45'; vlbaLA.elevation=1962; vlbaLA.name='vlbaLA'
#vlbaMK=ephem.Observer(); vlbaMK.lon='204:32:40.19'; vlbaMK.lat='19:48:04.97'; vlbaMK.elevation=3763; vlbaMK.name='vlbaMK'
#vlbaNL=ephem.Observer(); vlbaNL.lon='268:25:33.12'; vlbaNL.lat='41:46:17.13'; vlbaNL.elevation=222;  vlbaNL.name='vlbaNL'
#vlbaOV=ephem.Observer(); vlbaOV.lon='241:43:22.63'; vlbaOV.lat='37:13:53.95'; vlbaOV.elevation=1196; vlbaOV.name='vlbaOV'
#vlbaPT=ephem.Observer(); vlbaPT.lon='251:52:50.94'; vlbaPT.lat='34:18:03.61'; vlbaPT.elevation=2365; vlbaPT.name='vlbaPT'
#vlbaSC=ephem.Observer(); vlbaSC.lon='295:24:58.93'; vlbaSC.lat='17:45:23.68'; vlbaSC.elevation=-15;  vlbaSC.name='vlbaSC'

vlbaBR=create_ephem_observer('vlbaBR', '240:19:00.20', '48:07:52.42',  250, timezone='America/Los_Angeles')
vlbaFD=create_ephem_observer('vlbaFD', '256:03:18.66', '30:38:06.11', 1606, timezone='America/Chicago')
vlbaHN=create_ephem_observer('vlbaHN', '288:00:48.31', '42:56:00.99',  296, timezone='America/New_York')
vlbaKP=create_ephem_observer('vlbaKP', '248:23:15.28', '31:57:22.70', 1902, timezone='America/Phoenix')
vlbaLA=create_ephem_observer('vlbaLA', '253:45:15.85', '35:46:30.45', 1962, timezone='America/Denver')
vlbaMK=create_ephem_observer('vlbaMK', '204:32:40.19', '19:48:04.97', 3763, timezone='Pacific/Honolulu')
vlbaNL=create_ephem_observer('vlbaNL', '268:25:33.12', '41:46:17.13',  222, timezone='America/Chicago')
vlbaOV=create_ephem_observer('vlbaOV', '241:43:22.63', '37:13:53.95', 1196, timezone='America/Los_Angeles')
vlbaPT=create_ephem_observer('vlbaPT', '251:52:50.94', '34:18:03.61', 2365, timezone='America/Denver')
vlbaSC=create_ephem_observer('vlbaSC', '295:24:58.93', '17:45:23.68',  -15, timezone='America/St_Thomas')

#Observer.date defaults to current time when computer creates the object... need to update it before making any calculations!



#################################################################

#### Separations from standard fringe finders/BPcals  (also list some flux calibrators)
# Source positions from  /home/pjcigan/programs/sched/sched_11.6/catalogs/sources.gsfc2016
#SOURCE='3C84','J0319+4130','0316+413','J0319+41'       RA= 03:19:48.1600956 DEC= +41:30:42.104043
#SOURCE='0552+398','J0555+3948','DA193','J0555+39'      RA= 05:55:30.8056108 DEC= +39:48:49.164962
#SOURCE='4C39.25','J0927+3902','0923+392','J0927+39'    RA= 09:27:03.0139367 DEC= +39:02:20.851846
#SOURCE='3C273B','J1229+0203','3C273','1226+023','J1229+02'     RA= 12:29:06.6997310 DEC= +02:03:08.598201
#SOURCE='3C345','J1642+3948','1641+399','J1642+39'      RA= 16:42:58.8099666 DEC= +39:48:36.993986
#SOURCE='1921-293','J1924-2914','J1924-29'              RA= 19:24:51.0559525 DEC= -29:14:30.121070
#SOURCE='3C454.3','J2253+1608','2251+158','J2253+16'    RA= 22:53:57.7479385 DEC= +16:08:53.560908
#SOURCE='0234+285','J0237+2848','J0237+28'              RA= 02:37:52.4056765 DEC= +28:48:08.989988
#SOURCE='0528+134','J0530+1331','J0530+13'              RA= 05:30:56.4167463 DEC= +13:31:55.149483
#SOURCE='1758+388','J1800+3848','J1800+38'              RA= 18:00:24.7653628 DEC= +38:48:30.697484
#SOURCE='2007+777','J2005+7752','J2005+77'              RA= 20:05:30.9985258 DEC= +77:52:43.247541
#Flux calibrators: 3C84, 3C138, 3C147, 3C196, 3C286, 3C295
#3C84 (see above)
#SOURCE='3C138','J0521+1638','0518+165','J0521+16'      RA= 05:21:09.8859649 DEC= +16:38:22.051601
#SOURCE='3C147','J0542+4951','0538+498','J0542+49'      RA= 05:42:36.1378941 DEC= +49:51:07.233721
#3C196: no entry in sources.gsfc2016...  From NED: RA= 08:13:36.033, DEC= 48:13:02.56
#SOURCE='3C286','J1331+3030','1328+307','J1331+30'      RA= 13:31:08.2880613 DEC= +30:30:32.959300
#3C295: no entry in sources.gsfc2016...  From NED: RA= 14:11:20.519, DEC=52:12:09.97


### Fringe-finders
SRC_3C84=create_ephem_target('3C84','03:19:48.1600956','41:30:42.104043') 
SRC_DA193=create_ephem_target('DA193','05:55:30.8056108','39:48:49.164962') 
SRC_4C39p25=create_ephem_target('4C39.25','09:27:03.0139367','39:02:20.851846') 
SRC_3C273=create_ephem_target('3C273','12:29:06.6997310','02:03:08.598201') 
SRC_3C345=create_ephem_target('3C345','16:42:58.8099666','39:48:36.993986') 
SRC_1921m293=create_ephem_target('1921-293','19:24:51.0559525','-29:14:30.121070') 
SRC_3C454p3=create_ephem_target('3C454.3','22:53:57.7479385','16:08:53.560908') 
SRC_0234p285=create_ephem_target('0234+285','02:37:52.4056765','28:48:08.989988') 
SRC_0528p134=create_ephem_target('0528+134','05:30:56.4167463','13:31:55.149483') 
SRC_J1800p3848=create_ephem_target('J1800+3848','18:00:24.7653628','38:48:30.697484') #'1758+388'
SRC_2007p777=create_ephem_target('2007+777','20:05:30.9985258','77:52:43.247541') 
### Flux Calibrators
#SRC_3C84, defined above
SRC_3C138=create_ephem_target('3C138','05:21:09.8859649','16:38:22.051601') 
SRC_3C147=create_ephem_target('3C147','05:42:36.1378941','49:51:07.233721') 
SRC_3C196=create_ephem_target('3C196','08:13:36.033','48:13:02.56') 
SRC_3C286=create_ephem_target('3C286','13:31:08.2880613','30:30:32.959300') 
SRC_3C295=create_ephem_target('3C295','14:11:20.519','52:12:09.97') 




def nearest_from_target_list(obstarget,reference_source_list,verbose=False):
    """
    For a given observation target, determine the nearest from among a list of reference sky targets.
    
    Parameters
    ----------
    obstarget : ephem.FixedBody()
        The sky source of interest. 
    reference_source_list : array-like [list, tuple, np.array...]
        The list of reference objects, in ephem.FixedBody() format
    verbose : bool
        Set to True to print out the distance to each source in reference_source_list
    
    Returns
    -------
    src_nearest : str
        The name (ephem.FixedBody().name string) of the nearest Fringe Finder to the target
    
    Examples
    --------
    obs.nearest_from_target_list(obstarget,[obs.SRC_3C84,obs.SRC_3C286,obs.SRC_3C273],verbose=True)
    """
    angseps=[]
    for refsrc in reference_source_list:
        angseps.append(skysep_fixed_single(obstarget,refsrc))
    if verbose==True: 
        print('\nAngular separations on sky from %s:'%(obstarget.name))
        for i in range(len(reference_source_list)): 
            print('  %10s = %.2f deg'%(reference_source_list[i].name,angseps[i]))
    return reference_source_list[np.nanargmin(angseps)].name
    

def nearest_standard_fringefinder(obstarget,verbose=False):
    """
    For a given observation target, determine the nearest fringe finder, and optionally print the distance to standard fringe finders. 
    
    Parameters
    ----------
    obstarget : ephem.FixedBody()
        The sky source of interest. 
    verbose : bool
        Set to True to print out the distance to each fringe finder
    
    Returns
    -------
    src_nearest : str
        The name (ephem.FixedBody().name string) of the nearest Fringe Finder to the target
    
    Examples
    --------
    ngc3079=obs.create_ephem_target('NGC3079','10:01:57.80','55:40:47.24') \n
    obs.nearest_standard_fringefinder(ngc3079)   #--> '4C39.25'  \n
    obs.nearest_standard_fringefinder(ngc3079,verbose=True)   #--> prints individual seps \n
    ### \n
    # You can do this for the Sun/Moon as well, which first requires coords at some time \n
    suncoords=obs.compute_target_altaz_single(ephem.Sun(),obs.vlbaBR,'2021/01/01 00:00:00') \n
    obs.nearest_standard_fringefinder(suncoords)
    """
    
    #ffsrcs=[SRC_3C84,SRC_DA193,SRC_4C39p25,SRC_3C273,SRC_3C345,SRC_1921m293,SRC_3C454p3]
    #ffseps=[]
    #for ffsrc in ffsrcs:
    #    ffseps.append(skysep_fixed_single(obstarget,ffsrc))
    #if verbose==True: 
    #    print('\nAngular separations on sky from %s:'%(obstarget.name))
    #    for i in range(len(ffsrcs)): print('  %10s = %.2f deg'%(ffsrcs[i].name,ffseps[i]))
    #return ffsrcs[np.nanargmin(ffseps)].name
    
    src_nearest = nearest_from_target_list(obstarget, [SRC_3C84, SRC_DA193, SRC_4C39p25, SRC_3C273, SRC_3C345, SRC_1921m293, SRC_3C454p3, SRC_0234p285, SRC_0528p134, SRC_J1800p3848, SRC_2007p777], verbose=verbose)
    
    return src_nearest




### Todo:
# * bundle all of the time conversion functions into a single class that can just handle input timezones and do the conversions internally as needed.  [Probably convert to and store as pyephem.Date internally using any supplied tz else assume UTC, then for any time a local aware clock time is needed just return a utc_to_local datetime]


### For interferometer sensitivity equation ( VLBA-specific version)
#   --> Also want to make a function where you can just specify the standard VLBA receiver observing modes? 
#       (would need to hard-code those from the sched file, then call with a dict or something...)
#       That would put us well on the way to having our own python version of the sensitivity calculator...




