obsplanning -- a set of python utilities to aid in planning astronomical observations

version 1.0.1

API documentation at [https://obsplanning.readthedocs.io](https://obsplanning.readthedocs.io)





[![PyPI version](https://badge.fury.io/py/obsplanning.svg)](https://badge.fury.io/py/obsplanning)
[![Downloads](https://pepy.tech/badge/obsplanning)](https://pepy.tech/project/obsplanning)


Sharing/Customization: Please, play around!  (MIT License)

If you find this useful for your work, giving me (Phil Cigan) a nod in your acknowledgements would be greatly appreciated.  





-----------------------

# Elevator Pitch

Tools to plan your astronomical observations, from simple conversions of coordinates, to handlnig date/time objects, to tracking targets across the sky, to plotting functions and visualizations, to radio sensitivity calculations.  There's a little bit of everything.  I hope you find it useful, o intrepid observer.



-----------------------

# Dependencies

* numpy
* matplotlib 
* datetime
* pytz
* tzwhere
* astropy
* ephem
* astroquery
* scipy
* tqdm
* multicolorfits


- Tested in python 3.7 



# Installation

Install with pip
```console
pip install obsplanning
```

Alternatively, you can simply save a copy of obsplanning.py in a local working directory.  This is suitable for situations when a regular pip install is not feasible, or for simply testing it out as a standalone file of functions (running/loading the file within a script to access its functionality).


# Usage

Basic usage:
```python
import obsplanning as obs

#e.g., computing target altitude & azimuth from an observer station:
obs.compute_target_altaz(target,observer,t1,t2,nsteps=1000)

#etc...
```

See the examples in the following section, and the documentation at [https://obsplanning.readthedocs.io](https://obsplanning.readthedocs.io) for much more detail.


-----------------------

# Tutorials / Examples 

See the documentation at [readthedocs](https://obsplanning.readthedocs.io).  
Or, see examples on github [here](./examples.md).



# Features

Determine observability of objects in the sky from your observatory, and produce plots to help you prepare for your observations over the course of a session.  

- pyephem-based time conversions, observers (telescope), and target source (sky) objects
- Plotting/visualization tools for target elevation tracks and finder plots
- Transit time and target visibility calculations for single telescopes or multiple stations (VLBI/VLBA)
- Simple selection of observable targets from source lists, based on their coordinates. 
- Miscellaneous coordinate handling, radio astronomy info, and other helpful functions




-----------------------

# Motivation / Backstory

This collection of functions assembled somewhat organically over the years from a variety of different clusters of code I've written for different purposes.  (Early versions of some go back to around 2011!)  As I started getting more and more hands-on with calculations and plots for my own observations, I started writing more and more tools to make those processes more convenient.  There are certainly some very nice tools out in the world (like the venerable staralt), but I wanted some python tools to fit better with my own workflow.  Circa 2017 or so when I started getting serious about planning some observations at the WHT within my own ecosystem of code, I started a simple text file where I dumped all my general functions related to this.  None of the other packages I saw at that time really offered everything I was looking for in one place, so here we are.  As time went by, I added more bits and pieces as they came up, and now there are even some radio interferometry tools as well.  The current form here is not really object-oriented at all, but rather a collection of functions -- in the future I hope to have enough time and bandwidth to refactor these into classes and transform it to a more proper package.  Astropy and other resources definitely do have good tools for planning your observations now, but hopefully the functionality here will at least be complementary. 

Oh, and by the way - once you take your fancy observations with the help of obsplanning, how about making some fancy images using [multicolorfits](https://github.com/pjcigan/multicolorfits)?  Flag me down at a conference or wherever, I'd love to see what works of art you will make.  (Astronomy does have the best pictures in science, after all...)

Clear skies











