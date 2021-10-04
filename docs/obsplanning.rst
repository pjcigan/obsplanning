#####################
obsplanning 
#####################



=====================
Installation
=====================

* With pip:

.. code-block:: shell

    pip install obsplanning


* Standalone:

Alternatively, can download just the obsplanning.py to run in standalone mode -- this may be a better option if pip installation is not possible on your machine.  It's possible to simply run/execfile obsplanning.py in your python session to load its functions.  You could also add the file location to your local path:

(for .bashrc)

.. code-block:: shell

    export PYTHONPATH=$PYTHONPATH:/path/to/dir/with/obs/

(for .cshrc)

.. code-block:: shell

    setenv PYTHONPATH ${PYTHONPATH}:/path/to/dir/with/obs/


=====================
Basic Usage
=====================

.. code-block:: python

    import obsplanning as obs
    #e.g., computing target altitude & azimuth from an observer station:
    obs.compute_target_altaz(target,observer,t1,t2,nsteps=1000)  

Or, as a standalone script, download and cd to that directory, and type: 

.. code-block:: shell
    
    python obsplanning.py


=====================
Examples
=====================

See the main examples page on Github:

 `https://github.com/pjcigan/obsplanning/blob/master/examples.md <https://github.com/pjcigan/obsplanning/blob/master/examples.md>`_



