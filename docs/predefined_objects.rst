
Pre-defined objects
===================

Some telescope sites (pyephem.Observer objects) and astronomical sources (pyephem.FixedBody objects) are defined and callable within obsplanning. This list may grow in the future.


Telescopes
~~~~~~~~~~

.. code:: python

    obs.vlbaBR  # VLBA Brewster
    obs.vlbaFD  # VLBA Fort Davis
    obs.vlbaHN  # VLBA Hancock
    obs.vlbaKP  # VLBA Kitt Peak
    obs.vlbaLA  # VLBA Los Alamos
    obs.vlbaMK  # VLBA Mauna Kea
    obs.vlbaNL  # VLBA North Liberty
    obs.vlbaOV  # VLBA Owens Valley
    obs.vlbaPT  # VLBA Pietown
    obs.vlbaSC  # VLBA St Croix




Sources
~~~~~~~

.. code:: python

    ### Fringe finders
    obs.SRC_3C84        # 3C84'
    obs.SRC_DA193       # DA193
    obs.SRC_4C39p25     # 4C39.25
    obs.SRC_3C273       # 3C273
    obs.SRC_3C345       # 3C345
    obs.SRC_1921m293    # 1921-293
    obs.SRC_3C454p3     # 3C454.3
    obs.SRC_0234p285    # 0234+285
    obs.SRC_0528p134    # 0528+134
    obs.SRC_J1800p3848  # J1800+3848
    obs.SRC_2007p777    # 2007+777
    
    ### Flux Calibrators
    #[SRC_3C84, defined above]
    obs.SRC_3C138   # 3C138
    obs.SRC_3C147   # 3C147
    obs.SRC_3C196   # 3C19
    obs.SRC_3C286   # 3C286
    obs.SRC_3C295   # 3C295

