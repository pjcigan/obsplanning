
Pre-defined objects
===================

Some telescope sites (pyephem.Observer objects) and astronomical sources (pyephem.FixedBody objects) are defined and callable within obsplanning. This list may grow in the future.


Telescopes
~~~~~~~~~~

Optical / Infrared
------------------

.. code:: python
    
    obs.SIT_AAO       # Anglo-Australian Observatory, Siding Spring
    obs.SIT_CFHT      # Canada-France-Hawaii Telescope
    obs.SIT_CTIO      # Cerro Tololo Inter-American Observatory
    obs.SIT_GeminiN   # Gemini North, Hawaii
    obs.SIT_GeminiS   # Gemini South, Chile
    obs.SIT_GMT       # Giant Magellan Telescope, Chile
    obs.SIT_GTC       # Gran Telescopio Canarias, La Palma
    obs.SIT_HET       # Hobby-Eberly Telescope
    obs.SIT_Keck      # Keck Observatory, Hawaii
    obs.SIT_LDT       # Lowell Discovery Telescope (Prev. DCT)
    obs.SIT_Lick      # Lick Observatory, San Jose
    obs.SIT_Magellan  # Magellan twin 6.5 meters at LCO
    obs.SIT_Meudon    # Meudon great refractor
    obs.SIT_OHP       # Observatoire de Haute Provence
    obs.SIT_Palomar   # Palomar Observatory, San Diego
    obs.SIT_SALT      # South African Large Telescope
    obs.SIT_Sloan     # of SDSS fame, in Apache Point Obs. in Sunspot, NM
    obs.SIT_Subaru    # Subaru 8.2m on Mauna Kea, Hawaii
    obs.SIT_UKIRT     # UK Infrared Telescope, Mauna Kea, Hawaii
    obs.SIT_VLT       # Very Large Telescope, at ESO Cerro Paranal, Chile
    obs.SIT_WHT       # William Herschel Telescope, La Palma
    obs.SIT_WIYN      # WIYN 3.5m on Kitt Peak, Arizona
    obs.SIT_Yerkes    # Yerkes Observatory, Williams Bay, WI


Radio / submm
-------------

.. code:: python
    
    obs.SIT_ALMA        # Atacama Large Millimetre/submillimetre Array, Chile
    obs.SIT_Arecibo     # Arecibo Observatory, Puerto Rico
    obs.SIT_ASKAP       # Australian SKA Pathfinder
    obs.SIT_ATCA        # Australian Telescope Compact Array, Narrabri
    obs.SIT_Effelsberg  # Effelsberg 100m, Germany
    obs.SIT_FAST        # Five-hundred-meter Aperture Spherical Telescope, China
    obs.SIT_GBT         # Green Bank Telescope, West Virginia
    obs.SIT_GMRT        # Giant Metrewave Radio Telescope, Pune India
    obs.SIT_HARTRAO     # Hartebeesthoek Radio Astronomy Observatory
    obs.SIT_Haystack    # MIT Haystack, Massachussetts
    obs.SIT_JCMT        # James Clerk Maxwell Telescope, on Mauna Kea, Hawaii
    obs.SIT_Jodrell     # Jodrell Bank Observatory, England
    obs.SIT_Nancay      # Nancay Radio Observatory, France
    obs.SIT_NOEMA       # Northern Extended Millimeter Array (prev.Plateau de Bure), France 
    obs.SIT_Parkes      # CSIRO Parkes Observatory, Australia
    obs.SIT_VLA         # NRAO Karl G. Jansky Very Large Array, New Mexico
    obs.SIT_WSRT        # Westerbork Synthesis Radio Telescope, Netherlands
        
    ### VLBA
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
    
    ### IVS/VLBI stations, with coordinates from USNO quarterly solution
    obs.SIT_AGGO,     obs.SIT_AIRA,     obs.SIT_ALGOPARK, obs.SIT_AUSTINTX, 
    obs.SIT_AZORES,   obs.SIT_BADARY,   obs.SIT_BLKBUTTE, obs.SIT_BLOOMIND, 
    obs.SIT_BR_VLBA,  obs.SIT_BREST,    obs.SIT_CARNUSTY, obs.SIT_CARROLGA, 
    obs.SIT_CHICHI10, obs.SIT_CHLBOLTN, obs.SIT_CRIMEA,   obs.SIT_CTVASBAY, 
    obs.SIT_CTVASTJ,  obs.SIT_DEADMANL, obs.SIT_DSS13,    obs.SIT_DSS15, 
    obs.SIT_DSS26,    obs.SIT_DSS34,    obs.SIT_DSS36,    obs.SIT_DSS45, 
    obs.SIT_DSS56,    obs.SIT_DSS65,    obs.SIT_DSS65A,   obs.SIT_EFLSBERG, 
    obs.SIT_ELY,      obs.SIT_FD_VLBA,  obs.SIT_FLAGSTAF, obs.SIT_FORT_ORD, 
    obs.SIT_FORTLEZA, obs.SIT_FORTORDS, obs.SIT_FTD_7900, obs.SIT_GGAO12M,  
    obs.SIT_GGAO7108, obs.SIT_GILCREEK, obs.SIT_GOLDMARS, obs.SIT_GOLDVENU, 
    obs.SIT_GORF7102, obs.SIT_GRASSE,   obs.SIT_HALEAKAL, obs.SIT_HART15M,  
    obs.SIT_HARTRAO,  obs.SIT_HATCREEK, obs.SIT_HAYSTACK, obs.SIT_HN_VLBA,  
    obs.SIT_HOBART12, obs.SIT_HOBART26, obs.SIT_HOFN,     obs.SIT_HOHNBERG, 
    obs.SIT_HRAS_085, obs.SIT_ISHIOKA,  obs.SIT_JPL_MV1,  obs.SIT_KARLBURG, 
    obs.SIT_KASHIM11, obs.SIT_KASHIM34, obs.SIT_KASHIMA,  obs.SIT_KATH12M,  
    obs.SIT_KAUAI,    obs.SIT_KIRSBERG, obs.SIT_KODIAK,   obs.SIT_KOGANEI,  
    obs.SIT_KOKEE,    obs.SIT_KOKEE12M, obs.SIT_KP_VLBA,  obs.SIT_KUNMING,  
    obs.SIT_KWAJAL26, obs.SIT_LA_VLBA,  obs.SIT_LEONRDOK, obs.SIT_MACGO12M, 
    obs.SIT_MADRID64, obs.SIT_MAMMOTHL, obs.SIT_MARCUS,   obs.SIT_MARPOINT, 
    obs.SIT_MATERA,   obs.SIT_MCD_7850, obs.SIT_MEDICINA, obs.SIT_METSAHOV, 
    obs.SIT_METSHOVI, obs.SIT_MIAMI20,  obs.SIT_MILESMON, obs.SIT_MIZNAO10, 
    obs.SIT_MK_VLBA,  obs.SIT_MOJ_7288, obs.SIT_MOJAVE12, obs.SIT_MON_PEAK, 
    obs.SIT_NL_VLBA,  obs.SIT_NOBEY_6M, obs.SIT_NOME,     obs.SIT_NOTO,     
    obs.SIT_NRAO_140, obs.SIT_NRAO20,   obs.SIT_NRAO85_3, obs.SIT_NYALE13S, 
    obs.SIT_NYALES20, obs.SIT_OCOTILLO, obs.SIT_OHIGGINS, obs.SIT_ONSA13SW, 
    obs.SIT_ONSALA60, obs.SIT_OV_VLBA,  obs.SIT_OVR_7853, obs.SIT_OVRO_130, 
    obs.SIT_PARKES,   obs.SIT_PBLOSSOM, obs.SIT_PENTICTN, obs.SIT_PIETOWN,  
    obs.SIT_PINFLATS, obs.SIT_PLATTVIL, obs.SIT_PRESIDIO, obs.SIT_PT_REYES, 
    obs.SIT_PVERDES,  obs.SIT_QUINCY,   obs.SIT_RAEGSMAR, obs.SIT_RAEGYEB,  
    obs.SIT_RICHMOND, obs.SIT_ROBLED32, obs.SIT_SANPAULA, obs.SIT_SANTIA12, 
    obs.SIT_SC_VLBA,  obs.SIT_SEATTLE1, obs.SIT_SEJONG,   obs.SIT_SESHAN25, 
    obs.SIT_SINTOTU3, obs.SIT_SNDPOINT, obs.SIT_SOURDOGH, obs.SIT_SVETLOE,  
    obs.SIT_SYOWA,    obs.SIT_TIANMA65, obs.SIT_TIDBIN64, obs.SIT_TIGOCONC, 
    obs.SIT_TIGOWTZL, obs.SIT_TOULOUSE, obs.SIT_TROMSONO, obs.SIT_TRYSILNO, 
    obs.SIT_TSUKUB32, obs.SIT_URUMQI,   obs.SIT_VERAISGK, obs.SIT_VERAMZSW, 
    obs.SIT_VERNAL,   obs.SIT_VICTORIA, obs.SIT_VNDNBERG, obs.SIT_WARK12M,  
    obs.SIT_WESTFORD, obs.SIT_WETTZ13N, obs.SIT_WETTZ13S, obs.SIT_WETTZELL, 
    obs.SIT_WHTHORSE, obs.SIT_YAKATAGA, obs.SIT_YARRA12M, obs.SIT_YEBES,    
    obs.SIT_YEBES40M, obs.SIT_YELLOWKN, obs.SIT_YLOW7296, obs.SIT_YUMA,     
    obs.SIT_ZELENCHK, 


Other
-----

.. code:: python

    obs.SIT_ASTRON      # Netherlands Institute for Radio Astronomy
    obs.SIT_Goddard     # GSFC, Greenbelt Maryland
    obs.SIT_Greenwich   # Royal Observatory Greenwich
    obs.SIT_Lowell      # Lowell Observatory, Flagstaff Arizona
    obs.SIT_NAOJ        # National Observatory of Japan, in Mitaka
    obs.SIT_NRL         # Naval Research Laboratory
    obs.SIT_OPAR        # Paris Observatory
    obs.SIT_USNO        # United States Naval Observatory
    obs.SIT_Socorro     # NRAO Socorro




Sources
~~~~~~~

.. code:: python

    ### Fringe finders
    obs.SRC_3C84        # 3C84
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
    
    ### Other
    obs.SRC_Crab        # Crab Nebula / M1
    obs.SRC_SagAstar    # Sagittarius A*, the Galactic Center

