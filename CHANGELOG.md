# 1.0

2023-09-24

* Minor fixes and tweaks (1.0.2)
    - Merged pull request to remove numpy.warnings filter, due to it being removed from recent np versions and crashing obsplanning on import
    - package multicolorfits has been moved to 'extras' instead of a requirement.  Functions for plotting images (finder plots) will fail without it, but it requires QT - maintaining QT compatibility on various systems and software versions has become a headache, and the rest of the obsplanning package need not be held hostage by the image plotting dependencies.
    

2021-10-03

* Initial release
    - pyephem-based time conversions, observers (telescope), and target source (sky) objects
    - Plotting/visualization tools for target elevation tracks and finder plots
    - Transit time and target visibility calculations for single telescopes or multiple stations (VLBI/VLBA)
    - Miscellaneous coordinate handling, radio astronomy info, and other helpful functions
    - Preliminary versions use functional approach


