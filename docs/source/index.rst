Welcome to quasarscan's documentation!
======================================

Introduction
^^^^^^^^^^^^

Quasarscan is a python project which is designed to aid in analysis of the Circumgalactic Medium (CGM) of galaxy simulations. While there are a great number of ways to access information from these simulations and create plots and analysis (see yt), comparison with observations can be quite tricky, due to the fact that in the CGM we very rarely see full imaging of emission lines, only singular, scattered absorption lines. This project was created to replicate that experience, by running a random sample of several hundred absorption lines and tracking their properties, including both observable ones (e.g. column densities, metallicities) and non-observable ones (e.g. ion fractions, temperature bins, etc.)

This software creates those sightlines, calculates those values, and can plot the results in a variety of interesting and intuition-building ways. It is designed to work on any simulation code readable by ``yt``, and requires both ``yt`` and ``trident`` to operate. Importantly, (and unlike ``yt`` or ``trident``) it is designed to allow computation to take place on one computer (usually a supercomputer, where modern simulations tend to be stored), save small, easy-to-read output files, and then plotting and analysis can take place on another computer which does not need to have access to the original simulation files.

Installation
^^^^^^^^^^^^

To install quasarscan, visit https://github.com/claytonstrawn/quasarscan.git and clone the repository. 
    git clone https://github.com/claytonstrawn/quasarscan.git



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   gettingstarted
   creatingsightlines
   sortingfiltering
   plottingdata
   advancedplotting
   spectra

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
