Creating Sightlines in yt
=========================

Creating Sightline Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once we have all the necessary metadata for each simulation snapshot saved to the ``galaxy_catalogs`` folder, we're ready to run some sightlines. This is done in two steps. First, we create the endpoints for the sightlines as follows:

.. code-block:: python

 quasarscan.create_qso_endpoints(fullname, redshift, ions)

Where ``fullname`` is the simulation's name (which has to be the same as its name in the metadata catalog), ``redshift`` is the simulation's current redshift, and ``ions`` as a list of ions to track, e.g. ``["H I","O VI","Ne VIII"]`` 

By default, this creates 384 sightlines (384 = 32 * 12, see below), distributed in impact parameter ``r`` from 0 to 2 Rvir. The distribution is designed to be area-sampling, so further sightlines are more likely, though a small bias is added so that there is a nonzero chance of an ``r=0`` line. They are oriented according to discrete ``theta`` and ``phi`` spherical coordinates, corresponding to the polar and azimuthal angles relative to the halo's angular momentum (so ``theta = 0`` is face-on from the direction of the angular momentum, and ``theta = pi/2`` is edge-on). Each sightline randomly selects a starting point by sampling ``theta`` and ``phi``, and then identifying the coordinates of that point on a sphere of size ``6*Rvir``. Then, a midpoint of the sightline is chosen by creating a circle with the radius a random impact parameter as specified above, and choosing a random point along this circle. The line from the startpoint (on the outer sphere) to the midpoint (in the plane normal to the startpoint which goes through the galaxy center) is extended by a factor of 2 to reach an endpoint. This endpoint is close to but not exactly on the outer sphere (because r<<6Rvir).

The sightlines are immediately saved to a file called ``output/<fullname>coldensinfo.txt``. This is a text file, where each line represents a sightline and all its corresponding data. The first several columns are sightline position info, including the literal start and end in cartesian coordinates, the r, theta, and phi values, and the sightline number. The other columns represent the actual data: each ion's column density, ion fraction, and the fraction of that ion which falls into a variety of ``gasbins``, meaning "bins of different gas properties (temperature, density, PI or CI, etc.)". The final three columns are the sightline's overall density, metallicity, and (mass-weighted) temperature. However for now, all the data is filled with simply "-1" values, because we have not yet run the sightlines through the simulation, just determined their position.

Running Sightlines Through the Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most computation-intensive step of ``quasarscan`` is the command:

.. code-block:: python

 quasarscan.run_sightlines(outputfilename,\
						   save_after_num,\
						   parallel,\
						   simulation_dest)

where ``outputfilename`` is the file created with ``create_qso_endpoints``, ``output/<fullname>coldensinfo.txt``, and ``save_after_num`` and ``parallel`` are runtime parameters which will depend on the computer you run on. ``simulation_dest`` is the filepath of the simulation snapshot. If ``parallel`` is ``True``, then it will attempt to run in parallel on all processors, and will run ``save_after_num`` sightlines at once. It is your responsibility, if you have some large number of processors, to give a ``save_after_num`` value which divides your number of processors, otherwise you will not use all of them at all times. Some example scripts for running these jobs on NERSC are provided :doc:`Here (not implemented) <../samplescripts>`, though they will need to be changed to account for your own computers, files, and permissions. 

The filename ``output/<fullname>coldensinfo.txt`` will be overwritten in this process, and the empty "-1" values will be replaced with the corresponding ``float`` values from the sightlines. The different ionization species will be filled in according to the defaults in ``trident``. It's recommended that code groups test out ``trident`` functionality first, (looking at, for example, column density projection plots) before trying to run ``quasarscan``.

Current functionality of this does not create images of spectra, though that feature is in production. Current runtime of this is ~ 2-3 min/sightline. This means this finishes in a reasonable amount of time in parallel, but can be tedious to run without parallelization.

This is the second of two places where access to the simulation snapshots is needed, all further analysis can be done with the .txt files this creates, on home computers without many terabytes of simulation files.

Now that we have our sightlines, let's try to analyze them!
