Getting Started
===============

When you want to do a scan of a certain simulation, you'll need to start with a few pieces of information. We need to find the position of the galaxy in question, it's bulk velocity, it's angular momentum axis, and its virial radius, to start. It is also nice to know things like its star formation rate, total stellar mass, etc. because you can organize sightlines later by binning across these variables. We will call these properties of a single halo in a single snapshot it's "metadata".

This metadata is stored in a secondary folder, intended to be on the same level as quasarscan, called quasarscan_data. Inside it is a folder called galaxy_catalogs, and it will hold all metadata tables. These are organized by the name of the simulation (see below for format).

Creating Metadata
^^^^^^^^^^^^^^^^^

If you don't already have all the data or processing systems available to you to create this metadata, we can attempt to do it using yt. To do so, all you need is the simulation's name and its file path. We organize simulation names as follows:

.. code-block:: python

 fullname = <simname>_<version>_<code>_<simnum>

Where <simname> refers to the general project name (we have in the past used VELA and AGORA), <version> refers to the current version, if it has been redone multiple times (this is best used for comparing across different implementations of the same initial conditions, which has been one of the main uses of the code), <code> refers to the kind of code used, recognizable by yt, and simnum is an arbitrary numbering. If your simulation consists of one huge box with many halos, you might refer to each one with a different number in the same snapshot file, while if your simulation is of several smaller boxes with one halo each, each box might be numbered. You tell us how to sort the simulation halos, the only parameter with separate importance is the <code>, so please don't feed the wrong one.

To create the necessary metadata from the simulation snapshot, run:

.. code-block:: python

 quasarscan.create_metadata(fullname, filepath) 

By default, this will attempt to find the largest high-resolution halo in the simulation (meaning, the largest halo but with the constraint that it is only allowed to use halos which reach the simulation's minimum resolution at least once within the halo.) The list of metadata information it will save is:

* a (expansion parameter)
* Rvir
* Mvir
* center (x,y,z coordinates)
* L_x,L_y,L_z (angular momentum direction)
* Lmag (magnitude of angular momentum)
* Mstar (stellar mass within 0.1 Rvir)
* Mgas (gas mass within halo)
* Mdm (Dark Matter mass within halo)
* SFR (current star formation rate)
* bulk_velocity_x,bulk_velocity_y,bulk_velocity_z

This is the first of two places where ``quasarscan`` will require access to your saved simulation snapshots. Most of the analysis can be done without them.

Using Your Own Metadata
^^^^^^^^^^^^^^^^^^^^^^^

You may want to track additional information which we do not naturally track, or you may have already worked on calculating all these (very standard) quantities by your own methods. If so, you are encouraged to use your own metadata tables, as they may be more reliable. To add a metadata table, just add a text file with the name <simname>_<version>_<code>_<simnum>_metadata.txt to a folder called quasarscan_data/galaxy_catalogs/<simname>_<version>_<code>_<simnum>. 

The format of this text file is 

.. code-block:: python

 a, center_x, center_y, center_z, L_x, L_y, L_z, bulk_velocity_x, bulk_velocity_y, bulk_velocity_z, Rvir, Mvir, Lmag, sfr, Mstar, Mgas, compaction_stage
 0.1, 0.00523361, -0.0111065, -10.1503, -0.147827, -0.695122, 0.703529, 12.592, -14.2341, -0.147827, 4.75, 3041110000.0, 1.49325, 0.0201637, 9227770.0, 277344000.0, pre
 ...

In this case, we are tracking one additional parameter called the "compaction_stage" (see Huertas-Company et al 2018). The first line gives the names of the values in each column. The only requirement is that the first column be the expansion parameter a. It is required that we have the center (in "unitary" units, see yt) and the virial radius (in kpc) in order to run.

Once the metadata exists, you're ready to run some sightlines. 

Looking Up Metadata
^^^^^^^^^^^^^^^^^^^

Once your metadata exists, it can be nice to access it in python. You can do so at any time using

.. code-block:: python

 quasarscan.get_value(quantity,fullname,redshift)

Where ``quantity`` is which piece of metadata to access (e.g. ``'SFR'``, ``'Rvir'``, ``'Mstar'``), ``fullname`` is as defined above, and ``redshift`` is the current simulation redshift. It will return the data for the snapshot at the closest redshift to the one entered, as long as it is off by less than 0.2 in expansion factor ``a``. Currently, these are in units of kpc for Rvir, Msun for all masses, and Msun/yr for SFR, however it remains to be done to implement these as unitful quantities using ``unyt``.