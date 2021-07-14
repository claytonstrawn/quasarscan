Sorting and Filtering Data
==========================

After we have run ``run_sightlines`` and collected all our data, we need to be able to see the results! QUASARSCAN can automatically read all the data stored in your ``quasarscan_data`` folder and present it in a variety of useful ways in 2, 3, 4, or 5 dimensions, either for publishable papers or just data exploration. In addition to plotting the sightlines from the simulation, it is possible to load observational sightline data, or explore metadata as created by ``create_metadata``, to determine which simulations are worth analyzing in more detail.

First, we will discuss filtering your data so that plots can be maximally useful and only show what is needed to analyze a particular point.

Loading Data
^^^^^^^^^^^^

To begin analysis, the first step is to instantiate a ``MultiQuasarSpherePlotter`` object. This will automatically load all of the sightline data in your ``quasarscan_data`` folder. 

.. code-block:: python
    
 mq = quasarscan.create_mq(loadsim = "all",\
                           loadobs = 'all',\
                           loadempty = 'none',\
                           average = 'median')

Where ``loadsim`` and ``loadobs`` refer to which simulations/observations to load (specified by giving any part of the simulation/observation's ``fullname``), ``loadempty`` refers again to simulations, but loads simulations from their existing metadata, rather than because it has existing sightlines (thus it's "empty"). ``average`` refers to the default averaging functon. The options are:

* "mean"
* "median" (if "median" is given, one can also give either an integer percentile which represents the size of errorbars, default is 25)

The averaging function can also be specified during the ``plot`` step.

Constraining Data
^^^^^^^^^^^^^^^^^

While the software by default will load all saved simulations, it is best to constrain and sort by different quantities. For example, we can restrict to low/high redshift, low/high stellar mass, low/high SFR, etc. This is implemented in two ways. First, we can restrict the whole list and just cut out most of the data that doesn't fit our specifications. Second, we can sort the remaining data by another variable and put it into several bins. We will call each galaxy snapshot's "collection of lines" a QuasarSphere.

First, you can always see the current loaded QuasarSpheres and metadata. For example, to see all the star formation rates and stellar masses, run:

.. code-block:: python
    
 mq.list_all_quasar_spheres('SFR','Mstar')

To constrain the data, use the code below:

.. code-block:: python

 mq.constrain_current_quasar_array(constrain_criteria,
                                   bins=None,**kwargs)

where ``constrain_criteria`` refers to the metadata in question. For example, to restrict to redshifts between 0.5 and 1.0, use the code:

.. code-block:: python

 mq = quasarscan.create_mq()
 mq.constrain_current_quasar_array('redshift',[0.5,1.0])


You can also restrict by certain kinds of simulations or simulation number. To do this, run with a list of acceptable string values. For example, to restrict to the VELA simulation, number 1,2, and 3, we can run:

.. code-block:: python

 mq = quasarscan.create_mq()
 mq.constrain_current_quasar_array('simname',['VELA'])
 mq.constrain_current_quasar_array('simnum',['01','02','03'])

After doing this, you can always reset to the full list by running 

.. code-block:: python

 mq.reset_current_quasar_array()

The full list of arguments arguments for ``constrain_current_quasar_array`` is below.

* constrain_criteria: ``string``. Can refer to snapshot param (``Rvir``, ``Mvir``, ``redshift``, ``sfr`` etc.) or stringparam (``simname``, ``simnum``, ``version``, etc.)
* bins=None: ``list`` of either two numbers, which are low and high edges of bin, if ``constrain_criteria`` is a snapshot param, or multiple accepted strings, if ``constrain_criteria`` is a stringparam.
* qtype='all': ``string`` By default, sort all observations and empty QuasarSpheres (see "Advanced Plotting Techniques") alongside simulations. Can change to 'sim', 'obs', or 'empty' to only effect those lists. If ``constrain_criteria`` is a stingparam, this defaults to 'sim'.
* at_end=False: ``boolean`` or ``float`` If False, use the current value of the snapshot param. If ``float`` between 0 and 1, use the value of this simulation at expansion parameter a = ``at_end``. Simulations which do not run to that time are excluded.
* exclude=False: ``boolean``. If True, restrict to all values outside of bin, instead of inside. This is most useful to exclude a single simulation with a stringparam.
* split_even=False: ``boolean`` or ``string``. If False, use value in ``bins``. If ``split_even='high'``, create a bin of all simulations with ``constrain_criteria`` higher than the median and sort using that. If ``split_even='low'``, create a bin of all simulations with ``constrain_criteria`` lower than the median and sort using that.
* set_main_array=False: ``boolean``. If True, restrict the main array with this call, not just the current array. After running this, ``reset_current_quasar_array`` will no longer reset to before this line.

Sorting Data 
^^^^^^^^^^^^

After appropriately restricting your data, you will probably want to keep track of multiple bins of galaxy snapshots at once. This function is run in a very similar way, via a ``constrain_criteria`` such as mass, redshift, or SFR. In this case, one can give multiple bins. It returns a tuple of (0) a string describing the bins (which will be used as their label in the graphs below), (1) a list of bin edges, and (2) the list of arrays of QuasarSpheres. This is conventionally referred to as "lq" for "labels,bins,QuasarSpheres". The below will sort your data into three bins, galaxy snapshots with 

.. code-block:: python

 lq = mq.sort_by('sfr',[0.1,1.0,10.0,np.inf])

Unlike ``constrain_current_quasar_array``, ``sort_by`` does not effect the internal list of quasarspheres, it just distributes the existing list into multiple sublists and returns them. Note that any galaxies which do not fit in any bin, or have a ``nan`` for their ``criteria`` are simply not returned.

One useful keyword argument of ``sort_by`` is ``split_even=n``. This will split the list into ``n`` bins of equal size, without needing to specify the bins in advance. The bin edges will be thus somewhat arbitrary, but each bin all have a meaningful amount of data and will be useful for distinguishing low, medium, and high mass galaxies (for example):

.. code-block:: python

 lq = mq.sort_by('Mstar',split_even = 3)

The full list of arguments arguments for ``sort_by`` is below.

* criteria:``string``. Can refer to snapshot param (``Rvir``, ``Mvir``, ``redshift``, ``sfr`` etc.) or stringparam (``simname``, ``simnum``, ``version``, etc.)
* bins=[0,np.inf]:``list`` edges of bins, if snapshot param, or list of accepted bins, if stringparam. The default arg just checks that the value in question exists but does not filter for any value.
* at_end=False: ``boolean`` or ``float`` If False, use the current value of the snapshot param. If ``float`` between 0 and 1, use the value of this simulation at expansion parameter a = ``at_end``. Simulations which do not run to that time are not returned.
* split_even=False: ``boolean`` or ``int``. If False, use values in ``bins``. If ``int``, sort the simulation data into that many equal-sized bins
* reverse=False: ``boolean``, if ``True`` return the bins in reverse order (by default, they are returned low to high)
* sort_w_qtype='sim': ``string``, only used if ``split_even`` is ``True``. Split by putting ``qtype=('sim', 'obs', or 'empty')`` into equal sized bins.

To use the bins, we will keep this ``lq`` object and bring plug it into a plot function.





