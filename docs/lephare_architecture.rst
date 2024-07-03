Advanced Usage via Python Interface
--------------------------------

The best way to learn about the Python interface is through running the example :doc:`notebooks <notebooks>`.

Whether or not you are running via Python or the command line there are four key stages of a LePHARE run.

Core stages of lephare
======================

filter
++++++
This stage reads a configurable set of filters and store a representation of
them for later use.

sedtolib
++++++++
This stage reads a configurable set of SED, compute extinction corrections, and
stores the results into a binary library for later use.

mag_gal
+++++++
This stage uses the preceding outputs to compute expected magnitudes, applying
different rules for galaxies, stars, or QSO objects.

zphota
++++++
This final stage performs chisquare minimization in order to derive photometric
redshifts and other physical parameters.

<placeholder for diagram>

A complete run of lephare is typically composed of the preceeding 4 steps. However,
since these steps are often combined, we provide two convenience functions that
wrap the typical process.

prepare
+++++++

This function wraps the ``filter``, ``sedtolib``, and ``mag_gal`` stages.

process
+++++++

This function wraps the ``zphota`` stage.

Environmental variables
=======================

Lephare uses two environment variables internally - ``LEPHAREDIR`` and ``LEPHAREWORK``. 
``LEPHAREDIR`` defines a local data repository directory that contains the SEDs
and other input datasets that are not modified by lephare in any way.

``LEPHAREWORK`` defines a directory that contains the intermediate files produced
by lephare during processing. In general, the outputs of lephare should be
written to locations specified by the user. The user can specify the output location
either via the config file or as an input parameter to the python methods.

We set default values for the environment variables using the local cache.
*Generally using the defaults is preferred.* However these defaults can be overwritten
by the user to any location. If manually overwritten, care must be taken to ensure
consistency.


Data directory structure
========================

The ``LEPHAREDIR`` data repository contains filters, seds, opacities and other
required datasets. 
When using example config files, this directory structure is created automatically
behind the scenes.

If you wish to use user specified SEDs or filters, you'll need to add those to
the default location in the data repository.
See the auxiliary data repository for examples of where to place these files:
<https://github.com/lephare-photoz/lephare-data>

The ``LEPHAREWORK``  names intermediate files based on values set in the user
specified config file. In order to track run versions, we have created a method
to produce timestamped runs that can be referred back to at a later date.
By default, the `work` directory is symlinked to the most recent run.
See <this notebook> for an example of usage.


Required inputs - filters, seds, and configuration file
=======================================================

We do not expect a typical user to create their own filter or sed files.
Instead, they will be creating an input table and a matching config file.
The config file will specify which files are present in the input data and which
sed lists are being used to fit the data.

By default the lephare input table takes advantage of column ordering to determine
which columns correspond to which band.

Both the filter and sed files will be 2 column ascii tables of wavelength and flux.
For an example of well formatted input files, see the
:doc:`minimal_photoz_run notebook <../notebooks/Minimal_photoz_run>`.

The input configuration file is extensive - for a complete list of the available
options, see the :doc:`configuration file page <../keywords>`. 
As an example for further guidance, see the following in the lephare-data GitHub
repository:
https://github.com/lephare-photoz/lephare-data/blob/main/examples/COSMOS.para


Data retrieval tooling
======================

When using the default sed and filter files, we have built data retrieval
functionally to gather the required data files for a given run.

See the following notebook for a more thorough example of data retrieval with lephare:
:doc:`Data retrieval notebook <../notebooks/Data_retrieval>`

Intermediate files
==================

A large number of intermediate files will be produced by lephare and stored in
the directory defined by the ``LEPHAREWORK`` environmental variable. We generally
don't expect users to interact directly with these binary files. The intermediate
files are required by the C++ code and not easily visualized or modified.

Output files
============

The last stage of lephare, ``zphota``, will produce an output table. 
This table will contain columns specified in a secondary configuration file,
typically named ``output.para``, that lists the requested output columns.
An example can be found here
https://github.com/lephare-photoz/lephare-data/blob/main/examples/output.para

In general this output will contain the best estimate of the redshift alongside
other physical parameters.
