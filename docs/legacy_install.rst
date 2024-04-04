Legacy installation
========================================================================================

Some users may want to use the original binaries to run lephare. 
To do so requires making the c++ binaries using the MakeFile.

Requirements
------------

The C++ part of LePhare has no external dependency beyond OpenMP and standard compiling and linking libraries. The python part depends on cmake and several packages from the python scientific ecosystem (see below).
In order to run make you will need an appropriate compiler for instance gcc on a Mac which can be installed with ``brew install gcc``.

Legacy installation methods
---------------------------

Build only the C++ executables with make (historical method)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Follow the steps below:

.. code-block:: console

   >> git clone https://gitlab.lam.fr/Galaxies/LEPHARE
   >> cd LEPHARE
   >> cd src/lib
   >> make
   >> #Add the current directory to the path
   >> export PATH=$PATH:$(pwd)

Currently, with this method, the code just builds the executables in the *source* directory.

Build the C++ executables and the python module with cmake and setuptools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently, the following sequence allows you to build and install the C++ executables, together with the python package.
The pybind11 header files are provided in the *extern* directory, but after a `git clone` one needs to download these as well, with the command `git submodule update --init --recursive`.

.. code-block:: console

   >> conda create env -n <env_name> python=3.10
   >> conda activate <env_name>
   >> pip install -e .'[dev]'
   >> pre-commit install
   >> conda install pandoc

The C++ executables are installed in the bin sub directory. 
The python module is installed in the conventional `site-packages` area of the system. 
Note that you may need to create the LEPHAREWORK/filt, LEPHAREWORK/lib_bin and LEPHAREWORK/lib_mag by hand.


A user manual is provided through **LePhare_documentation.pdf** in the *doc* directory
If Doxygen is installed (see https://www.doxygen.nl/manual/install.html), code documentation is also available: execute ```python setup.py doc``` and an html sub-directory will be available for browsing in the *doc* directory.

Legacy code overview and main features
--------------------------------------

LePhare consists of a set of executables:

- ``filter``: read a configurable set of filters and store a representation of them for later use.

- ``sedtolib``: read a configurable set of SED, compute extinction corrections, and store the results into a binary library for later use.

- ``mag_gal``: use the preceding outputs to compute expected magnitudes, applying different rules for galaxies, stars, or QSO objects.

- ``zphota``: performs chisquare minimization in order to derive photometric redshifts and other physical parameters.

These executables are configurable via a set of parameters, that can be passed at the command line or through a configuration file.
The list of these parameters can be found  here <parameters>.

Quick Start
-----------

In order to run these commands you must first get the example data sets `available here <https://github.com/lephare-photoz/LEPHARE-data>`_. Next, a parameter file named ``COSMOS.para`` is shipped with LePhare, together with the input files needed for a quick analysis of the COSMOS field, under the examples directory.

.. code-block:: console

   >> git clone https://github.com/lephare-photoz/LEPHARE-data.git
   >> cd LEPHARE-data/examples
   >> #We must set the environment variables to find the required files
   >> export LEPHAREDIR=<path to current dir>
   >> export LEPHAREWORK='<path to where execution will save data>'

First, set OMP multithreaded to the value whished for, e.g. :

.. code-block:: console

   >> export OMP_NUM_THREADS='10'

The typical sequence of execution is then:
- read the filters and compile them into one file stored in ``$LEPHAREWORK/filt``

.. code-block:: console

   >> filter -c COSMOS.para

- read the galaxy templates (as used in Ilbert et al. 2013) and store them in ``$LEPHAREWORK/lib_bin``

.. code-block:: console

   >> sedtolib -c COSMOS.para -t G -GAL_SED COSMOS_MOD.list  -GAL_LIB LIB_VISTA

- use the galaxy templates and filters to derive a library of predicted magnitudes and store it in ``$WORK/lib_mag`` (the parameters correspond to enabling emission lines correlated to UV light + free factor in scaling these lines, more information in the original ``LePhare_documentation.pdf``)

.. code-block:: console

   >> mag_gal  -c COSMOS.para -t G -GAL_LIB_IN LIB_VISTA -GAL_LIB_OUT VISTA_COSMOS_FREE -MOD_EXTINC 18,26,26,33,26,33,26,33  -EXTINC_LAW SMC_prevot.dat,SB_calzetti.dat,SB_calzetti_bump1.dat,SB_calzetti_bump2.dat  -EM_LINES EMP_UV  -EM_DISPERSION 0.5,0.75,1.,1.5,2. -Z_STEP 0.04,0,6

- proceed in the same way for stellar and QSO templates

.. code-block:: console

   >> sedtolib -c COSMOS.para -t S -STAR_SED STAR_MOD_ALL.list
   >> mag_gal -c COSMOS.para -t S -LIB_ASCII YES -STAR_LIB_OUT ALLSTAR_COSMOS 
   >> # AGN models from Salvato 2009
   >> sedtolib -c COSMOS.para -t Q -QSO_SED  $LEPHAREDIR/sed/QSO/SALVATO09/AGN_MOD.list
   >> mag_gal -c COSMOS.para -t Q -MOD_EXTINC 0,1000  -EB_V 0.,0.1,0.2,0.3 -EXTINC_LAW SB_calzetti.dat -LIB_ASCII NO  -Z_STEP 0.04,0,6

- finally proceed to photometric redshift estimation

.. code-block:: console

   >> zphota -c COSMOS.para -CAT_IN COSMOS.in -CAT_OUT zphot_short.out -ZPHOTLIB VISTA_COSMOS_FREE,ALLSTAR_COSMOS,QSO_COSMOS  -ADD_EMLINES 0,100 -AUTO_ADAPT YES   -Z_STEP 0.04,0,6

- a python script is available to perform a quick diagnostics

.. code-block:: console

   >> python figuresLPZ.py zphot_short.out

.. toctree::