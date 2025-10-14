Advanced Usage via Command Line Interface
=========================================

The four principle command line tools are immediately available after installing 
with pip. They can also be called by building the c++ executables directly but 
that is generally not advised. The four basic operations are summarised below.

Running with the command line
-----------------------------

LePHARE consists of a set of four principle executables:

- ``filter``: read a configurable set of filters and store a representation of them for later use.

- ``sedtolib``: read a configurable set of SED, and store the results into a binary library in a common format for later use.

- ``mag_gal``: use the preceding outputs to compute expected magnitudes and compute extinction corrections, applying different rules for galaxies, stars, or QSO objects.

- ``zphota``: performs chisquare minimization in order to derive photometric redshifts and other physical parameters.

These executables are configurable via a set of keywords, that can be passed at the command line or through a configuration file.
The list of these keywords can be found :doc:`here <keywords>`.

`filter` can often be called with just the config file

.. code-block:: bash

  filter -c config_file.para

`sedtolib` will typically require the config and the object type. Updated library names which 
overwrite the config can also be passed to it.

.. code-block:: bash

  sedtolib -c config_file.para \
    -t G \ # type G for Galaxies
    --GAL_SED COSMOS_MOD.list \ # optional flag to overwrite
    --GAL_LIB LIB_VISTA 

`mag_gal` likewise can be run with just the config and type but is also highly customisable.

.. code-block:: bash

  mag_gal  -c config_file.para \
    -t G \
    --GAL_LIB_IN LIB_VISTA \
    --GAL_LIB_OUT VISTA_COSMOS_FREE \
    --MOD_EXTINC 18,26,26,33,26,33,26,33  \
    --EXTINC_LAW SMC_prevot.dat,SB_calzetti.dat,SB_calzetti_bump1.dat,SB_calzetti_bump2.dat  \
    --EM_LINES EMP_UV  \
    --EM_DISPERSION 0.5,0.75,1.,1.5,2.

Finally `zphota` can be run similarly.

.. code-block:: bash

  zphota -c config_file.para \
  --CAT_IN COSMOS.in \
  --CAT_OUT zphot_vista_adapt.out \
  --ZPHOTLIB VISTA_COSMOS_FREE,ALLSTAR_COSMOS,QSO_COSMOS  \
  --ADD_EMLINES 0,100 \
  --AUTO_ADAPT YES  

Most of the optional flags that are sent to the command line arguments correspond to 
keywords in the config file and effectively override them.

A more detailed example shell script which will run the COSMOS example can be found 
`here <https://github.com/lephare-photoz/lephare/blob/main/docs/test_suite/test_suite.sh>`_.

This `example <https://github.com/lephare-photoz/lephare-data/blob/main/examples/README_full>`_ 
also shows some more advanced features which can be accessed via the command line.

Build only the C++ executables (historical method)
--------------------------------------------------

In addition to the Command Line Interface that is available immediately following 
pip installation some users may want to use the original binaries to run lephare. 
this might be useful if you don't want to depend on the Python interface or
as a means to solve installation issues with an unsuported operating system.
To do so requires making the C++ binaries using the MakeFile. It is possible to
install these binaries using either cmake or make. This is not the recommended 
method but it remains a possibility.

The C++ part of LePhare has no external dependency beyond OpenMP and standard 
compiling and linking libraries. The python part depends on cmake and several 
packages from the python scientific ecosystem (see below).
In order to run make you will need an appropriate compiler for instance gcc on 
a Mac which can be installed with ``brew install gcc``.

Follow the steps below to build with make which may require updating the 
MakeFile for your system:

.. code-block:: bash

   git clone https://gitlab.lam.fr/Galaxies/LEPHARE
   cd LEPHARE
   cd src/lib
   make
   #Add the current directory to the path
   export PATH=$PATH:$(pwd)

Currently, with this method, the code just builds the executables in the *src/lib* directory.
To run code you would also need to download the additional data repository and 
set the environment variables to the correct location

.. code-block:: bash

   git clone https://github.com/lephare-photoz/LEPHARE-data.git
   # LEPHAREDIR must point to the additional data repository
   export LEPHAREDIR=$(pwd)/LEPHARE-data
   # The LEPHAREWORK directory can be anywhere useful
   mkdir lepharework
   export LEPHAREWORK=$(pwd)/LEPHARE-work