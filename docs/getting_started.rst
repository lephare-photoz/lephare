Getting Started
---------------



Installation
============
The simplest way to install LePHARE is with pip:

.. code-block:: bash
    
    pip install lephare

.. note::
    For existing users accustomed to using the command line arguments, those should 
    be immediately available after pip installation.


The following Python snippet is the most basic example of running lephare end-to-end. 
You can also get an example notebook running this code `here <https://github.com/lephare-photoz/lephare/blob/main/docs/notebooks/Minimal_photoz_run.ipynb>`_.


.. code-block:: python

    import lephare as lp
    from astropy.table import Table
    # The following config is highly dependent on your input data and science goals
    config=lp.all_types_to_keymap(lp.default_cosmos_config)
    lp.data_retrieval.get_auxiliary_data(keymap=config, additional_files=['examples/COSMOS.in'])
    lp.prepare(config)
    # The following example table is in the lephare input format.
    input_table=Table.read(f"{lp.LEPHAREDIR}/examples/COSMOS.in")
    output, pdfs, zgrid = lp.process(config, input_table)
    

This will take over ten minutes to run. To check that everything was successful, 
this workflow should produce a 1 to 1 relationship between the spectroscopic 
redshift `output['ZSPEC']` and predicted redshift `output['Z_BEST']`. Most users
will want to run lephare on their own data set so will have to change the config,
set the filters and modify the input tables.

Example Usage
*************

We have made an intermediate notebook to demonstrate how a new user can download
the required filters and run lephare on their own data. You can download the notebook
`here <https://github.com/lephare-photoz/lephare/blob/main/docs/notebooks/Intermediate_usage_notebook.ipynb>`_ 
which should directly run following installation of lephare via pip.
You can also see the executed notebook with all outputs :doc:`here <notebooks/Intermediate_usage_notebook>`.

.. note::
    Lephare can be used either via the Python interface or from the command line. 
    The use of the command line executables are generally for legacy purposes.
    An example of using the command line arguments can be found `here <https://github.com/lephare-photoz/lephare/blob/main/docs/historical_examples/test_suite.sh>`_.


Auxiliary Data and the Environment Variables
===========================================
LePHARE depends on auxiliary data sets such as spectral energy distributions,
filter transmission curves, and attenuation curves. In order to keep the pip
installation light these are now stored in a distinct repository called
`lephare-data <https://github.com/lephare-photoz/lephare-data>`_.


We have built some automatic machinery for downloading the required files 
for a given config `para` file. However, some users may prefer to simply clone
the entire directory. The automatic download functionality can also be used to
download all external data.

.. note::
    lephare uses environment variables to locate the external data and work files.
    These are set by default to your cache.
    If you want to download the external data to a specific location you must set the
    environment variable `LEPHAREDIR` to its location. This must be done prior to 
    importing lephare in a python session. If not lephare will use the default cache
    location on your system.

In the following snippet we show how you might set the `LEPHAREDIR` to a new location 
and download all the auxiliary data there:

.. code-block:: python

    import os
    os.environ['LEPHAREDIR']='/path/to/my/preferred/directory/'
    # You must import lephare after setting the variables
    import lephare as lp
    # If you do not set a config input to the following function in gets everything.
    lp.data_retrieval.get_auxiliary_data(clone=False)

* `LEPHAREDIR` is the location of the auxiliary input data.
* `LEPHAREWORK` is the location of the intermediate files produced during a lephare run.

Both can be set if preferred or left to the default location in the user cache.


Advanced Usage
==============

Taking advantage of the advanced capabilities of LePHARE will depend on a detailed
understanding of the configurations which can be specified by text file or via a dictionary 
in Python. In the later stages of the documentation we cover the various options
that can be specified via :doc:`keywords <keywords>`.

For an example text file see the COSMOS example `here <https://github.com/lephare-photoz/lephare-data/blob/main/examples/COSMOS.para>`_.

One way to set config values is to start with the default cosmos config 
dictionary which is shipped with the Python by default and to update those elements 
you want to change. In the following Python snippet we start with the default
COSMOS config and update the redshift grid using the `Z_STEP` keyword to a finer
grid which would increase accuracy but take longer to execute:

.. code-block:: python

    import lephare as lp
    config=lp.default_cosmos_config.copy()
    config.update({
        'Z_STEP': '0.001,0.,7.', # A very fine redshift grid
    })

Developer Guide
===============
Before installing any dependencies or writing code, it's a great idea to create 
a virtual environment. LINCC-Frameworks engineers primarily use conda to manage 
virtual environments. If you have conda installed locally, you can run the following 
to create and activate a new environment. We then recommend installing in 
editable mode with the `-e` option so that any changes are immediately propagated.

.. tabs::

    .. tab:: bash

        .. code-block:: bash

            conda create env -n <env_name> 
            conda activate <env_name>
            conda install cxx-compilers # May not be required for linux
            git clone https://github.com/lephare-photoz/lephare.git
            cd lephare
            git submodule update --init --recursive
            conda install -c conda-forge cxx-compiler
            pip install -e .'[dev]'

    .. tab:: OSX

        .. code-block:: bash

            conda create env -n <env_name> 
            conda activate <env_name>
            conda install cxx-compilers
            brew install llvm libomp
            git clone https://github.com/lephare-photoz/lephare.git
            cd lephare
            git submodule update --init --recursive
            conda install -c conda-forge cxx-compiler
            pip install -e .'[dev]'


At this stage running the tests is a good way to check everything is working:

.. code-block:: bash

    pytest tests

Once you have created a new environment, you can install precommit and pandoc 
which will help you to run precommit checks and create the documentation locally:

.. code-block:: bash

    pre-commit install
    conda install pandoc

Developers can also build the documentation in the following way:

.. code-block:: bash
    
    cd docs/
    pip install -r requirements.txt #install sphinx dependencies
    make html

The doc entry will then be located at `../_readthedocs/html/index.html`. The 
documentation includes a rendering of the notebooks, which thus need to be 
executed. You can bypass this stage by replacing `make html`` above by 
`make no-notebooks`. Executing `make` will list further options.


If you wish to incorporate your changes to the main branch, please make a fork of 
the repository and then create a pull request. 

If you are having problems with installations, there is a list of known issues `here <known_issues.rst>`_. 
If you can’t find a solution, feel free to `create an issue in the lephare repository 
<https://github.com/lephare-photoz/lephare/issues>`_.

Some developers who are familiar with the original version of the code may
want to have all the external data present in the same repository as the code
or some other preferred location. They could set the `LEPHAREDIR` to the code 
location and then use the automatic downloading functionality to put all
the auxiliary data there as it was in the previous versions.


.. note::
    The single quotes around `'[dev]'` may not be required for your operating system.

    `pre-commit install` will initialize pre-commit for this local repository, 
    so that a set of tests will be run prior to completing a local commit. For more 
    information, see the Python Project Template documentation on `pre-commit 
    <https://lincc-ppt.readthedocs.io/en/latest/practices/precommit.html>`_.

    Installing `pandoc` allows you to verify that automatic rendering of Jupyter 
    notebooks into documentation for ReadTheDocs works as expected. For more information, 
    see the Python Project Template documentation on `Sphinx and Python Notebooks 
    <https://lincc-ppt.readthedocs.io/en/latest/practices/sphinx.html#python-notebooks>`_.

    The environment variables `LEPHAREDIR` and `LEPHAREWORK` are set on import
    in Python. Care must be taken not to reset after importing.

    It remains possible to build the C++ code using either make or cmake directly.
    This is not recommended and will likely require OS specific changes. It may be 
    useful on unusual systems where we do not support compilation.
