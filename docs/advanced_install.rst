Advanced Installation with full clones of repositories
======================================================
LePHARE is split across two repositories one for the 
`code <https://github.com/lephare-photoz/lephare>`_ and one for the
auxiliary `data <https://github.com/lephare-photoz/lephare-data>`_. 
Advanced users or developers may want to edit the files in 
both. In this section we explain a developer install and clone of both 
git repositories in order to allow modification of the code and auxiliary data.


Auxiliary Data and the Environment Variables
********************************************

LePHARE depends on auxiliary data sets such as spectral energy distributions,
filter transmission curves, and attenuation curves. In order to keep the pip
installation light these are now stored in a distinct repository called
`lephare-data <https://github.com/lephare-photoz/lephare-data>`_.

LePHARE uses environment variables to locate the external data and work files. These are set by default to your cache. Both can be set if preferred. The two environment variables are the following:

* `LEPHAREDIR` is the location of the auxiliary data.
* `LEPHAREWORK` is the location of the intermediate files produced during a LePHARE run.



.. note::
    If you want to download the external data to a specific location you must set the
    environment variable `LEPHAREDIR` to its location. This must be done prior to 
    importing lephare in a python session. If not LePHARE will use the default cache
    location on your system.
    
Some users may prefer to simply clone the entire  auxiliary data directory:

.. code-block:: bash
    
    git clone https://github.com/lephare-photoz/lephare-data
    # Set the LEPHAREDIR to this data location
    export LEPHAREDIR=$PWD/lephare-data

    
We have also built some automatic machinery for downloading the
required files for a given config `para` file. The automatic download
functionality can also be used to download all external data. In the
following snippet we show how you might set the `LEPHAREDIR` to a new
location and download all the auxiliary data there:

.. code-block:: python

    import os
    os.environ['LEPHAREDIR']='/path/to/my/preferred/data_directory/'
    os.environ['LEPHAREWORK']='/path/to/my/preferred/working_directory/'
    # You must import lephare after setting the variables
    import lephare as lp
    # If you do not set a config input to the following function in gets everything.
    # Data will be put in $LEPHAREDIR
    lp.data_retrieval.get_auxiliary_data(clone=False)
    # Setting clone=True would use a git clone which may be faster but will only run 
    # on an empty directory.


Developer Installation
**********************
The developer install is required for editing the code but can also be useful
on systems that do not have PyPI binaries and for systems that are not well tested.
Before installing any dependencies or writing code, it's a great idea to create 
a virtual environment. LINCC-Frameworks engineers primarily use conda to manage 
virtual environments. If you have conda installed locally, you can run the following 
to create and activate a new environment. We then recommend installing in 
editable mode with the `-e` option so that any changes are immediately propagated.

.. tabs::

    .. tab:: bash

        .. code-block:: bash

            conda create -n <env_name> python=3.12
            conda activate <env_name>
            git clone https://github.com/lephare-photoz/lephare.git
            cd lephare
            git submodule update --init --recursive
            conda install -c conda-forge cxx-compiler
            pip install -e .'[dev]'

    .. tab:: OSX

        .. code-block:: bash

            conda create -n <env_name> python=3.12
            conda activate <env_name>
            brew install llvm libomp
            git clone https://github.com/lephare-photoz/lephare.git
            cd lephare
            git submodule update --init --recursive
            conda install -c conda-forge cxx-compiler
            pip install -e .'[dev]'


At this stage running the tests is a good way to check everything is working:

.. code-block:: bash

    python -m pytest tests

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
