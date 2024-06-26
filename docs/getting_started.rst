Getting Started
---------------

Installation
============
Lephare can be installed with pip:

.. code-block:: bash
    
    pip install lephare

.. note::
    For existing users accustomed to using the command line arguments, those should 
    be immediately available after pip installation.


Example Usage
*************
The following snippet is the most basic example of running lephare end-to-end. You can also get
an example notebook running this code `here <https://github.com/lephare-photoz/lephare/blob/main/docs/notebooks/Minimal_photoz_run.ipynb>`_.


.. code-block:: python

    import lephare as lp
    config=lp.all_types_to_keymap(lp.default_cosmos_config)
    lp.data_retrieval.get_auxiliary_data(keymap=config, additional_files=['examples/COSMOS.in'])
    lp.prepare(config)
    input=Table.read(f"{lp.LEPHAREDIR}/examples/COSMOS.in")
    output, pdfs, zgrid = lp.process(config, input)
    


To ensure that the everything was successful, this workflow should produce a 
a more or less 1 to 1 relationship between the spectroscopic redshift "output['ZSPEC']" 
and predicted redshift "output['Z_BEST']".



.. note

    Lephare can be used either via a Jupyter notebook or from the command line. 
    However, the use of the command line executables are generally for legacy purposes.
    An example of using the command line arguments can be found `here <https://github.com/lephare-photoz/lephare/blob/main/docs/historical_examples/test_suite.sh>`_.

External data and the environment variables
===========================================
LePHARE depends on external data sets such as spectral energy distributions,
filter transmission curves, and attenuation curves. In order to keep the pip
installation light these are now stored in a distinct repository called
`lephare-data <https://github.com/lephare-photoz/lephare-data>`_.



We have built some automatic machinery for downloading the required files 
for a given config `para` file. However, some users may prefer to simply clone
the entire directory. 

.. note::
    lephare uses environment variables to locate the external data and work files.
    If you want to download the full external data repository you must set the environment
    variable `LEPHAREDIR` to its location. If not lephare will use the default cache
    location.

Developer guide
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
or some other preferred location.
One way to acheive this is to use the automatic downloading functionality to 
put all the external data in that location after git cloning the main code:


.. code-block:: python
    
    # Set the LEPHAREWORK dir prior to importing lephare to use your preferred location.
    import os
    os.environ["LEPHAREDIR"]="/path/to/the/preferred/directory/"
    from lephare import data_retrieval as dr
    dr.get_auxiliary_data(clone=False)

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
