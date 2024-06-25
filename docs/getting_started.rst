Getting Started
---------------

Installation
============
Lephare can be installed with pip:

.. code-block:: bash
    
    >> pip install lephare

.. note::
    For existing users accustomed to using the command line arguments, those should 
    be immediately available after pip installation.


Example Usage
*************
This is the most basic example of running lephare end-to-end.


.. codeblock

    import lephare as lp
    config=lp.all_types_to_keymap(lp.default_cosmos_config)
    lp.data_retrieval.get_auxiliary_data(keymap=config, additional_files=['examples/COSMOS.in'])
    lp.prepare(config)
    input=Table.read(f"{lp.LEPHAREDIR}/examples/COSMOS.in")
    output, pdfs, zgrid = lp.process(config, input)
    


To ensure that the installation was successful, this workflow should produce a reasonable 
a more or less 1 to 1 relationship between the spectroscopic redshift "output['ZSPEC']" 
and predicted redshift "output['Z_BEST']".

For further example usage, see the Minimal_photoz_run notebook.

.. note

    Lephare can be used either via a Jupyter notebook or from the command line. 
    However, the use of the command line executables are generally for legacy purposes.


Developer guide
===============
Before installing any dependencies or writing code, it's a great idea to create 
a virtual environment. LINCC-Frameworks engineers primarily use conda to manage 
virtual environments. If you have conda installed locally, you can run the following 
to create and activate a new environment.

.. tabs::

    .. tab:: bash

        .. code-block:: bash

            >> conda create env -n <env_name> python=3.10
            >> conda activate <env_name>
            >> conda install cxx-compilers

    .. tab:: OSX

        .. code-block:: bash

            >> conda create env -n <env_name> python=3.10
            >> conda activate <env_name>
            >> conda install cxx-compilers
            >> brew install llvm libomp


Once you have created a new environment, you can install this project for local 
development using the following commands:

.. code-block:: bash

    >> pre-commit install
    >> conda install pandoc

If you wish to incorporate your changes to the main branch, please make a fork of 
the repository and then create a pull request. 

If you are having problems with installations, there is a list of known issues here. 
If you can’t find a solution, feel free to `create an issue in the lephare repository 
<https://github.com/lephare-photoz/lephare/issues>`_.

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
