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
The following example demonstrates how to use Lephare to perform a simple end-to-end flow.

.. 
    TODO - include Minimal_photoz_run notebook

.. 
    Also - do we include the little snippet?
    import lephare as lp
    lp.prepare()
    lp.process()
    matplotlib - plot something <>

Lephare can be used either via a Jupyter notebook or from the command line. However,
the use of the command line executables are generally for legacy purposes.

.. TODO - include a link to test_suite.sh?


Developer guide
===============
Before installing any dependencies or writing code, it's a great idea to create 
a virtual environment. LINCC-Frameworks engineers primarily use conda to manage 
virtual environments. If you have conda installed locally, you can run the following 
to create and activate a new environment.

.. code-block:: bash

    >> conda create env -n <env_name> python=3.10
    >> conda activate <env_name>
    >> conda install cxx-compilers

.. 
    TODO - tabs, and the following line for OSX
    brew install llvm libomp

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
    
    1. The single quotes around `'[dev]'` may not be required for your operating system.
    2. `pre-commit install` will initialize pre-commit for this local repository, 
    so that a set of tests will be run prior to completing a local commit. For more 
    information, see the Python Project Template documentation on `pre-commit 
    <https://lincc-ppt.readthedocs.io/en/latest/practices/precommit.html>`_.
    3. Installing `pandoc` allows you to verify that automatic rendering of Jupyter 
    notebooks into documentation for ReadTheDocs works as expected. For more information, 
    see the Python Project Template documentation on `Sphinx and Python Notebooks 
    <https://lincc-ppt.readthedocs.io/en/latest/practices/sphinx.html#python-notebooks>`_.
