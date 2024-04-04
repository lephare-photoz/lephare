.. lephare documentation main file.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to lephare's documentation!
========================================================================================
**WARNING: This is a work in progress, and should ne be used by anyone at this stage. The current LePHARE code is to be found at https://gitlab.lam.fr/Galaxies/LEPHARE .**

LePHARE (PHotometric Analysis for Redshift Estimation) is a Python package built on a complete rewrite in C++ of the `Fortran code <https://www.cfht.hawaii.edu/~arnouts/LEPHARE/acknowledgement.html>`_ LePhare.
LePHARE computes photometric redshifts and physical parameters by fitting spectral energy distributions (SED) to a dataset of photometric fluxes or apparent magnitudes.

Installation
------------

The simplest way to install lephare for most users is with pip:

.. code-block:: console

   >> pip install lephare

If you prefer to use binary executables from the command line you may wish to conduct a `legacy installation <https://gitlab.lam.fr/Galaxies/LEPHARE/>`_.  

Dev Guide - Getting Started
---------------------------

Before installing any dependencies or writing code, it's a great idea to create a
virtual environment. LINCC-Frameworks engineers primarily use `conda` to manage virtual
environments. If you have conda installed locally, you can run the following to
create and activate a new environment.

.. code-block:: console

   >> conda create env -n <env_name> python=3.10
   >> conda activate <env_name>


Once you have created a new environment, you can install this project for local
development using the following commands:

.. code-block:: console

   >> pip install -e .'[dev]'
   >> pre-commit install
   >> conda install pandoc


Notes:

1) The single quotes around ``'[dev]'`` may not be required for your operating system.
2) ``pre-commit install`` will initialize pre-commit for this local repository, so
   that a set of tests will be run prior to completing a local commit. For more
   information, see the Python Project Template documentation on
   `pre-commit <https://lincc-ppt.readthedocs.io/en/latest/practices/precommit.html>`_.
3) Installing ``pandoc`` allows you to verify that automatic rendering of Jupyter notebooks
   into documentation for ReadTheDocs works as expected. For more information, see
   the Python Project Template documentation on
   `Sphinx and Python Notebooks <https://lincc-ppt.readthedocs.io/en/latest/practices/sphinx.html#python-notebooks>`_.


.. toctree::
   :hidden:

   Home page <self>
   Notebooks <notebooks>
   Legacy installation <legacy_install>
   Keywords <keywords>
   C API Reference <doxygen_output/c_lib/library_root>
   Python API Reference <autoapi/index>

