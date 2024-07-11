.. lephare documentation main file.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

LePHARE (PHotometric Analysis for Redshift Estimation)
========================================================================================

.. warning::

   This is a work in progress, and should not be used by anyone at this stage. The 
   current LePHARE code is to be found at `https://gitlab.lam.fr/Galaxies/LEPHARE 
   <https://gitlab.lam.fr/Galaxies/LEPHARE>`_.


.. image:: https://avatars.githubusercontent.com/u/165841626?s=400&u=ff86bd4c19a9d36958cf1b47d84849dbe25c274a&v=4
   :alt: LePHARE logo
   :align: center
   :width: 150px

LePHARE computes photometric redshifts and physical parameters by fitting spectral energy distributions (SED) 
to a dataset of photometric fluxes or apparent magnitudes. LePHARE is a Python package built on 
a complete rewrite in C++ of the `Fortran code <https://www.cfht.hawaii.edu/~arnouts/LEPHARE>`_.

.. toctree::
   :hidden:

   Home page <self>
   Getting Started <getting_started>
   Advanced usage via Python interface <lephare_architecture>
   Advanced usage via Command Line Interface <legacy_install>
   All notebooks <notebooks>
   Detailed documentation <original>
   Keywords <keywords>
   C API Reference <doxygen_output/c_lib/library_root>
   Python API Reference <autoapi/index>
   Known issues <known_issues>
