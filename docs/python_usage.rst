Advanced Usage via Python Interface
-----------------------------------


The best way to learn about the Python interface is through running the 
example :doc:`notebooks`. We have made an intermediate notebook to demonstrate how a new user can download
the required filters and run lephare on their own data. You can download the notebook from
`here <https://github.com/lephare-photoz/lephare/blob/main/docs/notebooks/README.md>`_ 
which should directly run following installation of lephare via pip.
You can also see the executed notebook with all outputs :doc:`here <notebooks/Typical_use_case.ipynb>`.

.. note::
    Lephare can be used either via the Python interface or from the command line. 
    These both have different use cases. The Python interface is newer and has greater flexibility.

The notebooks are the best way to learn how to use the Python interface. 
The `Python API <https://lephare.readthedocs.io/en/latest/autoapi/index.html>`_ is useful for understanding all the individual classes and methods.
We hope that class and function doctstrings will be informative and intuitive. 
Please let us know via a GitHub issue if any functionality is unclear.

Notebooks
=========

We explain the Python interface mainly through the following example notebooks. 
They range from very simple examples largely to test the installation to 
demonstrating advanced access to the low level functionality. 
They should all be able to run independently following installation.

.. toctree::

    Minimal full run <notebooks/Minimal_photoz_run>
    Typical use case <notebooks/Typical_use_case>
    Typical use case: physical parameters <notebooks/Typical_use_case_physicalParameters>
    Detailed run <notebooks/Detailed_run>
    Developer example: Building a list of onesources <notebooks/Building_list_of_onesources>
    Developer example: Example of usage of magSvc <notebooks/Example_of_usage_of_magSvc>
    Developer example: Data retrieval <notebooks/Data_retrieval>
    Developer example: Testing fit of one object of the catalogue <notebooks/Testing_fit_of_one_object_of_the_catalogue>
    Developer example: Example of the cosmology class <notebooks/Example_cosmo>
    Developer example: Demonstration of SED objects <notebooks/Example_SED_manipulation>
