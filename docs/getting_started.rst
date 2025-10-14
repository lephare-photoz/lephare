Getting Started
---------------



Installation
============
LePHARE is distributed with the Python Package Index `(PyPI) <https://pypi.org/project/lephare/>`_, and 
thus the simplest way to install it is with pip:

.. code-block:: bash

    pip install lephare

    
We also reccommend using a conda 
environment to control Python version and isolate your installation:

.. code-block:: bash
    
    # We recommend using Python 3.12
    conda create -n <env_name> python=3.12
    conda activate <env_name>
    pip install lephare
    # We have prepared a number of introductory notebooks. In order to run them
    # you must install jupyterlab with the following commands.
    conda install -c conda-forge jupyterlab
    # And create a kernel which has access to this environment
    python -m ipykernel install --user --name <kernel_name>

A this stage, the following Python snippet should work:

.. code-block:: python

    import lephare as lp

    
.. warning:: 
    For previous users of who have set LEPHAREDIR and LEPHAREWORK
    in their environment these should be unset to avoid clashing versions. 

    Be aware that the structure of the repository has changed with respect to 
    all previous versions.

In order to run LePHARE we need to download auxiliary data such as filters, SEDs, 
and attenuation curves which are not shipped with the code. These are explained
in more detail below but for simplicity you can download everything via Python (~1.3Gb):

.. code-block:: python

    import lephare as lp
    lp.data_retrieval.get_auxiliary_data(clone=True)


The following Python snippet is the most basic example to test the installation has worked. 
This will generate intermediate files and outputs stored in the cache or to 
the user defined storage locations discussed later.
You can also get an example notebook running this code `here <https://github.com/lephare-photoz/lephare/blob/main/docs/notebooks/Minimal_photoz_run.ipynb>`_.


.. code-block:: python

    import lephare as lp
    from astropy.table import Table
    # The following config is highly dependent on your input data and science goals
    # You can change it for your own needs
    config=lp.default_cosmos_config.copy()
    lp.prepare(config)
    # The following example table is in the lephare input format.
    input_table=Table.read(f"{lp.LEPHAREDIR}/examples/COSMOS.in",format="ascii")
    # In the next command output is an astropy.table.Table object with the results
    output,_ = lp.process(config, input_table)
    # One can then inspect, for instance, the first 5 lines of output
    output[:5]
    

This will take over ten minutes to run. To check that everything was successful, 
this workflow should produce a 1 to 1 relationship between the spectroscopic 
redshift `output['ZSPEC']` and predicted redshift `output['Z_BEST']`. 



The Configuration Keywords
==========================

Most users will want to run lephare on their own data set so will have to change the config,
set the filters and modify the input tables.

Taking advantage of the full capabilities of LePHARE will depend on a detailed
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
        'Z_STEP': '0.01,0.,7.', # A very fine redshift grid
    })
