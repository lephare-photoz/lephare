"""
LePHARE is a code for computing photometric redshifts and physical parameters

For full documentation, see: https://lephare.readthedocs.io/

Since 1999, LePHARE has been a code for computing photometric redshifts and physical parameters
by fitting spectral energy distributions (SEDs) to a dataset of photometric fluxes or apparent
magnitudes. LePHARE was originally written in Fortran (Arnouts et al. 1999; Ilbert et al. 2006).
It has been completely rewritten in C++ with a Python interface, with a larger team of developers
including J. Cohen-Tanugi and R. Shirley.

In order to run LePHARE we need to download auxiliary data such as filters, SEDs,
and attenuation curves which are not shipped with the code. These are explained
in more detail below but for simplicity you can download everything via Python (~1.3 Gb):

    import lephare as lp
    lp.data_retrieval.get_auxiliary_data(clone=True)

The following Python snippet is the most basic example to test that the installation has worked.
This will generate intermediate files and outputs stored in the cache or in user-defined storage
locations. You can also get an example notebook running this code here.

    import lephare as lp
    from astropy.table import Table

    # The following config is highly dependent on your input data and science goals
    # You can change it for your own needs
    config = lp.default_cosmos_config.copy()
    lp.prepare(config)

    # The following example table is in the lephare input format
    input_table = Table.read(f"{lp.LEPHAREDIR}/examples/COSMOS.in", format="ascii")

    # In the next command output is an astropy.table.Table object with the results
    output, _ = lp.process(config, input_table)

    # One can then inspect, for instance, the first 5 lines of output
    output[:5]

This workflow may take over ten minutes to run. To check that everything was successful,
this example should produce a 1-to-1 relationship between the spectroscopic redshift
output['ZSPEC'] and predicted redshift output['Z_BEST'].
"""
# ruff: noqa: E402
# ruff: noqa: F403

# Why is this global?
global LEPHAREDIR

from .data_manager import DataManager

dm = DataManager()
dm.configure_directories()  # noqa: F405
LEPHAREDIR = dm.LEPHAREDIR

from ._set_omp_num_threads import _set_omp_num_threads

_set_omp_num_threads()

from ._lephare import *

# import explicitly the internal variables and functions
# that we need to expose for testing and documentation
from ._lephare import (  # noqa: F401
    _closeAge,
    _emission_lines,
    _empirical_ratio,
    _empirical_ratio_ori,
    _ga_2q_val,
    _ga_H_val,
    _ga_HeI_val,
    _ga_lamb,
    _ga_total,
    _read_ages_from_file,
)
from ._version import *

# make LEPHAREDIR and LEPHAREWORK avaliable to the C++ codes
get_lephare_env()  # noqa: F405

from ._flt import *
from ._onesource import *
from ._pdf import *
from ._photoz import *
from ._plot_utils import *
from ._spec import *
from .data_retrieval import *
from .default_cosmos_config import *
from .filter import *
from .filter_extinc import *
from .filterSvc import *
from .mag_gal import *
from .magSvc import *
from .prepare import *
from .process import *
from .reddening import *
from .runner import *
from .sedtolib import *
from .zphota import *
