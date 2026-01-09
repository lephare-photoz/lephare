# ruff: noqa: E402
# ruff: noqa: F403

# Why is this global?
global LEPHAREDIR

from .data_manager import DataManager

dm = DataManager()
dm.configure_directories()  # noqa: F405
LEPHAREDIR = dm.LEPHAREDIR

from ._lephare import *

# import explicitly the internal variables that we
# need to expose for testing and documentation
from ._lephare import (  # noqa: F401
    _closeAge,
    _emission_lines,
    _empirical_ratio,
    _empirical_ratio2,
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
from .runner import *
from .sedtolib import *
from .zphota import *
