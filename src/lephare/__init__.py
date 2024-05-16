# ruff: noqa: E402
# ruff: noqa: F403

# Why is this global?
global LEPHAREDIR

from .data_manager import DataManager

dm = DataManager()
dm.configure_directories()  # noqa: F405
LEPHAREDIR = dm.LEPHAREDIR


from ._lephare import *

# from lephare._lephare import  get_lephare_env
# make LEPHAREDIR and LEPHAREWORK avaliable to the C++ codes
get_lephare_env()  # noqa: F405

from ._flt import *
from ._pdf import *
from ._photoz import *
from ._spec import *
from .data_retrieval import *
from .filter import *
from .filterSvc import *
from .mag_gal import *
from .magSvc import *
from .prepare import *
from .process import *
from .runner import *
from .sedtolib import *
from .zphota import *
