# ruff: noqa: E402
# ruff: noqa: F403

import os

global LEPHAREDIR

try:
    LEPHAREDIR = os.environ["LEPHAREDIR"]
except KeyError:
    raise RuntimeError("Environment variable LEPHAREDIR has not been set")  # noqa: B904

try:
    os.mkdir(os.environ["LEPHAREWORK"])
    os.mkdir(os.path.join(os.environ["LEPHAREWORK"], "filt"))
    os.mkdir(os.path.join(os.environ["LEPHAREWORK"], "lib_bin"))
    os.mkdir(os.path.join(os.environ["LEPHAREWORK"], "lib_mag"))
    os.mkdir(os.path.join(os.environ["LEPHAREWORK"], "zphota"))
except FileExistsError:
    pass
except KeyError:
    raise RuntimeError("Environment variable LEPHAREWORK has not been set")  # noqa: B904


from ._lephare import *

# from lephare._lephare import  get_lephare_env
# make LEPHAREDIR and LEPHAREWORK avaliable to the C++ codes
get_lephare_env()  # noqa: F405

from ._flt import *
from ._pdf import *
from ._photoz import *
from ._spec import *
from .filter import *
from .filterSvc import *
from .mag_gal import *
from .magSvc import *
from .sedtolib import *
from .zphota import *
