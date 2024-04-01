# ruff: noqa: E402
# ruff: noqa: F403

import datetime
import os

# Why is this global?
global LEPHAREDIR
# Get the default cache location
if os.name == "nt":
    this_os_cache = os.getenv("LOCALAPPDATA")
else:
    this_os_cache = os.getenv("XDG_CACHE_HOME", os.path.expanduser("~/Library/Caches"))

# Check if user defined environment variables present
LEPHAREDIR = os.environ.get("LEPHAREDIR", None)
# If not set set to default and make directories if not present
if LEPHAREDIR is None:
    # Default location in system cache
    LEPHAREDIR = this_os_cache + "/lephare/data"
    # Set environment variable to default
    os.environ["LEPHAREDIR"] = LEPHAREDIR
    if os.path.isdir(LEPHAREDIR):
        pass
    else:
        print(
            f"""Lephare default cache directory being created at
            {LEPHAREDIR}. More than 1Gb may be written there."""
        )
        os.makedirs(LEPHAREDIR)
else:
    print(
        f"""User defined LEPHAREDIR is set. Code runs depend on all required
        auxilliary data being present at {LEPHAREDIR}."""
    )

LEPHAREWORK = os.environ.get("LEPHAREWORK", None)
# If not set thenset to default and make directories if not present
if LEPHAREWORK is None:
    # Default location in system cache
    LEPHAREWORK = this_os_cache + "/lephare/work"
    # Set environment variable to default
    os.environ["LEPHAREWORK"] = LEPHAREWORK
    if os.path.isdir(LEPHAREWORK):
        print(
            f"""Lephare default working directory already exists at
            {LEPHAREWORK}."""
        )
    else:
        # Make a timestamped 'run' directory which is symlinked from default
        now = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
        RUNDIR = f"{this_os_cache}/lephare/runs/{now}"
        print(
            f"""Lephare default working directory being created at
            {LEPHAREWORK}."""
        )
        os.makedirs(RUNDIR)
        os.symlink(RUNDIR, LEPHAREWORK)
        os.makedirs(os.path.join(os.environ["LEPHAREWORK"], "filt"))
        os.makedirs(os.path.join(os.environ["LEPHAREWORK"], "lib_bin"))
        os.makedirs(os.path.join(os.environ["LEPHAREWORK"], "lib_mag"))
        os.makedirs(os.path.join(os.environ["LEPHAREWORK"], "zphota"))
else:
    print(
        f"""User defined LEPHAREWORK is set. All intermediate files will
         be written to {LEPHAREWORK}."""
    )


from ._lephare import *

# from lephare._lephare import  get_lephare_env
# make LEPHAREDIR and LEPHAREWORK avaliable to the C++ codes
get_lephare_env()  # noqa: F405

from ._flt import *
from ._pdf import *
from ._photoz import *
from ._spec import *
from ._version import *
from .filter import *
from .filterSvc import *
from .mag_gal import *
from .magSvc import *
from .sedtolib import *
from .zphota import *
