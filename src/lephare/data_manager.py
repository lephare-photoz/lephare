import os
import datetime
from platformdirs import user_cache_dir


class DataManager:
    def __init__(self):
        self.LEPHAREDIR = None
        self.LEPHAREWORK = None

    def configure_directories(self):
        this_os_cache = user_cache_dir("lephare", ensure_exists=True)

        # Check if user defined environment variables present
        self.LEPHAREDIR = os.getenv("LEPHAREDIR", None)
        # If not set set to default and make directories if not present
        if self.LEPHAREDIR is None:
            # Default location in system cache
            self.LEPHAREDIR = this_os_cache + "/data"
            # Set environment variable to default
            os.environ["LEPHAREDIR"] = self.LEPHAREDIR
            if os.path.isdir(self.LEPHAREDIR):
                pass
            else:
                print(
                    f"""Lephare default cache directory being created at
                    {self.LEPHAREDIR}. More than 1Gb may be written there."""
                )
                os.makedirs(self.LEPHAREDIR)
        else:
            print(
                f"""User defined LEPHAREDIR is set. Code runs depend on all required
                auxiliary data being present at {self.LEPHAREDIR}."""
            )

        self.LEPHAREWORK = os.getenv("LEPHAREWORK", None)
        # If not set then set to default and make directories if not present
        if self.LEPHAREWORK is None:
            # Default location in system cache
            self.LEPHAREWORK = this_os_cache + "/work"
            # Set environment variable to default
            os.environ["LEPHAREWORK"] = self.LEPHAREWORK
            if os.path.isdir(self.LEPHAREWORK):
                print(
                    f"""Lephare default working directory already exists at
                    {self.LEPHAREWORK}."""
                )
            else:
                # Make a timestamped 'run' directory which is symlinked from default
                now = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
                RUNDIR = f"{this_os_cache}/lephare/runs/{now}"
                print(
                    f"""Lephare default working directory being created at
                    {self.LEPHAREWORK}."""
                )
                os.makedirs(RUNDIR)
                os.symlink(RUNDIR, self.LEPHAREWORK)
                os.makedirs(os.path.join(os.environ["LEPHAREWORK"], "filt"))
                os.makedirs(os.path.join(os.environ["LEPHAREWORK"], "lib_bin"))
                os.makedirs(os.path.join(os.environ["LEPHAREWORK"], "lib_mag"))
                os.makedirs(os.path.join(os.environ["LEPHAREWORK"], "zphota"))
        else:
            print(
                f"""User defined LEPHAREWORK is set. All intermediate files will
                be written to {self.LEPHAREWORK}."""
            )
