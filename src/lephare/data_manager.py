import os
import datetime
from platformdirs import user_cache_dir
import warnings


class DataManager:
    def __init__(self, lephare_dir=None, lephare_work_dir=None):
        self.lephare_dir = lephare_dir
        self.lephare_work_dir = lephare_work_dir

    @property
    def LEPHAREDIR(self):
        return self.lephare_dir

    @property
    def LEPHAREWORK(self):
        return self.lephare_work_dir

    def configure_directories(self):
        this_os_cache = user_cache_dir("lephare", ensure_exists=True)

        # Check if user defined environment variables present
        self.lephare_dir = os.getenv("LEPHAREDIR", None)
        # If not set set to default and make directories if not present
        if self.lephare_dir is None:
            # Default location in system cache
            self.lephare_dir = this_os_cache + "/data"
            # Set environment variable to default
            os.environ["LEPHAREDIR"] = self.lephare_dir
            if os.path.isdir(self.lephare_dir):
                pass
            else:
                print(
                    f"""Lephare default cache directory being created at
                    {self.lephare_dir}. More than 1Gb may be written there."""
                )
                os.makedirs(self.lephare_dir)
        else:
            print(
                f"""User defined LEPHAREDIR is set. Code runs depend on all required
                auxiliary data being present at {self.lephare_dir}."""
            )

        self.lephare_work_dir = os.getenv("LEPHAREWORK", None)
        # If not set then set to default and make directories if not present
        if self.lephare_work_dir is None:
            # Default location in system cache
            self.lephare_work_dir = this_os_cache + "/work"
            # Set environment variable to default
            os.environ["LEPHAREWORK"] = self.lephare_work_dir
            if os.path.isdir(self.lephare_work_dir):
                print(
                    f"""Lephare default working directory already exists at
                    {self.lephare_work_dir}."""
                )
            else:
                # Make a timestamped 'run' directory which is symlinked from default
                now = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
                run_directory = f"{this_os_cache}/runs/{now}"
                print(
                    f"""Lephare default working directory being created at
                    {self.lephare_work_dir}."""
                )
                os.makedirs(run_directory)
                os.symlink(run_directory, self.lephare_work_dir)

                work_sub_directories = ["filt", "lib_bin", "lib_mag", "zphota"]
                for sub_dir in work_sub_directories:
                    os.makedirs(os.path.join(self.lephare_work_dir, sub_dir))
        else:
            print(
                f"""User defined LEPHAREWORK is set. All intermediate files will
                be written to {self.lephare_work_dir}."""
            )
