import datetime
import os

from platformdirs import user_cache_dir


class DataManager:
    def __init__(self):
        self.lephare_dir = None
        self.lephare_work_dir = None

    @property
    def LEPHAREDIR(self):
        return self.lephare_dir

    @property
    def LEPHAREWORK(self):
        return self.lephare_work_dir

    def configure_directories(self):
        # Check if user defined environment variables present
        self.lephare_dir = os.getenv("LEPHAREDIR", None)
        # If not, set set to default and make subdirectories
        if self.lephare_dir is None:
            default_os_cache = user_cache_dir("lephare", ensure_exists=True)
            # Default location in system cache
            self.lephare_dir = f"{default_os_cache}/data"
            # Set environment variable to default
            os.environ["LEPHAREDIR"] = self.lephare_dir
            # If this directory does not exist, create it
            if not os.path.isdir(self.lephare_dir):
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
        # TODO Add a check for existing /work/* subdirectories and /runs directory
        if self.lephare_work_dir is None:
            default_os_cache = user_cache_dir("lephare", ensure_exists=True)
            # Default location in system cache
            self.lephare_work_dir = f"{default_os_cache}/work"
            # Set environment variable to default
            os.environ["LEPHAREWORK"] = self.lephare_work_dir
            if os.path.isdir(self.lephare_work_dir):
                print(
                    f"""Lephare default working directory already exists at
                    {self.lephare_work_dir}."""
                )
            else:
                self.create_new_run()
        else:
            print(
                f"""User defined LEPHAREWORK is set. All intermediate files will
                be written to {self.lephare_work_dir}."""
            )

    # TODO Define a method to create a new run directory
    def create_new_run(self):
        default_os_cache = user_cache_dir("lephare", ensure_exists=True)
        now = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
        run_directory = f"{default_os_cache}/runs/{now}"
        print(
            f"""Lephare default working directory being created at
            {self.lephare_work_dir}."""
        )
        os.makedirs(run_directory)
        os.symlink(run_directory, self.lephare_work_dir)

        work_sub_directories = ["filt", "lib_bin", "lib_mag", "zphota"]
        for sub_dir in work_sub_directories:
            #! Should this be creating the subdirectories in the run directory, `run_directory`?
            os.makedirs(os.path.join(self.lephare_work_dir, sub_dir))
