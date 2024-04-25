import datetime
import os

from platformdirs import user_cache_dir

from lephare._lephare import get_lephare_env


class DataManager:
    def __init__(self):
        self.lephare_dir = os.getenv("LEPHAREDIR", None)
        self.lephare_work_dir = os.getenv("LEPHAREWORK", None)

    @property
    def LEPHAREDIR(self):  # noqa: N802
        return self.lephare_dir

    @property
    def LEPHAREWORK(self):  # noqa: N802
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

        # If LEPHAREWORK is not set then set to default and make directories
        if self.lephare_work_dir is None:
            # get <default_cache> locations
            default_os_cache = user_cache_dir("lephare", ensure_exists=True)

            # Location of work linked dir <default_cache>/work and the timestamped directory
            symlink_work_directory = f"{default_os_cache}/work"

            # first remove the symlink if it already exists
            if os.path.islink(symlink_work_directory):
                # os.unlink(symlink_work_directory) #We no longer make a new link
                print(
                    f"""Default work cache at {symlink_work_directory}
                    is already linked. This is linked to the run directory:
                    {os.readlink(symlink_work_directory)}"""
                )
            else:
                # create a runs directory in the default cache locations
                os.makedirs(f"{default_os_cache}/runs", exist_ok=True)

                # create a timestamped directory under runs
                now = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
                run_directory = f"{default_os_cache}/runs/{now}"
                os.makedirs(run_directory, exist_ok=True)

                # create the subdirectories in the new run directory
                self.create_work_subdirectories(run_directory)
                os.symlink(run_directory, symlink_work_directory)

            # set the LEPHAREWORK environment variable to the <default_cache>/work symlink
            os.environ["LEPHAREWORK"] = symlink_work_directory

            # set the instance variable to the <default_cache>/work symlink
            self.lephare_work_dir = symlink_work_directory
        else:
            # the lephare work dir env var is set, create subdirectories if needed
            self.create_work_subdirectories(self.lephare_work_dir)
            print(
                f"""User defined LEPHAREWORK is set. All intermediate files will
                be written to {self.lephare_work_dir}."""
            )

    def create_new_run(self):
        """Create a timestamped directory to contain the output from the current run.
        The newly created timestamped directory is symlinked to the path defined
        by the LEPHAREWORK environment variable."""

        lephare_work_dir = os.getenv("LEPHAREWORK", None)

        # if the LEPHAREWORK environment variable is not a symlink,
        # then this directory structure is unusual and we cannot create a new run.
        # We'll raise an exception and direct the user to the documentation.
        if not os.path.islink(f"{lephare_work_dir}"):
            # TODO include link to documentation
            raise RuntimeError(
                """The current directory structure does not support the
                               creation of new runs. Please refer to the documentation for
                               information on how to set up the directory structure."""
            )

        # given that LEPHAREWORK is a symlink, create a new timestamped run directory
        now = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
        run_directory = os.path.realpath(f"{lephare_work_dir}/../{now}")
        print(f"Creating new run directory at {run_directory}.")
        os.makedirs(run_directory, exist_ok=True)

        # remove the existing `work` symlink and create a new one
        if os.path.islink(lephare_work_dir):
            os.unlink(lephare_work_dir)
        os.symlink(run_directory, lephare_work_dir)

        # create the subdirectories in the new run directory
        self.create_work_subdirectories(run_directory)

    def create_work_subdirectories(self, parent_dir):
        """Creates the required work subdirectories in the parent directory if they
        are not present. No action is taken if the subdirectories already exist."""
        work_sub_directories = ["filt", "lib_bin", "lib_mag", "zphota"]
        for sub_dir in work_sub_directories:
            os.makedirs(os.path.join(parent_dir, sub_dir), exist_ok=True)


def check_lephare_directories():
    dm = DataManager()
    if dm.LEPHAREDIR is None or dm.LEPHAREWORK is None:
        dm.configure_directories()
    get_lephare_env()
    return dm.LEPHAREDIR, dm.LEPHAREWORK
