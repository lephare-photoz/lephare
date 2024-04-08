import os
import tempfile
import time

import pytest
from lephare import data_manager as dm
from platformdirs import user_cache_dir


def test_data_manager_default_directories(unset_env_vars):
    """Test DataManager default configuration"""
    new_dm = dm.DataManager()
    assert new_dm.LEPHAREDIR is None
    assert new_dm.LEPHAREWORK is None


def test_data_manager_with_predefined_directories(unset_env_vars):
    """Test DataManager with predefined directories"""
    os.environ["LEPHAREDIR"] = "/tmp1"
    os.environ["LEPHAREWORK"] = "/tmp2"
    new_dm = dm.DataManager()
    assert new_dm.LEPHAREDIR == "/tmp1"
    assert new_dm.LEPHAREWORK == "/tmp2"


def test_data_manager_configure_directories_with_default(unset_env_vars):
    """Test DataManager configure_directories method with default directories"""
    new_dm = dm.DataManager()
    new_dm.configure_directories()
    expected_data_directory = user_cache_dir("lephare", ensure_exists=True) + "/data"
    expected_work_directory = user_cache_dir("lephare", ensure_exists=True) + "/work"
    assert os.getenv("LEPHAREDIR") == expected_data_directory
    assert os.getenv("LEPHAREWORK") == expected_work_directory
    assert expected_data_directory == new_dm.LEPHAREDIR
    assert expected_work_directory == new_dm.LEPHAREWORK


def test_create_new_run():
    """Make sure that new run directory structure is created"""
    with pytest.raises(RuntimeError):
        with tempfile.TemporaryDirectory() as tmpdir:
            os.environ["LEPHAREWORK"] = tmpdir
            new_dm = dm.DataManager()
            new_dm.create_new_run()


def test_create_new_run_with_existing_run_directory(unset_env_vars):
    """Make sure that no exceptions are raised if the run directory already exists"""
    new_dm = dm.DataManager()
    new_dm.configure_directories()
    time.sleep(1)
    new_dm.create_new_run()


def test_create_work_subdirectories():
    """Make sure that the expected subdirectories are created"""
    with tempfile.TemporaryDirectory() as tmpdir:
        new_dm = dm.DataManager()
        new_dm.create_work_subdirectories(tmpdir)
        assert os.path.isdir(os.path.join(tmpdir, "filt"))
        assert os.path.isdir(os.path.join(tmpdir, "lib_bin"))
        assert os.path.isdir(os.path.join(tmpdir, "lib_mag"))
        assert os.path.isdir(os.path.join(tmpdir, "zphota"))


def test_create_work_subdirectories_with_existing_directories():
    """Makes sure that no exceptions are raised if the directories already exist
    and ensure that existing files are not overwritten or deleted."""

    # Create a temp directory
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create the `filt` directory and add a test file to it
        os.makedirs(os.path.join(tmpdir, "filt"))
        with open(os.path.join(tmpdir, "filt", "test.txt"), "w") as file:
            file.write("test")

        # Create the other directories
        new_dm = dm.DataManager()
        new_dm.create_work_subdirectories(tmpdir)
        assert os.path.isdir(os.path.join(tmpdir, "filt"))
        assert os.path.isdir(os.path.join(tmpdir, "lib_bin"))
        assert os.path.isdir(os.path.join(tmpdir, "lib_mag"))
        assert os.path.isdir(os.path.join(tmpdir, "zphota"))

        # Ensure the test file is still there
        with open(os.path.join(tmpdir, "filt", "test.txt"), "r") as file:
            assert file.read() == "test"
