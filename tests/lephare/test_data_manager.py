import os
import tempfile
import time
from unittest.mock import patch

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


def test_data_manager_with_predefined_dir_directory(unset_env_vars):
    """Test DataManager with predefined dir directory"""
    os.environ["LEPHAREDIR"] = "/tmp1"
    new_dm = dm.DataManager()
    assert new_dm.LEPHAREDIR == "/tmp1"
    assert new_dm.LEPHAREWORK is None


def test_data_manager_with_predefined_work_directory(unset_env_vars):
    """Test DataManager with predefined work directory"""
    os.environ["LEPHAREWORK"] = "/tmp2"
    new_dm = dm.DataManager()
    assert new_dm.LEPHAREDIR is None
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


def test_create_new_run_with_unlinked_work_dir():
    """Make sure that a runtime error is raised when we attempt to create a new
    run directory structure within a pre-existing and un-linked work directory"""
    with pytest.raises(RuntimeError):
        with tempfile.TemporaryDirectory() as tmpdir:
            os.environ["LEPHAREWORK"] = tmpdir
            new_dm = dm.DataManager()
            new_dm.create_new_run()


def test_create_new_run_in_mock_cache_dir(unset_env_vars):
    """Make sure that the expected directory structure is created in a mock cache
    directory"""
    with tempfile.TemporaryDirectory() as tmpdir:
        with patch("lephare.data_manager.user_cache_dir", return_value=tmpdir):
            new_dm = dm.DataManager()
            new_dm.configure_directories()
            assert len(os.listdir(os.path.join(tmpdir, "runs"))) == 1


def test_create_new_run_creates_symlink(unset_env_vars):
    """Make sure that the expected symlinked directory structure is created"""
    with tempfile.TemporaryDirectory() as tmpdir:
        with patch("lephare.data_manager.user_cache_dir", return_value=tmpdir):
            # Set up data manager and make a run
            new_dm = dm.DataManager()
            new_dm.configure_directories()
            new_dm.create_new_run()
            # Check work directory exists
            assert os.path.isdir(os.path.join(tmpdir, "work"))
            # Check the symlink exists
            assert os.path.islink(os.path.join(tmpdir, "work"))
            assert os.readlink(os.path.join(tmpdir, "work")) == os.path.join(
                tmpdir, "runs", os.listdir(os.path.join(tmpdir, "runs"))[0]
            )


def test_create_new_run_overwrites_symlink(unset_env_vars):
    """Make sure that the expected symlinked directory structure is created, and
    the previous symlink is overwritten"""
    with tempfile.TemporaryDirectory() as tmpdir:
        with patch("lephare.data_manager.user_cache_dir", return_value=tmpdir):
            # Set up data manager and make a run
            new_dm = dm.DataManager()
            new_dm.configure_directories()
            new_dm.create_new_run()
            # Sleep for a second to ensure the timestamp is different
            time.sleep(1)
            # Create a new run
            new_dm.create_new_run()
            # Check the symlink has been moved from the original directory to new directory
            sorted_timestamped_dirs = sorted(os.listdir(os.path.join(tmpdir, "runs")))
            assert os.readlink(os.path.join(tmpdir, "work")) != sorted_timestamped_dirs[0]
            assert os.readlink(os.path.join(tmpdir, "work")) != sorted_timestamped_dirs[-1]


def test_create_work_subdirectories():
    """Make sure that the expected subdirectories are created"""
    with tempfile.TemporaryDirectory() as tmpdir:
        new_dm = dm.DataManager()
        new_dm.create_work_subdirectories(tmpdir)
        assert os.path.isdir(os.path.join(tmpdir, "filt"))
        assert os.path.isdir(os.path.join(tmpdir, "lib_bin"))
        assert os.path.isdir(os.path.join(tmpdir, "lib_mag"))
        assert os.path.isdir(os.path.join(tmpdir, "zphota"))
        assert len(os.listdir(tmpdir)) == 4


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


def test_remove_empty_run_directories(unset_env_vars):
    """Make sure that the method to remove empty run directories does remove
    empty directories."""
    with tempfile.TemporaryDirectory() as tmpdir:
        with patch("lephare.data_manager.user_cache_dir", return_value=tmpdir):
            runs_dir = os.path.join(tmpdir, "runs")
            # Set up data manager and make a run
            new_dm = dm.DataManager()
            new_dm.configure_directories()
            new_dm.create_new_run()
            # Assert there's just one run directory
            assert len(os.listdir(runs_dir)) == 1
            # Make more (sleeping to ensure timestamps are different)
            time.sleep(1)
            new_dm.create_new_run()
            time.sleep(1)
            new_dm.create_new_run()
            # Assert there are now three run directories
            assert len(os.listdir(runs_dir)) == 3
            # Remove the empty run directories
            new_dm.remove_empty_run_directories()
            # Assert there's just one run directory (the new one created to preserve symlink)
            assert len(os.listdir(runs_dir)) == 1


def test_remove_empty_skips_nonempty_symlinked_run(unset_env_vars):
    """Make sure that the method to remove empty run directories does not remove
    non-empty directories. This is the case where the non-empty directory is the
    symlinked work directory, so we do not need to generate a new run directory
    to preserve the symlink."""
    with tempfile.TemporaryDirectory() as tmpdir:
        with patch("lephare.data_manager.user_cache_dir", return_value=tmpdir):
            # Set up data manager and make three runs (sleeping for diff timestamps)
            new_dm = dm.DataManager()
            new_dm.configure_directories()
            new_dm.create_new_run()
            time.sleep(1)
            new_dm.create_new_run()
            time.sleep(1)
            new_dm.create_new_run()
            # Assert there are now three run directories
            assert len(os.listdir(os.path.join(tmpdir, "runs"))) == 3
            # Create a file in the most recent run directory
            most_recent_run = sorted(os.listdir(os.path.join(tmpdir, "runs")))[-1]
            with open(os.path.join(tmpdir, "runs", most_recent_run, "filt", "test.txt"), "w") as file:
                file.write("not empty!")
            # Remove the empty run directories
            new_dm.remove_empty_run_directories()
            # Assert there's just one run directory
            assert len(os.listdir(os.path.join(tmpdir, "runs"))) == 1


def test_remove_empty_skips_nonempty_non_symlinked_run(unset_env_vars):
    """Make sure that the method to remove empty run directories does not remove
    non-empty directories. This is the case where the non-empty directory is not
    the symlinked work directory, so we need to generate a new run directory to
    preserve the symlink."""
    with tempfile.TemporaryDirectory() as tmpdir:
        with patch("lephare.data_manager.user_cache_dir", return_value=tmpdir):
            # Set up data manager and make three runs (sleeping for diff timestamps)
            new_dm = dm.DataManager()
            new_dm.configure_directories()
            new_dm.create_new_run()
            time.sleep(1)
            new_dm.create_new_run()
            time.sleep(1)
            new_dm.create_new_run()
            # Assert there are now three run directories
            assert len(os.listdir(os.path.join(tmpdir, "runs"))) == 3
            # Create a file in the first run directory
            first_run = sorted(os.listdir(os.path.join(tmpdir, "runs")))[0]
            with open(os.path.join(tmpdir, "runs", first_run, "filt", "test.txt"), "w") as file:
                file.write("not empty!")
            # Remove the empty run directories
            new_dm.remove_empty_run_directories()
            # Assert there's two run directories: the non-empty one and the new one
            assert len(os.listdir(os.path.join(tmpdir, "runs"))) == 2
