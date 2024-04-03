import os

from platformdirs import user_cache_dir

from lephare import data_manager as dm


def test_data_manager_default_directories():
    """Test DataManager default configuration"""
    new_dm = dm.DataManager()
    assert new_dm.LEPHAREDIR is None
    assert new_dm.LEPHAREWORK is None


def test_data_manager_configure_directories_with_default(unset_env_vars):
    """Test DataManager configure_directories method with default directories"""
    new_dm = dm.DataManager()
    new_dm.configure_directories()
    expected_data_directory = user_cache_dir("lephare", ensure_exists=True) + "/data"
    assert os.getenv("LEPHAREDIR") == expected_data_directory
    assert expected_data_directory == new_dm.LEPHAREDIR
