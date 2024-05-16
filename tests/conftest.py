import os

import pytest

# Set "LEPHAREDIR" locally for tests if it is not already set.
if "LEPHAREDIR" not in os.environ:
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "data")

# Set "LEPHAREWORK" locally for tests if it is not already set.
if "LEPHAREWORK" not in os.environ:
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "tmp")


@pytest.fixture
def test_data_dir():
    test_dir = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(test_dir, "data")


@pytest.fixture
def data_registry_file(test_data_dir):
    return os.path.join(test_data_dir, "test_data_registry.txt")


@pytest.fixture
def unset_env_vars():
    os.environ.pop("LEPHAREDIR", None)
    os.environ.pop("LEPHAREWORK", None)


@pytest.fixture
def set_env_vars():
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "tmp")
