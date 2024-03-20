import os

import pytest

# Set "LEPHAREDIR" locally for tests if it is not already set.
if "LEPHAREDIR" not in os.environ:
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "..")

# Set "LEPHAREWORK" locally for tests if it is not already set.
if "LEPHAREWORK" not in os.environ:
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "tmp")


@pytest.fixture
def test_data_dir():
    test_dir = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(test_dir, "data")
