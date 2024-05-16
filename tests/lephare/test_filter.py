import os

import matplotlib
import numpy as np
import pytest
from lephare import (
    Filter,
    flt,  # noqa: E402
)
from lephare.filterSvc import FilterSvc  # noqa: E402

matplotlib.use("Agg")


@pytest.fixture
def filter_file(test_data_dir):
    return os.path.join(test_data_dir, "filt/subaru/IB527.pb")


def test_flt_class(filter_file):
    tophat = flt(100.0, 200.0, 50)
    assert tophat.name == ""
    # this is wrong but an effect of improper C++ code
    assert tophat.width() == 101.0
    #
    f = flt(-1, filter_file, trans=1, calib=1)
    assert f.width() == pytest.approx(241.9479, 1.0e-4)
    assert f.lambdaMean() == pytest.approx(5262.2831, 1.0e-4)


def test_filtersvc(test_data_dir, filter_file):
    f = FilterSvc.from_file(filter_file, trans=1, calib=0)
    assert f.width() == pytest.approx(241.9479, 1.0e-4)
    assert f.lambdaMean() == pytest.approx(5262.2831, 1.0e-4)

    g = FilterSvc.from_svo(-1, "Subaru/Suprime.IB527", "AB")
    if g is not None:
        assert np.allclose(f.data()[1], g.data()[1])
    else:
        print("SVO not tested due to exception raised : server not reachable?")
    # filters in COSMOS.para not available to the unit tests
    fltvec = FilterSvc.from_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    assert len(fltvec) == 2


def test_filter_base_class(test_data_dir, set_env_vars):
    """Simple test to ensure that we can create an instance of a Filter object."""
    config_file_path = os.path.join(test_data_dir, "examples/COSMOS.para")
    filter = Filter(config_file=config_file_path)
    filter.run()
    assert len(filter.keymap)


def test_filter_with_kwargs(test_data_dir, set_env_vars):
    """Simple test to ensure that we can create an instance of a Filter object
    and that we can pass kwargs when instantiating the object."""
    config_file_path = os.path.join(test_data_dir, "examples/COSMOS.para")
    input_args = {"verbose": True, "TRANS_TYPE": 42}
    filter = Filter(config_file=config_file_path)
    filter.run(**input_args)
    assert len(filter.keymap)
    assert filter.verbose
    assert filter.keymap["TRANS_TYPE"].value == "42"


def test_flt_plot_filter_curve():
    """Simple test that exercises the plot_filter_curve method."""
    tophat = flt(100.0, 200.0, 50)
    tophat.plot_filter_curve()
    tophat.plot_filter_curve(normed=True)

    # this is assertion doesn't mean anything, just here to assert _something_.
    assert tophat.name == ""
