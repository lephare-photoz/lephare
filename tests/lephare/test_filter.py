import os
import warnings

import matplotlib
import numpy as np
import pytest
import requests
from lephare import (
    Filter,
    flt,  # noqa: E402
)
from lephare.filterSvc import FilterSvc, SVOUnavailableError  # noqa: E402

matplotlib.use("Agg")


@pytest.fixture
def filter_file(test_data_dir):
    return os.path.join(test_data_dir, "filt/subaru/IB527.pb")


def test_flt_class(filter_file):
    tophat = flt(100.0, 200.0, 50)
    assert tophat.name == "Heavy"
    print(tophat.lmin(), tophat.lmax())
    assert tophat.width() == 100.0
    #
    f = flt(-1, filter_file, trans=1, calib=1)
    assert f.width() == pytest.approx(241.9479, 1.0e-4)
    assert f.lambdaMean() == pytest.approx(5262.2831, 1.0e-4)
    # testing relative path
    tests_dir = os.path.abspath(os.path.curdir)
    flt_dir = os.path.dirname(filter_file)
    os.chdir(flt_dir)
    print(os.path.abspath(os.path.curdir))
    f2 = flt(-1, "./IB527.pb", trans=1, calib=1)
    assert f2.width != 0.0
    f3 = flt(-1, "IB527.pb", trans=1, calib=1)
    assert f3.width != 0.0
    os.chdir(tests_dir)


def test_filtersvc(test_data_dir, filter_file):
    f = FilterSvc.from_file(filter_file, trans=1, calib=0)
    assert f.width() == pytest.approx(241.9479, 1.0e-4)
    assert f.lambdaMean() == pytest.approx(5262.2831, 1.0e-4)

    try:
        g = FilterSvc.from_svo(-1, "Subaru/Suprime.IB527", "AB")
    except SVOUnavailableError as err:
        # The SVO server is sometimes down: do not fail, but make it visible
        warnings.warn(f"SVO filter service not tested: {err}", UserWarning, stacklevel=1)
    else:
        assert np.allclose(f.data()[1], g.data()[1])
    # filters in COSMOS.para not available to the unit tests
    fltvec = FilterSvc.from_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    assert len(fltvec) == 2


class _FakeResponse:
    def __init__(self, status_code=200, content=b""):
        self.status_code = status_code
        self.content = content
        self.headers = {}


@pytest.mark.parametrize(
    "failure",
    [
        requests.exceptions.ConnectionError("refused"),
        requests.exceptions.ConnectTimeout("timeout"),
        requests.exceptions.ReadTimeout("timeout"),
        requests.exceptions.SSLError("ssl"),
        _FakeResponse(status_code=503, content=b"Service Unavailable"),
        _FakeResponse(status_code=200, content=b"<html>not a votable"),
        # well formed document, but not a VOTable
        _FakeResponse(status_code=503, content=b"<html><body>Service Unavailable</body></html>"),
    ],
)
def test_svo_server_unavailable(monkeypatch, failure):
    """An unreachable or failing SVO server raises SVOUnavailableError."""

    def fake_get(*args, **kwargs):
        if isinstance(failure, Exception):
            raise failure
        return failure

    monkeypatch.setattr(requests, "get", fake_get)
    with pytest.raises(SVOUnavailableError):
        FilterSvc.from_svo(-1, "Subaru/Suprime.IB527", "AB")


@pytest.mark.parametrize("status_code", [200, 400])
def test_svo_invalid_filter_id_raises(monkeypatch, status_code):
    """A VOTable with QUERY_STATUS other than OK is an invalid filter id, whatever the HTTP status."""
    votable = (
        b'<?xml version="1.0"?><VOTABLE><RESOURCE><INFO name="QUERY_STATUS" value="ERROR">'
        b"unknown filter</INFO></RESOURCE></VOTABLE>"
    )
    monkeypatch.setattr(requests, "get", lambda *a, **k: _FakeResponse(status_code, votable))
    with pytest.raises(AssertionError, match="QUERY_STATUS"):
        FilterSvc.from_svo(-1, "Not/AFilter", "AB")


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
    input_args = {"VERBOSE": "YES", "TRANS_TYPE": 42}
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
    assert tophat.name == "Heavy"
