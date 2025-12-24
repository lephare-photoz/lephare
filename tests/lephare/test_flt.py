import os

import lephare as lp
import pytest

TESTDIR = os.path.abspath(os.path.dirname(__file__))


def test_read_filters_from_file():
    with pytest.raises(ValueError):
        lp.read_filters_from_file("non_existing_file")


def test_flt_constructor():
    with pytest.raises(ValueError):
        f = lp.flt(0, "non_existing_file", 0, 0)

    filename = os.path.join(TESTDIR, "..", "data", "filt", "subaru", "IB527.pb")
    f = lp.flt(0, filename, 0, 0)
    assert len(f.lamb_trans) > 0
