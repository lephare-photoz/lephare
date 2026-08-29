import os

import lephare as lp
import pytest

TESTDIR = os.path.abspath(os.path.dirname(__file__))


def test_read_filters_from_file():
    with pytest.raises(ValueError):
        lp.read_filters_from_file("non_existing_file")

    flts = lp.read_filters_from_file(os.path.join(TESTDIR, "..", "data", "filt", "LSST_FILTERS.dat"))
    assert len(flts) == 6
    lp.write_output_filter("tmp.dat", "tmp.doc", flts)
    assert os.path.exists("tmp.dat")
    assert os.path.exists("tmp.doc")
    # need to build an absolute path
    flts2 = lp.read_doc_filters(os.path.join(os.environ["PWD"], "tmp.doc"))
    assert len(flts2) == 6
    os.remove("tmp.dat")
    os.remove("tmp.doc")
    # usual case with only the basename expected in LEPHAREWORK/filt
    flts3 = lp.read_doc_filters("filter_cosmos")
    assert len(flts3) == 2


def test_flt_constructor():
    with pytest.raises(ValueError):
        f = lp.flt(0, "non_existing_file", 0, 0)

    filename = os.path.join(TESTDIR, "..", "data", "filt", "subaru", "IB527.pb")
    f = lp.flt(0, filename, 0, 0)
    assert len(f.lamb_trans) > 0
