import matplotlib

matplotlib.use("Agg")

import os

import lephare as lp

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


def test_spec_plotspec():
    lp._spec.plotspec(os.path.join(TESTDATADIR, "example.spec"))
    assert True
