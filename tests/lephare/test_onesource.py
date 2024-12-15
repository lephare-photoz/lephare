import numpy as np
from lephare import onesource


def test_onesource_creation():
    src = onesource()
    assert src.spec == "1"
    assert src.zs == -99.9
    assert src.cont == 0
    assert src.closest_red == 0
    assert np.array_equal(src.zmin, [-99.9, -99.9, -99.9])
    assert np.array_equal(src.chimin, [1.0e9, 1.0e9, 1.0e9])
    assert np.array_equal(src.indmin, [-99, -99, -99])
    assert np.array_equal(src.imasmin, [-99, -99, -99])

    src = onesource(31, [0, 0.5, 1])
    assert src.pos == 31
    assert src.pdfmap[9].size() == 3
    assert src.spec == "1"
    assert src.zs == -99.9
    assert src.cont == 0
    assert src.closest_red == 0
    assert np.array_equal(src.zmin, [-99.9, -99.9, -99.9])
    assert np.array_equal(src.chimin, [1.0e9, 1.0e9, 1.0e9])
    assert np.array_equal(src.indmin, [-99, -99, -99])
    assert np.array_equal(src.imasmin, [-99, -99, -99])


def test_verbosity():
    src = onesource()
    src.set_verbosity(True)
    assert src.get_verbosity()
    src.set_verbosity(False)
    assert not src.get_verbosity()


def test_onesource_set_priors():
    src = onesource()
    src.setPriors([0.0, 1000.0], [0, 1000])
    assert np.array_equal(src.priorLib, [0.0, 0.0, 1000.0, 1000.0])


def test_validlib():
    s = onesource()
    s.closest_red = 0.15
    assert s.validLib([0, 0.05, 0.15, 0.2], True) == [
        2,
    ]
    assert s.validLib([0, 0.05, 0.15, 0.2], False) == [0, 1, 2, 3]
