import numpy as np
from lephare import onesource


def test_onesource_creation():
    src = onesource()
    assert src.spec == "1"
    assert src.zs == -99.9
    assert src.cont == 0
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


def test_readsource():
    # Instantiate a source
    src = onesource(101, [0, 0.1, 1])
    # read the source, change Id, attribute flux/err, ...
    src.readsource(
        "10", [-6.414e-32, 1.3182e-31, 1.6905e-31], [1.1022e-31, 9.8579e-32, 5.8665e-32], 6, 2.1, "add"
    )
    assert src.spec == "10"
    assert src.pdfmap[9].size() == 11
    assert np.isclose(src.zs, 2.1)
    assert src.cont == 6
    assert np.isclose(src.ab[0] * 1.0e32, -6.414)
    assert np.isclose(src.sab[0] * 1.0e31, 1.1022)


def test_fltused():
    src = onesource(101, [0, 0.1, 1])
    # upper limit in band 3 using negative error
    src.readsource(
        "10", [-6.414e-32, 1.3182e-31, 1.6905e-31], [1.1022e-31, 9.8579e-32, -5.8665e-32], 7, 2.1, "add"
    )
    # Test without global or forbiden context
    src.fltUsed(-1, -1, 3)
    assert src.cont == 7
    assert np.array_equal(src.busnorma, [1, 1, 0])
    assert np.array_equal(src.busul, [0, 0, 1])
    assert src.nbused == 3
    assert src.nbul == 1
    # Test with a forbidden context removing the first band
    src.fltUsed(-1, 1, 3)
    assert src.cont == 7
    assert np.array_equal(src.busnorma, [0, 1, 0])
    assert np.array_equal(src.busul, [0, 0, 1])
    assert src.nbused == 2
    assert src.nbul == 1
    # Test with a global context using only the second band
    src.fltUsed(2, -1, 3)
    assert src.cont == 2
    assert np.array_equal(src.busnorma, [0, 1, 0])
    assert np.array_equal(src.busul, [0, 0, 0])
    assert src.nbused == 1
    assert src.nbul == 0


def test_convert_mag():
    src = onesource(101, [0, 0.1, 1])
    # upper limit in band 3 using negative error
    src.readsource(
        "10", [-6.414e-32, 1.3182e-31, 1.6905e-31], [1.1022e-31, 9.8579e-32, -5.8665e-32], 7, 2.1, "add"
    )
    src.convertMag()
    assert np.testing.assert_almost_equal(src.mab, [1000, 28.6000467, 28.3299621]) is None
    assert np.testing.assert_almost_equal(src.msab, [1000, 0.81214379, -1]) is None


def test_adapt_mag():
    src = onesource(101, [0, 0.1, 1])
    # upper limit in band 3 using negative error
    src.readsource(
        "10", [6.414e-32, 1.3182e-31, 1.6905e-31], [1.1022e-31, 9.8579e-32, 5.8665e-32], 7, 2.1, "add"
    )
    # keep original flux before offset
    src.keepOri()
    src.convertMag()
    assert np.testing.assert_almost_equal(src.ab_ori, [6.414e-32, 1.3182e-31, 1.6905e-31]) is None
    src.adapt_mag([0, 0, 0])
    assert np.testing.assert_almost_equal(src.ab, [6.414e-32, 1.3182e-31, 1.6905e-31]) is None
    # Offset of +2.5 mag multiplies the flux by 10 (like dividing model by 10)
    src.adapt_mag([2.5, 2.5, 2.5])
    assert np.testing.assert_almost_equal(src.ab, [6.414e-31, 1.3182e-30, 1.6905e-30]) is None
