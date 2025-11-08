import os
import tempfile

import lephare as lp
import numpy as np
import pytest
import scipy.integrate as sciint
from lephare import QSOSED, SED, GalSED, StarSED, flt

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


def test_string_to_object():
    for t in ["s", "S", "sOap", "STAR"]:
        a = lp.SED.string_to_object(t)
        assert a == lp.object_type.STAR
    for t in ["g", "G", "GAGA", "GAL"]:
        a = lp.SED.string_to_object(t)
        assert a == lp.object_type.GAL
    for t in ["q", "Q", "QUASI", "QSO"]:
        a = lp.SED.string_to_object(t)
        assert a == lp.object_type.QSO
    with pytest.raises(ValueError):
        _ = lp.SED.string_to_object("wrong")


def sed_get_object_type():
    sed = SED("dummy", 0, "GAL")
    assert sed.get_object_type == lp.object_type.GAL
    sed = SED("dummy", 0, "QSO")
    assert sed.get_object_type == lp.object_type.QSO
    sed = SED("dummy", 0, "STAR")
    assert sed.get_object_type == lp.object_type.STAR


def test_sed_constructors():
    sed = SED("dummy", 10, "GAL")
    assert sed.name == "dummy"
    assert sed.nummod == 10
    assert sed.is_gal()
    sed2 = StarSED(sed)
    assert sed2.is_star()
    assert sed2.name == "dummy"
    assert sed2.nummod == 10
    sed2 = QSOSED(sed)
    assert sed2.is_qso()
    assert sed2.name == "dummy"
    assert sed2.nummod == 10
    sed2 = GalSED(sed)
    assert sed2.is_gal()
    assert sed2.name == "dummy"
    assert sed2.nummod == 10


def test_sed_fromfile():
    sed = SED("toto", 10, "GAL")
    sed.read(os.path.join(TESTDATADIR, "sed/o5v.sed.ext"))
    assert sed.size() == 4860
    d = sed.data()
    y = d[1]
    sed.rescale(2.0)
    y2 = sed.data()[1]
    assert np.allclose(y2, 2 * y)


def test_sed_readwrite():
    star = StarSED("star_test")
    star.read(os.path.join(TESTDATADIR, "sed/WDgd71.sed.ext"))

    with tempfile.TemporaryDirectory() as dir_name:
        file_path = os.path.join(dir_name, "test_star.fits")
        star.writeSED(file_path + ".bin", file_path + ".phys", file_path + ".doc")

        newsed = StarSED("newsed")
        newsed.readSEDBin(file_path + ".bin")
        np.testing.assert_array_equal(newsed.data()[0], star.data()[0])
        np.testing.assert_array_equal(newsed.data()[1], star.data()[1])


def test_integration():
    sed = SED("toto", 10, "GAL")
    # test on a case of exact computation
    x = np.linspace(100, 500, 101)
    sed.set_vector(x, x)
    for lmin, lmax in [(20, 50), (20, 150), (200, 700), (400, 800), (500, 800)]:
        assert sed.integrate(lmin, lmax) == lp.INVALID_VAL

    for lmin, lmax in [
        (200, 300),
        (199, 301),
        (254, 256),
        (100, 500),
    ]:
        result = sed.integrate(lmin, lmax)
        print(result)
        assert result == 0.5 * (lmax**2 - lmin**2)


def test_integration2():
    # test on a more realistic case
    sed = SED("toto", 10, "GAL")
    sed.read(os.path.join(TESTDATADIR, "sed/o5v.sed.ext"))
    x, y = sed.data()

    def interp(z):
        return np.interp(z, x, y, 0, 0)

    for lmin, lmax in [(5000, 7000), (5001, 6999), (4996, 7003)]:
        result = sed.integrate(lmin, lmax)
        #
        hat = flt(lmin, lmax, 100)
        result2 = sed.integrateSED(hat)
        #
        v, ev = sciint.quad(interp, lmin, lmax, limit=500, epsabs=1.0e-4, epsrel=1.0e-4)
        print(result, result2[3])
        print(v, ev)
        assert result == pytest.approx(v, ev)

def test_sedproperties():
    # test on a more realistic case
    sed = GalSED("toto", 10)
    sed.read(os.path.join(TESTDATADIR, "sed/o5v.sed.ext"))
    [luv, lopt, lnir, d4000, ltir] = sed.compute_luminosities()
    sed.SEDproperties()
    assert luv == pytest.approx(sed.luv, 1.e-2)
    assert lopt == pytest.approx(sed.lopt, 1.e-2)
    assert lnir == pytest.approx(sed.lnir, 1.e-2)
    assert d4000 == pytest.approx(sed.d4000, 1.e-2)
    assert ltir ==  lp.INVALID_VAL and sed.ltir == lp.INVALID_VAL
    sed.read(os.path.join("../lephare-data/sed/GAL/DALE/dale_1.sed"))
    [luv, lopt, lnir, d4000, ltir] = sed.compute_luminosities()
    sed.SEDproperties()
    assert luv == pytest.approx(sed.luv, 1.e-2)
    assert lopt == pytest.approx(sed.lopt, 1.e-2)
    assert lnir == pytest.approx(sed.lnir, 1.e-2)
    assert d4000 == pytest.approx(sed.d4000, 1.e-2)
    assert ltir == pytest.approx(sed.ltir, 1.e-2)
