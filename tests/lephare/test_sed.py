import os
import tempfile

import lephare as lp
import numpy as np
import pytest
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


def test_emplace_back():
    sed = SED("toto", 10, "GAL")
    sed.emplace_back(100, 1)
    assert sed.size() == 1
    assert sed.lamb_flux[0].lamb == 100
    assert sed.lamb_flux[0].val == 1
    assert sed.lamb_flux[0].ori == 1


def test_set_vector():
    sed = SED("toto", 10, "GAL")
    x = np.linspace(100, 500, 1000)
    with pytest.raises(RuntimeError):
        # not the same size
        sed.set_vector(x, np.linspace(0, 1, 500))
    sed.set_vector(x, np.ones_like(x))
    assert len(x) == sed.size()
    for el in sed.lamb_flux:
        assert el.val == 1
        assert el.ori == 1
    # test clearing
    sed.set_vector(x, x)
    assert len(x) == sed.size()


def test_integration():
    sed = SED("toto", 10, "GAL")
    # Tophat filter
    lmin = 200
    lmax = 300
    hat = flt(lmin, lmax, 100)
    # flat SED ov val=1
    x = np.linspace(100, 500, 1000)
    sed.set_vector(x, np.ones_like(x))
    result = sed.integrateSED(hat)
    c = 2.99792458e18
    true_res = [
        lmax - lmin,
        c * (1 / lmin - 1 / lmax),
        0.5 * (lmax**2 - lmin**2),
        lmax - lmin,
        0.5 * (lmax**2 - lmin**2),
        1 / lmin - 1 / lmax,
    ]
    print(result)
    print(true_res)
    assert np.allclose(np.array(result), np.array(true_res), 1.1e-2)
