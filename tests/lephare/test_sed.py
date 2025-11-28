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
    x, y = sed.data()

    def interp(z):
        return np.interp(z, x, y, 0, 0)

    tmp = 4 * np.pi * 100 * (3.086e18) ** 2 / 2.99792458e18
    sed.compute_luminosities()
    # for lmin,lmax in [(2100, 2500), (5500, 6500), (21000, 23000), ]
    res = sciint.quad(interp, 2100, 2500, limit=500, epsabs=1.0e-4, epsrel=1.0e-4)
    assert np.log10(res[0] * 2300**2 / 400 * tmp) == pytest.approx(sed.luv, res[1])

    res = sciint.quad(interp, 5500, 6500, limit=500, epsabs=1.0e-4, epsrel=1.0e-4)
    assert np.log10(res[0] * 6000**2 / 1000 * tmp) == pytest.approx(sed.lopt, res[1])

    res = sciint.quad(interp, 21000, 23000, limit=500, epsabs=1.0e-4, epsrel=1.0e-4)
    assert np.log10(res[0] * 22000**2 / 2000 * tmp) == pytest.approx(sed.lnir, res[1])

    res1 = sciint.quad(interp, 3750, 3950, limit=500, epsabs=1.0e-3, epsrel=1.0e-3)
    res2 = sciint.quad(interp, 4050, 4250, limit=500, epsabs=1.0e-4, epsrel=1.0e-4)
    assert res2[0] / res1[0] == pytest.approx(sed.d4000, res1[1])

    assert sed.ltir == lp.INVALID_VAL

    # need to turn qi as a vector to extract qi[2]
    # edge = 12398.42/13.60
    # res = sciint.quad(interp, sed.lamb_flux[0].lamb, edge, limit=500, epsabs=1.0e-4, epsrel=1.0e-4)
    # print(sed.qi)
    # assert res[0] == pytest.approx(sed.qi[2], res[1])

    sed.read(os.path.join(TESTDATADIR, "sed/dale_1.sed"))
    x, y = sed.data()
    sed.compute_luminosities()
    res = sciint.quad(interp, 80000, 1.0e6, limit=500, epsabs=1.0e-4, epsrel=1.0e-4)
    assert np.log10(res[0] * tmp * 2.99792458e18 / 3.826e33) == pytest.approx(sed.ltir, res[1])

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


def test_resample():
    x1 = np.linspace(0, 10, 11)
    y1 = np.ones_like(x1)
    x2 = x1 + 0.5
    y2 = np.zeros_like(x2)
    z1 = []
    z2 = []
    for i in range(10):
        z1.append(lp.oneElLambda(x1[i], y1[i], 0))
        z2.append(lp.oneElLambda(x2[i], y2[i], 1))
    z = lp.concatenate_and_sort(z1, z2)
    print("z:")
    print([e.lamb for e in z])
    print([e.val for e in z])
    print([e.ori for e in z])
    # resample z1 at the position of z2
    res = lp.SED.resample(z, 0, 0, 10)
    print("res 0:")
    print([e.lamb for e in res])
    print([e.val for e in res])
    print([e.ori for e in res])
    res2 = SED.resample(z, 1, 0, 10)
    print("res 1:")
    print([e.lamb for e in res2])
    print([e.val for e in res2])
    print([e.ori for e in res2])
    for e in res[:-1]:
        assert e.ori == 0
        assert e.val == 1
    # resample z2 at the position of z1
    for e in res2[1:]:
        assert e.ori == 1
        assert e.val == 0


def test_resample2():
    v = np.array(
        [lp.oneElLambda(1, 1, 1), lp.oneElLambda(2, 0, 0), lp.oneElLambda(3, 0, 0), lp.oneElLambda(4, 1, 1)]
    )

    res = lp.SED.resample(v, 1, 1, 5)
    print([e.lamb for e in res])
    print([e.val for e in res])
    print([e.ori for e in res])
    for e in res:
        assert e.val == 1
        assert e.ori == 1

    res = lp.SED.resample(v, 0, 1, 5)
    print([e.lamb for e in res])
    print([e.val for e in res])
    print([e.ori for e in res])

    for e in res[1:-1]:
        assert e.val == 0
        assert e.ori == 0
    for e in (res[0], res[-1]):
        assert e.ori == -99


def test_calc_ph():
    sed = GalSED("", 0)
    hc = 12398.42
    x = np.linspace(0, 1500, 10000)
    # with such a sed, calc_ph just computes the integral of lambda
    # from 0 to the 4 edges
    y = np.ones_like(x) * (hc * 1.6022e-12)
    sed.set_vector(x, y)
    assert sed.qi == [0.0, 0.0, 0.0, 0.0]
    sed.calc_ph()
    for i, e in enumerate([54.42, 24.59, 13.60, hc / 1108.7]):
        print(sed.qi[i])
        assert sed.qi[i] == pytest.approx(0.5 * (hc / e) ** 2, 1.0e-3)
