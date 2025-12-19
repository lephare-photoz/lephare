import os
import tempfile

import lephare as lp
import numpy as np
import pytest
import scipy.integrate as sciint
from lephare import QSOSED, SED, GalSED, StarSED, flt

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


def test_apply_ext():
    onext = lp.ext("", 0)
    x1 = np.linspace(0, 50, 20)
    y1 = np.ones_like(x1) / 0.4
    onext.set_vector(x1, y1)

    sed = SED()
    x2 = np.linspace(10, 20, 10)
    y2 = x2

    sed.set_vector(x2, y2)
    sed.apply_extinction(1.0, onext)
    res1 = sed.lamb_flux
    print([el.lamb for el in res1])
    res1y = [el.val for el in res1]
    print(res1y)
    assert np.allclose(res1y, np.array(y2) / 10.0)


def test_generate_lines():
    sed = lp.GalSED("", 0)
    mnuv_int = -3
    nuvr = 3
    scalefac = 1
    sed.generateEmEmpUV(mnuv_int, nuvr)
    if nuvr < 4:
        scalefac = 10 ** (-0.4 * mnuv_int - 6.224) / 2.85
    vals = [ell.val for ell in sed.fac_line]
    assert np.allclose(np.array(vals), scalefac * np.array(lp.empirical_ratio))
    sfr = 10**50.0
    sed.generateEmEmpSFR(sfr, nuvr)
    if nuvr < 4:
        scalefac = 10 ** (np.log10(sfr) + 41.27 - np.log10(4 * np.pi * 100 * 3.08568**2) - 36) / 2.85
    vals = [ell.val for ell in sed.fac_line]
    print(vals)
    print(scalefac * np.array(lp.empirical_ratio2))
    assert np.allclose(np.array(vals), scalefac * np.array(lp.empirical_ratio2))


# def test_apply_ext_lines():
#     onext = lp.ext("", 0)
#     x1 = np.linspace(1000, 17000, 20)
#     y1 = np.ones_like(x1) / 0.4
#     onext.set_vector(x1, y1)

#     sed = SED()
#     x2 = np.linspace(10, 20, 10)
#     y2 = x2
#     sed.set_vector(x2, y2)
#     sed.fac_line = lp._ga_total

#     sed.apply_extinction_to_lines(.1, onext)
#     res1 = sed.fac_line
#     assert np.allclose(res1, np.array(lp._ga_total) / 10.0)


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


def test_sumspectra():
    sed1 = GalSED("toto", 10)
    sed1.read(os.path.join(TESTDATADIR, "sed/o5v.sed.ext"))
    x1, y1 = sed1.data()
    sed2 = GalSED("toto", 11)
    sed1.sumSpectra(sed2, 0)  # case of empty additional sed
    for xx1, yy1, o1 in zip(x1, y1, sed1.lamb_flux):
        assert o1.lamb == xx1
        assert o1.val == yy1
    sed2.read(os.path.join(TESTDATADIR, "sed/o5v.sed.ext"))
    sed1.sumSpectra(sed2, 0)  # case of nul rescale factor
    for xx1, yy1, o1 in zip(x1, y1, sed1.lamb_flux):
        assert o1.lamb == xx1
        assert o1.val == yy1

    # easy case : sum the same SED with a rescal factor of 2
    # here in practice there is no interpolation and so the
    # the input and output have the same size
    sed1 = GalSED("toto", 10)
    sed1.read(os.path.join(TESTDATADIR, "sed/o5v.sed.ext"))
    sed2 = GalSED("toto", 11)
    sed2.read(os.path.join(TESTDATADIR, "sed/o5v.sed.ext"))
    sed1.sumSpectra(sed2, 2)
    x2, y2 = sed2.data()
    x1, y1 = sed1.data()
    print(len(x1), len(x2))
    assert np.allclose(np.array(y1), 3 * np.array(y2))

    # harder case : different SEDs
    sed1 = GalSED("toto", 10)
    sed1.read(os.path.join(TESTDATADIR, "sed/o5v.sed.ext"))
    sed2 = GalSED("toto", 11)
    sed2.read(os.path.join(TESTDATADIR, "sed/GAL/COSMOS_SED/Ell1_A_0.sed"))
    x2, y2 = sed2.data()
    x1, y1 = sed1.data()
    sed1.sumSpectra(sed2, 1)
    x3, y3 = sed1.data()

    def f1(x):
        return np.interp(x, x1, y1, 0, 0)

    def f2(x):
        return np.interp(x, x2, y2, 0, 0)

    newy3 = f1(x3) + f2(x3)
    assert np.allclose(newy3, y3)
