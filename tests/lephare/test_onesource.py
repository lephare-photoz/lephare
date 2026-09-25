import numpy as np
import pytest
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
    src = onesource()
    vals = [30.9393, 29.4864, 28.102, 27.1517, 26.8568, 26.6285]
    err_vals = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    err_vals_wrong = [0.01, 0.01, 0.01, 0.01, 0.01]
    with pytest.raises(ValueError):
        src.readsource("65", vals, err_vals_wrong, 0, 0.65, "test")

    src.readsource("65", vals, err_vals, 1, 0.65, "test")
    assert src.ab == vals
    assert src.sab == err_vals
    assert src.spec == "65"
    assert src.cont == 1
    assert src.zs == 0.65
    assert src.str_inp == "test"


def test_readsource2():
    # Instantiate a source
    src = onesource(101, [0, 0.1, 1])
    # read the source, change Id, attribute flux/err, ...
    src.readsource(
        "10", [-6.414e-32, 1.3182e-31, 1.6905e-31], [1.1022e-31, 9.8579e-32, 5.8665e-32], 6, 2.1, "add"
    )
    assert src.spec == "10"
    assert src.pdfmap[9].size() == 11
    assert src.zs == pytest.approx(2.1)
    assert src.cont == 6
    assert src.ab[0] * 1.0e32 == pytest.approx(-6.414)
    assert src.sab[0] * 1.0e31 == pytest.approx(1.1022)


def test_fltused():
    src = onesource(101, [0, 0.1, 1])
    # upper limit in band 3 using negative error
    src.readsource(
        "10", [-6.414e-32, 1.3182e-31, 1.6905e-31], [1.1022e-31, 9.8579e-32, -5.8665e-32], 7, 2.1, "add"
    )
    # Test without global or forbiden context
    src.fltUsed(-1, -1)
    assert src.cont == 7
    assert np.array_equal(src.busnorma, [1, 1, 0])
    assert np.array_equal(src.busul, [0, 0, 1])
    assert src.nbused == 3
    assert src.nbul == 1
    # Test with a forbidden context removing the first band
    src.fltUsed(-1, 1)
    assert src.cont == 7
    assert np.array_equal(src.busnorma, [0, 1, 0])
    assert np.array_equal(src.busul, [0, 0, 1])
    assert src.nbused == 2
    assert src.nbul == 1
    # Test with a global context using only the second band
    src.fltUsed(2, -1)
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


def _make_src_for_error_rescaling():
    src = onesource(101, [0, 0.1, 1])
    ab = [10.0, 20.0, 30.0]
    sab = [1.0, 2.0, 3.0]  # 10% relative flux error in every band
    src.readsource("1", ab, sab, 7, 0.5, "test")
    return src, np.array(ab), np.array(sab)


def test_rescale_flux_errors_per_band():
    # one min_err/fac_err value per band: previously the only branch tested
    src, ab, sab = _make_src_for_error_rescaling()
    min_err = np.array([0.05, 0.02, 0.0])
    fac_err = np.array([1.0, 2.0, 0.5])
    src.rescale_flux_errors(list(min_err), list(fac_err))

    # replicate the documented formula (fractional error added in quadrature
    # to min_err, in magnitude space, then rescaled back to flux and by
    # fac_err)
    frac_err = 1.086 * sab / ab
    frac_err = np.sqrt(frac_err**2 + min_err**2)
    expected = np.abs(ab) * frac_err / 1.086
    expected = expected * fac_err
    assert np.allclose(src.sab, expected)


def test_rescale_flux_errors_scalar():
    # a single-element min_err/fac_err is meant to apply to every band alike;
    # this branch had zero test coverage before this addition
    src, ab, sab = _make_src_for_error_rescaling()
    min_err_scalar = 0.03
    fac_err_scalar = 1.5
    src.rescale_flux_errors([min_err_scalar], [fac_err_scalar])

    frac_err = 1.086 * sab / ab
    frac_err = np.sqrt(frac_err**2 + min_err_scalar**2)
    expected = np.abs(ab) * frac_err / 1.086
    expected = expected * fac_err_scalar
    assert np.allclose(src.sab, expected)

    # a scalar and a uniform per-band array with the same value must agree
    src2, _, _ = _make_src_for_error_rescaling()
    src2.rescale_flux_errors([min_err_scalar, min_err_scalar, min_err_scalar], [fac_err_scalar] * 3)
    assert np.allclose(src.sab, src2.sab)


def test_rescale_flux_errors_size_mismatch_is_noop():
    # neither a single value nor one per band: lephare can't apply the
    # correction and must leave sab untouched (and not crash) rather than
    # silently guessing; this is the kind of defensive branch that matters
    # because a silent wrong-size mismatch would otherwise corrupt fits
    src, ab, sab = _make_src_for_error_rescaling()
    src.rescale_flux_errors([0.05, 0.05], [1.0, 1.0])
    assert np.allclose(src.sab, sab)
