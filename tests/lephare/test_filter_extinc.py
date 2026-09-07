import os

import lephare as lp
import numpy as np
import pytest
import scipy.integrate as sciint


def test_filter_extinc():
    """Test the filter_extinc runner makes a file and has good values"""
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    out_file = os.path.join(test_dir, "../tmp", "filter_extinc.dat")
    options = {
        "FILTER_FILE": os.path.join(test_dir, "../data", "filt", "LSST_FILTERS.dat"),
        "EXT_ATMOS_CURVE": "SB_calzetti.dat",
        "EXT_MW_CURVE": "CARDELLI",
        "OUTPUT": out_file,
    }
    runner = lp.FiltExt(config_keymap=lp.all_types_to_keymap(lp.default_cosmos_config), **options)
    all_filters, aint, albdav, albd = runner.run()

    # Check it made the file
    assert os.path.exists(out_file)
    with open(out_file, "r") as f:
        contents = f.read()
    assert float(contents.split()[-1]) == pytest.approx(1.3038307748211582)

    # check computation
    atmoext = lp.ext("atmo", 0)
    atmoext.read(os.path.join(os.environ["LEPHAREDIR"], "ext", "SB_calzetti.dat"))
    x1 = [o.lamb for o in atmoext.lamb_ext]
    y1 = [o.val for o in atmoext.lamb_ext]

    def ext(x):
        return np.interp(x, x1, y1, 0, 0)  # 0, 0 are the extrapolated default values

    aint2 = []
    # errs = []
    for f in all_filters:
        x2 = [o.lamb for o in f.lamb_trans]
        y2 = [o.val for o in f.lamb_trans]

        def filt(x):
            # 0, 0 are the extrapolated default values
            return np.interp(x, x2, y2, 0, 0)  # noqa: B023

        def fe(x):
            return ext(x) * filt(x)

        n = sciint.quad(fe, f.lmin(), f.lmax(), limit=200, epsabs=1.0e-3, epsrel=1.0e-4)
        d = sciint.quad(filt, f.lmin(), f.lmax(), limit=200, epsabs=1.0e-3, epsrel=1.0e-4)
        res = n[0] / d[0]
        aint2.append(res)
        # errs.append(sqrt(res**2 *(n[1]**2/n[0]**2 + d[1]**2/d[0]**2)) )

    np.testing.assert_array_almost_equal(aint, aint2, 5.0e-4)


@pytest.fixture
def lsst_filters(set_env_vars):
    return os.path.join(os.environ["LEPHAREDIR"], "filt", "LSST_FILTERS.dat")


def test_atmospheric_curve_none(lsst_filters):
    """With EXT_ATMOS_CURVE=NONE the atmospheric extinction is flagged as 99."""
    all_filters, aint, albdav, albd = lp.calculate_extinction_values(lsst_filters, "NONE", "CARDELLI")

    assert len(all_filters) == 6
    # 99 is the "not computed" sentinel, one per filter
    assert aint == [99.0] * len(all_filters)
    # The galactic values are still computed
    np.testing.assert_array_less(0, albdav)
    np.testing.assert_allclose(albd, np.array(albdav) * 3.1)


def test_cardelli_matches_per_filter_call(lsst_filters):
    """The CARDELLI branch is just cardelli_ext applied filter by filter."""
    all_filters, _, albdav, albd = lp.calculate_extinction_values(lsst_filters, "NONE", "CARDELLI")

    expected = np.array([lp.cardelli_ext(f) for f in all_filters])
    np.testing.assert_allclose(albdav, expected)
    # A(lbd)/E(B-V) = Rv * A(lbd)/Av, with Rv = 3.1 hardcoded for Cardelli
    np.testing.assert_allclose(albd, expected * 3.1)
    # Extinction must decrease from u through y
    np.testing.assert_array_less(albdav[1:], albdav[:-1])


@pytest.mark.parametrize(
    "galec,rv",
    [
        ("SB_calzetti.dat", 4.05),
        ("SMC_prevot.dat", 2.72),
        ("MW_seaton.dat", 3.1),
        ("LMC_Fitzpatrick.dat", 3.1),
    ],
)
def test_tabulated_galactic_curve_rv(lsst_filters, galec, rv, capsys):
    """A tabulated MW curve is given in k(lbd), and Rv depends on the law."""
    all_filters, _, albdav, albd = lp.calculate_extinction_values(lsst_filters, "NONE", galec, verbose=True)

    # The chosen Rv is reported when verbose
    assert f"assuming Rv={rv}" in capsys.readouterr().out
    # Tabulated curves give A(lbd)/E(B-V) directly; A(lbd)/Av divides by Rv
    np.testing.assert_allclose(albdav, np.array(albd) / rv)
    np.testing.assert_array_less(0, albd)
    # Cross-check against the underlying per-filter convolution
    galactic_ext = lp.ext(galec, 1)
    galactic_ext.read(os.path.join(os.environ["LEPHAREDIR"], "ext", galec))
    expected = [lp.compute_filter_extinction(f, galactic_ext) for f in all_filters]
    np.testing.assert_allclose(albd, expected)


def test_absolute_curve_paths(lsst_filters):
    """Absolute curve paths are used as given, not resolved against LEPHAREDIR."""
    ext_dir = os.path.join(os.environ["LEPHAREDIR"], "ext")
    relative = lp.calculate_extinction_values(lsst_filters, "SB_calzetti.dat", "SMC_prevot.dat")
    absolute = lp.calculate_extinction_values(
        lsst_filters,
        os.path.join(ext_dir, "SB_calzetti.dat"),
        os.path.join(ext_dir, "SMC_prevot.dat"),
    )
    np.testing.assert_allclose(relative[1], absolute[1])
    np.testing.assert_allclose(relative[3], absolute[3])


def test_classic_extinction_values(lsst_filters):
    """classic_extinction_values pulls the same parameters out of a config keymap."""
    config = {
        "FILTER_FILE": lsst_filters,
        "EXT_ATMOS_CURVE": "SB_calzetti.dat",
        "EXT_MW_CURVE": "CARDELLI",
    }
    all_filters, aint, albdav, albd = lp.classic_extinction_values(config)

    direct = lp.calculate_extinction_values(lsst_filters, "SB_calzetti.dat", "CARDELLI")
    assert len(all_filters) == len(direct[0])
    np.testing.assert_allclose(aint, direct[1])
    np.testing.assert_allclose(albdav, direct[2])
    np.testing.assert_allclose(albd, direct[3])


def test_classic_extinction_values_defaults(lsst_filters):
    """Empty curve keywords fall back to NONE (atmospheric) and CARDELLI (galactic)."""
    config = {"FILTER_FILE": lsst_filters, "EXT_ATMOS_CURVE": "", "EXT_MW_CURVE": ""}
    _, aint, albdav, albd = lp.classic_extinction_values(config)

    assert aint == [99.0] * 6
    np.testing.assert_allclose(albd, np.array(albdav) * 3.1)
