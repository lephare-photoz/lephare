import os

import numpy as np
import pytest

import lephare as lp


def _build_cosmos_lsst_photoz(z_method):
    """Build a PhotoZ object on the standard LSST/COSMOS library used
    throughout test_fit_onesource.py, with a given Z_METHOD."""
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")

    config_file = os.path.expandvars("$LEPHAREDIR/examples/COSMOS.para")
    config = lp.read_config(config_file)
    fltstr = "lsst/total_u.pb,lsst/total_g.pb,lsst/total_r.pb,lsst/total_i.pb,lsst/total_z.pb,lsst/total_y.pb"
    config.update(
        {
            "VERBOSE": "NO",
            "FILTER_LIST": fltstr,
            "FILTER_FILE": "filters_lsst",
            "STAR_SED": "$LEPHAREDIR/sed/STAR/STAR_MOD_ALL.list",
            "QSO_SED": "$LEPHAREDIR/sed/QSO/SALVATO09/AGN_MOD.list",
            "GAL_SED": "$LEPHAREDIR/sed/GAL/COSMOS_SED/COSMOS_MOD.list",
            "LIB_ASCII": "YES",
            "AUTO_ADAPT": "NO",
            "Z_STEP": "0.05,0,1",
            "ZFIX": "NO",
            "EB_V": "0.,0.1,0.2,0.3",
            "MOD_EXTINC": "0,100",
            "ADD_EMLINES": "0,0",
            "EM_DISPERSION": "1.",
            "ERR_SCALE": " 0.0",
            "ERR_FACTOR": " 1.",
            "Z_INTERP": "NO",
            "MAG_ABS": "-24,-5",
            "MAG_REF": "2",
            "MABS_METHOD": "0",
            "MABS_CONTEXT": "0",
            "Z_METHOD": z_method,
        }
    )
    lp.prepare(config)
    keymap = lp.all_types_to_keymap(config)
    return lp.PhotoZ(keymap)


def _noisy_source(photz):
    # A source with a small, fixed amount of noise added to otherwise
    # perfectly-matching photometry: enough to make the marginalized-PDF
    # median (zgmed[0]) differ measurably from the chi2-minimum (zmin[0]),
    # without changing the chi2-minimum grid point itself, so the two
    # Z_METHOD choices are guaranteed to disagree on this fixture.
    rng = np.random.default_rng(42)
    base_mag = np.array([30.9393, 29.4864, 28.102, 27.1517, 26.8568, 26.6285])
    emag = np.array([0.05] * 6)
    noisy_mag = base_mag + rng.normal(0, emag)

    src = lp.onesource(101, photz.gridz)
    src.readsource("65", list(noisy_mag), list(emag), 0, 0.65, "test")
    src.convertFlux("AB", photz.allFilters)
    photz.prep_data(src)
    a0 = photz.compute_offsets([])
    photz.fit(src, a0)
    photz.fit_uncertainties(src)
    return src


def test_z_method_best_uses_chi2_minimum():
    photz = _build_cosmos_lsst_photoz("BEST")
    src = _noisy_source(photz)
    # sanity check on the fixture: the chi2-minimum sits exactly on the
    # z=0.65 grid point, while the marginalized-PDF median does not
    assert src.zmin[0] == pytest.approx(0.65)
    assert src.zgmed[0] != pytest.approx(0.65, abs=1e-3)

    photz.physical_parameters(src)
    assert src.consiz == pytest.approx(src.zmin[0])


def test_z_method_med_uses_pdf_median():
    # Z_METHOD had zero test coverage anywhere in the suite before this
    # addition: no test ever set this keyword to anything but its "BEST"
    # default.
    photz = _build_cosmos_lsst_photoz("MED")
    src = _noisy_source(photz)
    assert src.zmin[0] == pytest.approx(0.65)
    assert src.zgmed[0] != pytest.approx(0.65, abs=1e-3)

    photz.physical_parameters(src)
    # consiz must now follow the marginalized-PDF median, not the chi2 min
    assert src.consiz == pytest.approx(src.zgmed[0])
    assert src.consiz != pytest.approx(src.zmin[0], abs=1e-3)


def test_z_method_med_and_best_give_different_consiz():
    # Direct A/B comparison on the exact same (seeded, reproducible) noisy
    # source: this is the behavioural difference an end user actually cares
    # about when choosing Z_METHOD.
    photz_best = _build_cosmos_lsst_photoz("BEST")
    src_best = _noisy_source(photz_best)
    photz_best.physical_parameters(src_best)

    photz_med = _build_cosmos_lsst_photoz("MED")
    src_med = _noisy_source(photz_med)
    photz_med.physical_parameters(src_med)

    assert src_best.consiz != pytest.approx(src_med.consiz, abs=1e-3)
