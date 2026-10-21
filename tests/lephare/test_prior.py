import os

import lephare as lp
import numpy as np


def test_prior():
    # Instantiate a source (Id, gridz)
    src = lp.onesource(1, np.arange(0, 6.1, 0.1))

    bp = [4, 5]
    prior_blue_low = []
    prior_blue_high = []
    prior_red_low = []
    prior_red_high = []

    # loop over the mag to test 19, 21, 23, 25,27, 29, 31
    test_ab = [
        9.1201084e-28,
        1.4454398e-28,
        2.2908677e-29,
        3.6307805e-30,
        5.7543994e-31,
        9.1201084e-32,
        1.4454398e-32,
    ]
    for k in range(0, 7):
        # define the source
        vals = [-99, -99, -99, -99, test_ab[k], -99]
        err_vals = [-99, -99, -99, -99, 0.01, -99]
        src.readsource("1", vals, err_vals, 0, -99, "test")
        src.fltUsed(-1, -1)
        # Compute the N(z) prior
        prior_blue_low.append(src.nzprior(10, 10, 0.1, bp))
        prior_blue_high.append(src.nzprior(10, 10, 4, bp))
        prior_red_low.append(src.nzprior(6, 10, 0.1, bp))
        prior_red_high.append(src.nzprior(6, 10, 4, bp))

    # 1. prior_blue_low decreasing in [0,1]
    assert all(
        prior_blue_low[i] >= prior_blue_low[i + 1] for i in range(len(prior_blue_low) - 1)
    ), "prior_blue_low not decreasing"
    assert all(0 <= val <= 1 for val in prior_blue_low), "prior_blue_low n'est pas dans [0,1]"
    # 2. prior_blue_high croissant et dans [0,1]
    assert all(
        prior_blue_high[i] <= prior_blue_high[i + 1] for i in range(len(prior_blue_high) - 1)
    ), "prior_blue_high not increasing"
    assert all(0 <= val <= 1 for val in prior_blue_high), "prior_blue_high n'est pas dans [0,1]"
    # 3. prior_red_low décroissant et dans [0,1]
    assert all(
        prior_red_low[i] >= prior_red_low[i + 1] for i in range(len(prior_red_low) - 1)
    ), "prior_red_low not decreasingt"
    assert all(0 <= val <= 1 for val in prior_red_low), "prior_red_low n'est pas dans [0,1]"
    # 4. prior_red_high lower than 0.0001
    assert all(val < 0.0001 for val in prior_red_high), "prior_red_high not lower then 0.0001"

    # case with the first filter solution not available. Should switch to next one
    prior_sbc = []
    prior_scd = []
    # loop over the mag to test 21, 23, 25
    test_ab = [1.4454398e-28, 2.2908677e-29, 3.6307805e-30]
    for k in range(0, 3):
        # define the source
        vals = [-99, -99, -99, -99, -99, test_ab[k]]
        err_vals = [-99, -99, -99, -99, -99, 0.01]
        src.readsource("1", vals, err_vals, 0, -99, "test")
        src.fltUsed(-1, -1)
        # Compute the N(z) prior
        # test the comparison between mod Sbc and Scd
        prior_sbc.append(src.nzprior(10, 11.6, 2, bp))
        prior_scd.append(src.nzprior(10, 10.8, 2, bp))

    # The prior should be higher for Sdc galaxies
    for sbc, scd in zip(prior_sbc, prior_scd):
        assert sbc < scd, f"prior_Sbc ({sbc}) not lower than prior_Scd ({scd})"


def test_prior_no_band():
    """No band available for the prior: the weight is 1 whatever the template."""
    src = lp.onesource(1, np.arange(0, 6.1, 0.1))
    bp = [4, 5]
    # nothing observed in the two bands the prior is conditioned on
    src.readsource("1", [-99] * 6, [-99] * 6, 0, -99, "test")
    src.fltUsed(-1, -1)
    for luv, lnir, z in [(10, 10, 0.1), (6, 10, 4.0), (10, 11.6, 2.0)]:
        assert src.nzprior(luv, lnir, z, bp) == 1.0


def test_prior_intermediate_magnitude():
    """20 <= I < 22: the prior is computed, except for z > 2 where it is 0."""
    src = lp.onesource(1, np.arange(0, 6.1, 0.1))
    bp = [4, 5]
    ab = 10 ** (-0.4 * (21.0 + 48.6))
    src.readsource("1", [-99, -99, -99, -99, ab, -99], [-99, -99, -99, -99, 0.01, -99], 0, -99, "test")
    src.fltUsed(-1, -1)
    for luv, lnir in [(6, 10), (10, 11.6), (10, 10.8), (10, 10)]:  # the 4 types
        low = src.nzprior(luv, lnir, 0.5, bp)
        assert np.isfinite(low) and low > 0
        assert src.nzprior(luv, lnir, 3.0, bp) == 0.0


def test_prior_in_fit(tmp_path, monkeypatch):
    """Fit with the N(z) prior enabled, for sources in the three magnitude
    regimes of the prior (I < 20, 20 <= I < 22, I >= 22)."""
    test_dir = os.path.abspath(os.path.dirname(__file__))
    data = os.path.join(test_dir, "..", "data")
    monkeypatch.setenv("LEPHAREDIR", data)
    monkeypatch.setenv("LEPHAREWORK", str(tmp_path))
    for sub in ("filt", "lib_bin", "lib_mag"):
        (tmp_path / sub).mkdir()

    config = lp.read_config(os.path.join(data, "examples/COSMOS.para"))
    fltstr = ",".join(f"lsst/total_{b}.pb" for b in "ugrizy")
    config.update(
        {
            "VERBOSE": "NO",
            "FILTER_LIST": fltstr,
            "FILTER_FILE": "filters_lsst",
            "STAR_SED": "$LEPHAREDIR/sed/STAR/ONE_SED.list",
            "QSO_SED": "$LEPHAREDIR/sed/QSO/ONE_SED.list",
            "GAL_SED": "$LEPHAREDIR/sed/GAL/ONE_SED.list",
            "LIB_ASCII": "NO",
            "CAT_FMT": "MMEE",
            "INP_TYPE": "M",
            "CAT_MAG": "AB",
            "CAT_TYPE": "LONG",
            "GLB_CONTEXT": "-1",
            "AUTO_ADAPT": "NO",
            "Z_STEP": "0.05,0,3",
            "ZFIX": "NO",
            "EB_V": "0.,0.1,0.2",
            "MOD_EXTINC": "1,1",
            "EXTINC_LAW": "SB_calzetti.dat",
            "EM_LINES": "NO",
            "ADD_EMLINES": "0,0",
            "Z_INTERP": "NO",
            "MAG_ABS": "-24,-5",
            "MAG_REF": "3",
            "MABS_METHOD": "0",
            "MABS_CONTEXT": "0",
            "ADDITIONAL_MAG": "none",
            "APPLY_MW_EXTINCTION": "NONE",
            "EXT_MW_CURVE": "CARDELLI",
            "MW_EBV_FILE": "NONE",
            "MW_GLOBAL_EBV": "-1",
            "SPEC_OUT": "NO",
            "PDZ_OUT": "NONE",
            "CAT_OUT": str(tmp_path / "out.dat"),
        }
    )
    lp.prepare(config)

    def run(nz_prior, i_mag):
        config["NZ_PRIOR"] = nz_prior
        photz = lp.PhotoZ(lp.all_types_to_keymap(config))
        # galaxy template at z~1, scaled to the requested magnitude in i
        gal = [g for g in photz.fullLib if isinstance(g, lp.GalSED) and abs(g.red - 1.0) < 1e-6][0]
        mags = np.array(gal.mag) + i_mag - gal.mag[3]
        objs = []
        for k in range(3):
            obj = lp.onesource(k, photz.gridz)
            obj.readsource(str(k), list(mags), [0.05] * 6, 63, -99.0, "")
            photz.prep_data(obj)
            objs.append(obj)
        photz.run_photoz(objs, photz.compute_offsets([]))
        return objs[0]

    for i_mag in (19.0, 21.0, 24.0):
        ref = run("-1", i_mag)
        withprior = run("4,3", i_mag)
        # a solution is found with and without the prior
        assert ref.chimin[0] < 1.0e8
        assert withprior.chimin[0] < 1.0e8
        if i_mag < 20.0:
            # brighter than I=20 the prior forbids z>1, the solution is at z~1
            assert withprior.zmin[0] <= 1.0 + 1e-6
        else:
            # the weight enters the chi2
            assert withprior.chimin[0] != ref.chimin[0]
