import os

import lephare as lp
import pytest


def _base_cosmos_config():
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
            "CAT_IN": str(os.path.expandvars("$LEPHAREWORK/externalz_mag.in")),
            "CAT_FMT": "MMEE",
            "INP_TYPE": "M",
            "CAT_MAG": "AB",
            "CAT_TYPE": "LONG",
            "GLB_CONTEXT": "-1",
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
            "APPLY_MW_EXTINCTION": "NONE",
        }
    )
    return config


def test_read_externalz():
    # PhotoZ::read_externalz (EXTERNALZ_FILE keyword) had zero test coverage:
    # no test in the whole suite even sets this keyword to a real file. It
    # lets the user override the catalogue spec-z from a separate file,
    # matched by source Id, e.g. to fold in a spec-z compilation gathered
    # after the photometric catalogue was built.
    config = _base_cosmos_config()

    # A handful of synthetic sources; the actual magnitudes don't matter here
    # since we are only checking that `zs` gets correctly overridden.
    mag_sources = [
        [24.5493, 23.1701, 22.5265, 22.2859, 22.1366, 22.0255],
        [30.2765, 30.1974, 30.126, 29.6699, 29.4879, 29.4514],
        [23.3172, 22.7789, 22.4013, 21.9102, 21.7947, 21.0578],
    ]
    emag_sources = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    zs_in = [0.65, 0.9, 0.5]

    cat_in = os.path.expandvars(config["CAT_IN"])
    with open(cat_in, "w") as f:
        for idline, (mags, zsin) in enumerate(zip(mag_sources, zs_in), start=1):
            f.write(
                f"{idline} {' '.join(map(str, mags))}  "
                f"{' '.join(map(str, emag_sources))} 63 {zsin} -99 \n"
            )

    # Override the spec-z of sources 1 and 3 only; source 2 is intentionally
    # absent from the external file and must keep its original catalogue zs.
    externalz_file = os.path.expandvars("$LEPHAREWORK/externalz.in")
    externalz_override = {1: 1.234, 3: 0.789}
    with open(externalz_file, "w") as f:
        for src_id, z in externalz_override.items():
            f.write(f"{src_id} {z}\n")
    config["EXTERNALZ_FILE"] = externalz_file

    lp.prepare(config)
    photz = lp.PhotoZ(lp.all_types_to_keymap(config))
    sources = photz.read_photoz_sources()

    assert len(sources) == 3
    by_id = {int(s.spec): s for s in sources}
    assert by_id[1].zs == pytest.approx(externalz_override[1])
    assert by_id[3].zs == pytest.approx(externalz_override[3])
    # source 2 was not listed in the external file: untouched
    assert by_id[2].zs == pytest.approx(zs_in[1])


def test_read_externalz_none_is_a_noop():
    # sanity check: the default "NONE" value must leave zs untouched, exactly
    # as if EXTERNALZ_FILE had never been read
    config = _base_cosmos_config()
    mag_sources = [[24.5493, 23.1701, 22.5265, 22.2859, 22.1366, 22.0255]]
    emag_sources = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    zs_in = [0.65]

    cat_in = os.path.expandvars(config["CAT_IN"])
    with open(cat_in, "w") as f:
        for idline, (mags, zsin) in enumerate(zip(mag_sources, zs_in), start=1):
            f.write(
                f"{idline} {' '.join(map(str, mags))}  "
                f"{' '.join(map(str, emag_sources))} 63 {zsin} -99 \n"
            )
    config["EXTERNALZ_FILE"] = "NONE"

    lp.prepare(config)
    photz = lp.PhotoZ(lp.all_types_to_keymap(config))
    sources = photz.read_photoz_sources()
    assert sources[0].zs == pytest.approx(zs_in[0])
