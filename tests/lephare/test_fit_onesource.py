import os

import lephare as lp
import numpy as np


def test_onefit():
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    print(test_dir)

    # Read the config file.
    config_file = os.path.expandvars("$LEPHAREDIR/examples/COSMOS.para")
    # config_file = os.path.join(test_data_dir, "examples/COSMOS.para")
    config = lp.read_config(config_file)
    print(config_file)
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
        }
    )
    keymap = lp.all_types_to_keymap(config)

    # Run preparation tasks (libraries)
    lp.prepare(config)
    print("Done reading libraries")
    assert os.path.exists(os.path.expandvars("$LEPHAREWORK/lib_mag/CE_COSMOS.doc"))

    # Instantiate the photoz object
    photz = lp.PhotoZ(keymap)
    print("Done instantiate photoz")
    assert len(photz.gridz) == 21
    assert photz.gridz[0] == 0
    assert photz.gridz[20] == 1

    # Instantiate a source (Id, gridz)
    src = lp.onesource(101, photz.gridz)
    # read the source, change Id, attribute mag/err, ...
    cont = 0
    src.readsource(
        "65",
        [30.9393, 29.4864, 28.102, 27.1517, 26.8568, 26.6285],
        [0.01, 0.01, 0.01, 0.01, 0.01],
        cont,
        0.65,
        "test",
    )
    src.convertFlux("AB", photz.allFilters)
    print("Done creating source")
    assert src.spec == "65"
    assert (
        np.testing.assert_almost_equal(
            -2.5 * np.log10(src.ab) - 48.60, [30.9393, 29.4864, 28.102, 27.1517, 26.8568, 26.6285]
        )
        is None
    )

    photz.prep_data(src)

    a0 = photz.compute_offsets([])
    src.adapt_mag(a0)
    print("Done with offsets")
    assert len(a0) == 6

    photz.fit_onesource(src)
    print("Done with fit")
    assert np.isclose(src.zmin[0], 0.65)
    assert np.isclose(src.dmmin[0], 1.0000, atol=1e-04)
    assert np.isclose(src.imasmin[0], 1)
    assert src.chimin[0] < src.chimin[1] < src.chimin[2]

    photz.uncertainties_onesource(src)
    print("Done with uncertainties")
    max_position = np.argmax(src.pdfmap[11].vPDF)
    assert max_position == 13
    assert len(src.pdfmap[11].xaxis) == 21
    assert np.isclose(src.zgmed[0], 0.65, atol=1e-02)
    assert np.isclose(src.zgmode[0], 0.65, atol=1e-02)
    assert np.isclose(src.zgmin[0], 0.65, atol=1e-02)
    assert src.zgmin[3] < src.zgmin[1] < src.zgmin[0] < src.zgmin[2] < src.zgmin[4]
    assert src.zgmed[3] < src.zgmed[1] < src.zgmed[0] < src.zgmed[2] < src.zgmed[4]
    assert src.zgmode[3] < src.zgmode[1] < src.zgmode[0] < src.zgmode[2] < src.zgmode[4]

    photz.physpara_onesource(src)
    print("Done with physical parameters")
    assert np.isclose(src.mabs[5], 26.6285 - 42.9504 - src.kap[5], atol=3e-2)
    assert np.isclose(src.mabs[0], 30.9393 - 42.9504 - src.kap[0], atol=3e-2)
    assert np.testing.assert_almost_equal(src.absfilt, [0, 1, 2, 3, 4, 5]) is None

    minl = 3000.0
    maxl = 13000.0
    gal1 = photz.besttemplate_onesource(src, 0, minl, maxl)
    gal2 = photz.besttemplate_onesource(src, 1, minl, maxl)
    fir = photz.besttemplate_onesource(src, 2, minl, maxl)
    qso = photz.besttemplate_onesource(src, 3, minl, maxl)
    star = photz.besttemplate_onesource(src, 4, minl, maxl)
    print("Done with templates generation")
    assert len(gal2[0]) == 0
    assert len(fir[0]) == 0
    assert all(2999.0 < x < 13001.0 for x in gal1[0])
    assert all(10.0 < x < 40.0 for x in gal1[1])
    assert all(2999.0 < x < 13001.0 for x in qso[0])
    assert all(10.0 < x < 40.0 for x in qso[1])
    assert all(2999.0 < x < 13001.0 for x in star[0])
    assert all(10.0 < x < 40.0 for x in star[1])
