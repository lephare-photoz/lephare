import os
import time

import lephare as lp
import numpy as np
import pytest


def test_photoz_cosmos():
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    print(test_dir)

    # Read the config file.
    config_file = os.path.expandvars("$LEPHAREDIR/examples/COSMOS.para")
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
            "INP_TYPE": "M",
            "CAT_MAG": "AB",
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
            "ADDITIONAL_MAG": "filters_lsst",
            "SPEC_OUT": str(os.path.expandvars("$LEPHAREWORK/spec")),
        }
    )
    print("SPEC ", os.path.expandvars("$LEPHAREWORK/spec"))
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

    allsources = []
    mag_sources = [
        [30.9393, 29.4864, 28.102, 27.1517, 26.8568, 26.6285],  # same test as one source z=0.65
        [24.5493, 23.1701, 22.5265, 22.2859, 22.1366, 22.0255],  # mod 1, no attenuation, z=0.1
        [30.2765, 30.1974, 30.126, 29.6699, 29.4879, 29.4514],  # mod 30, ebv=0.2, z=0.9
        [0.656911, -0.0506009, 0.148528, 0.357171, 0.483391, 0.548042],  # star, mod 24
        [23.3172, 22.7789, 22.4013, 21.9102, 21.7947, 21.0578],
    ]  # ebv=0.3, z=0.5, pl_TQSO1_template_norm.sed
    print(mag_sources[0])
    for k in range(0, 5):
        # Instantiate a source (Id, gridz)
        src = lp.onesource(101, photz.gridz)
        # read the source, change Id, attribute mag/err, ...
        cont = 0
        src.readsource(
            str(k),
            mag_sources[k],
            [0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
            cont,
            -99,
            "test",
        )
        allsources.append(src)
    print("Done creating source")

    photz.prep_data(allsources)

    a0 = photz.compute_offsets([])
    src.adapt_mag(a0)
    print("Done with offsets")
    assert len(a0) == 6

    photz.run_photoz(allsources, a0)
    print("Done with fit")

    # Check library
    assert (allsources[0].chimin[0] < allsources[0].chimin[1]) and (
        allsources[0].chimin[0] < allsources[0].chimin[2]
    )
    assert (allsources[1].chimin[0] < allsources[1].chimin[1]) and (
        allsources[1].chimin[0] < allsources[1].chimin[2]
    )
    assert (allsources[2].chimin[0] < allsources[2].chimin[1]) and (
        allsources[2].chimin[0] < allsources[2].chimin[2]
    )
    assert (allsources[3].chimin[2] < allsources[3].chimin[0]) and (
        allsources[3].chimin[2] < allsources[3].chimin[1]
    )
    assert (allsources[4].chimin[1] < allsources[4].chimin[0]) and (
        allsources[4].chimin[1] < allsources[4].chimin[0]
    )
    assert allsources[0].chimin[0] < 1.0e-3
    assert allsources[1].chimin[0] < 1.0e-3
    assert allsources[2].chimin[0] < 1.0e-3
    assert allsources[3].chimin[2] < 1.0e-3
    assert allsources[4].chimin[1] < 1.0e-3

    # Check minimisation source 1, same test as one source z=0.65
    assert allsources[0].zmin[0] == pytest.approx(0.65)
    assert allsources[0].dmmin[0] == pytest.approx(1.0000, abs=1e-04)
    assert allsources[0].imasmin[0] == pytest.approx(1)
    # Check minimisation source 2, mod 1, no attenuation, z=0.1
    assert allsources[1].zmin[0] == pytest.approx(0.1)
    assert allsources[1].dmmin[0] == pytest.approx(1.0000, abs=1e-04)
    assert allsources[1].imasmin[0] == pytest.approx(1)
    # Check minimisation source 3, mod 30, ebv=0.2, z=0.9
    assert allsources[2].zmin[0] == pytest.approx(0.9)
    assert allsources[2].dmmin[0] == pytest.approx(1.0000, abs=1e-04)
    assert allsources[2].imasmin[0] == pytest.approx(30)
    # Check minimisation source 4,# star, mod 24
    assert allsources[3].dmmin[2] == pytest.approx(1.0000, abs=1e-04)
    assert allsources[3].imasmin[2] == pytest.approx(24)
    # Check minimisation source 5, ebv=0.3, z=0.5, pl_TQSO1_template_norm.sed
    assert allsources[4].zmin[1] == pytest.approx(0.5)
    assert allsources[4].dmmin[1] == pytest.approx(1.0000, abs=1e-04)
    assert allsources[4].imasmin[1] == pytest.approx(30)

    # Check les mag predites
    # Since we use the same filters as additional, we should have the input magnitude
    # True only for galaxies
    for k in range(0, 3):
        assert np.isclose(allsources[k].magPred, mag_sources[k], atol=0.01).all()

    ti1 = int(time.time())
    photz.write_outputs(allsources, ti1)
    print("Done with writting")

    for k in range(0, 5):
        file_test = os.path.expandvars("$LEPHAREWORK/spec/") + "Id" + str(k) + ".spec"
        assert os.path.exists(file_test)
