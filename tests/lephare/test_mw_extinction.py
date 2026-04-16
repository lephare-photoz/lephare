import os

import lephare as lp
import numpy as np
import pytest


def test_mw_classic():
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
            "CAT_IN": str(os.path.expandvars("$LEPHAREWORK/mag.in")),
            "CAT_FMT": "MMEE",
            "INP_TYPE": "M",
            "CAT_MAG": "AB",
            "CAT_TYPE": "LONG",
            "GLB_CONTEXT": "-1",
            "AUTO_ADAPT": "YES",
            "ADAPT_LIM": "1.5,26.0",
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
            "APPLY_MW_EXTINCTION": "NONE",
            "EXT_MW_CURVE": "CARDELLI",
            "MW_EBV_FILE": str(os.path.expandvars("$LEPHAREWORK/mw_ebv.in")),
            "MW_GLOBAL_EBV": "-99",
            "MW_REFERENCE_MODEL": "sed/STAR/PICKLES/b5i.sed",
        }
    )

    # Run preparation tasks (libraries)
    lp.prepare(config)
    print("Done reading libraries")

    # Instantiate the photoz object with no MW attenuation
    keymap = lp.all_types_to_keymap(config)
    photz = lp.PhotoZ(keymap)
    print("Instantiate photoz with APPLY_MW_EXTINCTION at none")

    # check that nothing change if NONE
    assert len(photz.mw_classic_extinction_values) == 6
    # No correction proposed for the MW since none
    assert np.isclose(photz.mw_classic_extinction_values, 0.0, atol=0.001).all()

    # Instantiate the photo-z with classic and Cardelli
    config.update(
        {
            "APPLY_MW_EXTINCTION": "CLASSIC",
            "EXT_MW_CURVE": "CARDELLI",
        }
    )
    alambda_cardelli = [4.774, 3.666, 2.694, 2.050, 1.587, 1.304]

    # Instantiate the photoz object
    photz = lp.PhotoZ(lp.all_types_to_keymap(config))
    print("Done instantiate photoz")

    assert len(photz.mw_classic_extinction_values) == 6
    assert np.isclose(photz.mw_classic_extinction_values, alambda_cardelli, atol=0.001).all()
    print("Alambda classic MW correction with Cardelli ", photz.mw_classic_extinction_values)

    # Instantiate the photo-z with classic and Fitzpatrick
    config.update(
        {
            "APPLY_MW_EXTINCTION": "CLASSIC",
            "EXT_MW_CURVE": "LMC_Fitzpatrick.dat",
        }
    )
    alambda_fitzpatrick = [4.814, 3.681, 2.611, 1.962, 1.546, 1.245]

    # Instantiate the photoz object
    photz = lp.PhotoZ(lp.all_types_to_keymap(config))
    print("Done instantiate photoz")

    assert len(photz.mw_classic_extinction_values) == 6
    assert np.isclose(photz.mw_classic_extinction_values, alambda_fitzpatrick, atol=0.001).all()
    print("Alambda classic MW correction with LMC Fitzpatrick ", photz.mw_classic_extinction_values)

    # Create the input ascii file
    mag_sources = [
        [30.9393, 29.4864, 28.102, 27.1517, 26.8568, 26.6285],  # same test as one source z=0.65
        [24.5493, 23.1701, 22.5265, 22.2859, 22.1366, 22.0255],  # mod 1, no attenuation, z=0.1
        [30.2765, 30.1974, 30.126, 29.6699, 29.4879, 29.4514],  # mod 30, ebv=0.2, z=0.9
        [0.656911, -0.0506009, 0.148528, 0.357171, 0.483391, 0.548042],  # star, mod 24
        [23.3172, 22.7789, 22.4013, 21.9102, 21.7947, 21.0578],  # ebv=0.3, z=0.5, pl_TQSO1_template_norm.sed
    ]
    print(mag_sources[0])
    emag_sources = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    mw_ebv_sources = [0.01, 0.05, 0.1, 0.2, 0.3]
    zs_in = [0.65, 0.1, 0.9, 0.000, 0.5]

    # Temporary input file with the observed magnitudes
    fil = os.path.expandvars("$LEPHAREWORK/mag.in")
    with open(fil, "w") as f:
        for idline in range(1, 6):
            mags = mag_sources[idline - 1][0:6]
            zsin = zs_in[idline - 1]
            line = (
                f"{idline} {' '.join(map(str, mags))}  "
                f"{' '.join(map(str, emag_sources))} 63 "
                f"{str(zsin)} -99 \n"
            )
            f.write(line)

    # Temporary input file with the MW EBV
    fil = os.path.expandvars("$LEPHAREWORK/mw_ebv.in")
    print(fil)
    with open(fil, "w") as f:
        for idline in range(1, 6):
            line = f"{idline}  " f"{str(mw_ebv_sources[idline - 1])} \n"
            f.write(line)

    # Compute offsets depending on the AUTO_ADAPT and APPLY_SYSSHIFT options (0 if none)
    adapt_srcs = photz.read_autoadapt_sources()
    # Only 2 sources satisfy the mag/redshift criteria
    assert len(adapt_srcs) == 2
    photz.read_mw_ebv(adapt_srcs)
    a0 = photz.compute_offsets(adapt_srcs)
    assert len(a0) == 6

    # Check that the observed magnitudes have been correctly changed in the auto-adapt
    # First source
    assert adapt_srcs[0].mw_ebv == pytest.approx(0.05)
    expected_corr = -1.0 * np.array(alambda_fitzpatrick) * 0.05
    applied_corr = adapt_srcs[0].mab - np.array([24.5493, 23.1701, 22.5265, 22.2859, 22.1366, 22.0255])
    assert np.isclose(expected_corr, applied_corr, atol=0.01).all()
    # Second source
    assert adapt_srcs[1].mw_ebv == pytest.approx(0.3)
    expected_corr = -1.0 * np.array(alambda_fitzpatrick) * 0.3
    applied_corr = adapt_srcs[1].mab - np.array([23.3172, 22.7789, 22.4013, 21.9102, 21.7947, 21.0578])
    assert np.isclose(expected_corr, applied_corr, atol=0.01).all()

    # read ascii table with sources
    allsources = photz.read_photoz_sources()
    photz.read_mw_ebv(allsources)
    print("Done creating sources")
    # Check that input file is well read
    for k in range(0, 5):
        # check the the MW EBV read by the code corresponds to the input one
        assert allsources[k].mw_ebv == pytest.approx(mw_ebv_sources[k])

    # run the fit
    photz.run_photoz(allsources, [])
    print("Done with fit")

    # Check that the observed magnitudes have been correctly changed
    for k in range(0, 5):
        # Check magnitudes
        expected_corr = -1.0 * np.array(alambda_fitzpatrick) * mw_ebv_sources[k]
        print("expected corr ", expected_corr)
        applied_corr = allsources[k].mab - np.array(mag_sources[k])
        assert np.isclose(expected_corr, applied_corr, atol=0.01).all()


def test_mw_galametz_single_ebv():
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
            "CAT_IN": str(os.path.expandvars("$LEPHAREWORK/mag.in")),
            "CAT_FMT": "MMEE",
            "INP_TYPE": "M",
            "CAT_MAG": "AB",
            "CAT_TYPE": "LONG",
            "GLB_CONTEXT": "-1",
            "AUTO_ADAPT": "YES",
            "ADAPT_LIM": "1.5,26.0",
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
            "APPLY_MW_EXTINCTION": "NONE",
            "EXT_MW_CURVE": "LMC_Fitzpatrick.dat",
            "MW_EBV_FILE": str(os.path.expandvars("$LEPHAREWORK/mw_ebv.in")),
            "MW_GLOBAL_EBV": "0.1",
            "MW_REFERENCE_MODEL": "sed/STAR/PICKLES/b5i.sed",
        }
    )

    # Run preparation tasks (libraries)
    lp.prepare(config)
    print("Done reading libraries")

    print("Instantiate photoz with no MW attenuation")
    photz = lp.PhotoZ(lp.all_types_to_keymap(config))

    no_mw_flux = photz.flux[0][0]

    print("Instantiate photoz with GALAMETZ + single MW EBV")
    config.update(
        {
            "APPLY_MW_EXTINCTION": "GALAMETZ",
            "EXT_MW_CURVE": "LMC_Fitzpatrick.dat",
            "MW_GLOBAL_EBV": "0.1",
        }
    )
    lp.prepare(config)
    print("Done reading libraries")
    photz = lp.PhotoZ(lp.all_types_to_keymap(config))

    # Value taken for the first QSO, first redshift and band (first item of the library)
    mw_flux = photz.flux[0][0]
    assert no_mw_flux == pytest.approx(10 ** (0.4 * (4.7932 / 0.7906 * 0.9099) * 0.1) * mw_flux)


def test_mw_galametz_multiple_ebv():
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
            "CAT_IN": str(os.path.expandvars("$LEPHAREWORK/mag.in")),
            "CAT_FMT": "MMEE",
            "INP_TYPE": "M",
            "CAT_MAG": "AB",
            "CAT_TYPE": "LONG",
            "GLB_CONTEXT": "-1",
            "AUTO_ADAPT": "YES",
            "ADAPT_LIM": "1.5,26.0",
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
            "APPLY_MW_EXTINCTION": "NONE",
            "EXT_MW_CURVE": "LMC_Fitzpatrick.dat",
            "MW_GLOBAL_EBV": "-1.",
            "MW_EBV_FILE": str(os.path.expandvars("$LEPHAREWORK/mw_ebv.in")),
            "MW_REFERENCE_MODEL": "sed/STAR/PICKLES/b5i.sed",
        }
    )

    # Run preparation tasks (libraries)
    lp.prepare(config)
    print("Done reading libraries")

    print("Instantiate photoz with no MW attenuation")
    photz = lp.PhotoZ(lp.all_types_to_keymap(config))

    no_mw_flux = photz.flux[0][0]

    print("Instantiate photoz with GALAMETZ + multiple MW EBV")
    # Create the input ascii file
    mag_sources = [
        [30.9393, 29.4864, 28.102, 27.1517, 26.8568, 26.6285],  # same test as one source z=0.65
        [24.5493, 23.1701, 22.5265, 22.2859, 22.1366, 22.0255],  # mod 1, no attenuation, z=0.1
        [30.2765, 30.1974, 30.126, 29.6699, 29.4879, 29.4514],  # mod 30, ebv=0.2, z=0.9
        [0.656911, -0.0506009, 0.148528, 0.357171, 0.483391, 0.548042],  # star, mod 24
        [23.3172, 22.7789, 22.4013, 21.9102, 21.7947, 21.0578],  # ebv=0.3, z=0.5, pl_TQSO1_template_norm.sed
    ]
    print(mag_sources[0])
    emag_sources = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    mw_ebv_sources = [0.01, 0.05, 0.1, 0.2, 0.3]
    zs_in = [0.65, 0.1, 0.9, 0.000, 0.5]

    # Temporary input file with the observed magnitudes
    fil = os.path.expandvars("$LEPHAREWORK/mag.in")
    with open(fil, "w") as f:
        for idline in range(1, 6):
            mags = mag_sources[idline - 1][0:6]
            zsin = zs_in[idline - 1]
            line = (
                f"{idline} {' '.join(map(str, mags))}  "
                f"{' '.join(map(str, emag_sources))} 63 "
                f"{str(zsin)} -99 \n"
            )
            f.write(line)

    # Temporary input file with the MW EBV
    fil = os.path.expandvars("$LEPHAREWORK/mw_ebv.in")
    print(fil)
    with open(fil, "w") as f:
        for idline in range(1, 6):
            line = f"{idline}  " f"{str(mw_ebv_sources[idline - 1])} \n"
            f.write(line)

    # create the library with the right MW extnction option for Galametz
    config.update(
        {
            "APPLY_MW_EXTINCTION": "GALAMETZ",
            "EXT_MW_CURVE": "LMC_Fitzpatrick.dat",
            "MW_GLOBAL_EBV": "-99",
            "MW_EBV_FILE": str(os.path.expandvars("$LEPHAREWORK/mw_ebv.in")),
        }
    )
    lp.prepare(config)

    print("Done reading libraries")
    photz = lp.PhotoZ(lp.all_types_to_keymap(config))

    # Value taken for the first QSO, first redshift and band (first item of the library)
    # value don't change when multiple MW EBV
    mw_flux = photz.flux[0][0]
    assert (no_mw_flux / mw_flux) == pytest.approx(1, abs=1e-3)

    # read ascii table with sources
    allsources = photz.read_photoz_sources()
    photz.read_mw_ebv(allsources)
    print("Done creating sources")
    # Check that input file is well read
    for k in range(0, 5):
        # check the the MW EBV read by the code corresponds to the input one
        assert allsources[k].mw_ebv == pytest.approx(mw_ebv_sources[k])

    # run the fit
    photz.run_photoz(allsources, [])
    print("Done with fit")

    # Value taken for the first QSO, first redshift and band (first item of the library)
    # flux should contain the value for the last mw ebv
    mw_flux = photz.flux[0][0]
    assert (no_mw_flux / mw_flux) == pytest.approx(10 ** (0.4 * (4.7932 / 0.7906 * 0.9099) * 0.3), abs=1e-3)
