import os

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
        [0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
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
    print("Done with the spectra")


def test_rm_discrepant():
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
            "LIB_ASCII": "NO",
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
            "RM_DISCREPANT_BD": "500",
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
        [30.9393, 29.4864, 28.102, 22.1517, 26.8568, 26.6285],  # Add fake offset to the i-band
        [0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        cont,
        0.65,
        "test",
    )
    src.convertFlux("AB", photz.allFilters)
    print("Done creating source")
    assert src.spec == "65"
    assert (
        np.testing.assert_almost_equal(
            -2.5 * np.log10(src.ab) - 48.60, [30.9393, 29.4864, 28.102, 22.1517, 26.8568, 26.6285]
        )
        is None
    )

    photz.prep_data(src)
    photz.fit_onesource(src)
    assert src.nbused == 5
    assert np.isclose(src.zmin[0], 0.65)
    assert np.isclose(src.dmmin[0], 1.0000, atol=1e-04)
    assert np.isclose(src.imasmin[0], 1)


def test_physicalpara_bc03():
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    print(test_dir)

    # Read the config file.
    config_file = os.path.expandvars("$LEPHAREDIR/examples/COSMOS.para")
    config = lp.read_config(config_file)
    print(config_file)
    lsst_files = [
        "lsst/total_u.pb",
        "lsst/total_g.pb",
        "lsst/total_r.pb",
        "lsst/total_i.pb",
        "lsst/total_z.pb",
        "lsst/total_y.pb",
    ]
    spitzer_files = ["spitzer/irac_3.pb", "spitzer/irac_4.pb"]
    herschel_files = ["herschel/PACS_100.pb", "herschel/PACS_160.pb"]
    fltstr = ",".join(lsst_files + spitzer_files + herschel_files)
    config.update(
        {
            "VERBOSE": "NO",
            "STAR_SED": "$LEPHAREDIR/sed/STAR/STAR_MOD_ALL.list",
            "QSO_SED": "$LEPHAREDIR/sed/QSO/SALVATO09/AGN_MOD.list",
            "FILTER_LIST": fltstr,
            "FILTER_FILE": "filters_fir",
            "FILTER_CALIB": "0,0,0,0,0,0,1,1,1,1",
            "ZPHOTLIB": "BC03_IR",
            "AUTO_ADAPT": "NO",
            "Z_STEP": "0.2,0,1",
            "ZFIX": "YES",
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

    # Filters
    filterlib = lp.Filter(config_keymap=keymap)
    filterlib.run()

    # Generate a test BC03 template list
    pathlist_bc03 = os.path.expandvars("$LEPHAREWORK/testBC03_MOD.list")
    print("BC03 file ", pathlist_bc03)
    lines = [
        "BC03_CHAB/bc2003_lr_m62_chab_tau1_dust00.ised_ASCII  BC03\n",
        "BC03_CHAB/bc2003_lr_m62_chab_tau15_dust00.ised_ASCII  BC03\n",
    ]
    with open(os.path.expandvars(pathlist_bc03), "w") as file:
        file.writelines(lines)
    # Generate a test BC03 age list
    pathage_bc03 = os.path.expandvars("$LEPHAREWORK/agesBC03.txt")
    print("BC03 ages ", pathage_bc03)
    lines = ["0.360203008\n", "2\n", "7\n"]
    with open(os.path.expandvars(pathage_bc03), "w") as file:
        file.writelines(lines)

    # Download large files if needed
    lp.data_retrieval.get_auxiliary_data(
        keymap=lp.default_cosmos_config.copy(),
        additional_files=[
            "sed/GAL/BC03_CHAB/bc2003_lr_m62_chab_tau1_dust00.ised_ASCII",
            "sed/GAL/BC03_CHAB/bc2003_lr_m62_chab_tau15_dust00.ised_ASCII",
        ],
    )

    sedlib = lp.Sedtolib(config_keymap=keymap)
    sedlib.run(
        typ="GAL",
        gal_lib="LIB_BC03",
        gal_sed=pathlist_bc03,
        sel_age=pathage_bc03,
    )

    maglib = lp.MagGal(config_keymap=keymap)
    maglib.run(
        typ="GAL",
        lib_ascii="YES",
        gal_lib_in="LIB_BC03",
        gal_lib_out="BC03_IR",
        eb_v="0.,0.25,0.5",
        mod_extinc="0,12",
        extinc_law="SB_calzetti.dat",
        em_lines="PHYS",
        em_dispersion="1.",
    )

    print("Done creating libraries")
    assert os.path.exists(os.path.expandvars("$LEPHAREWORK/lib_mag/BC03_IR.doc"))

    # Instantiate the photoz object
    photz = lp.PhotoZ(keymap)
    print("Done instantiate photoz")
    assert len(photz.gridz) == 6
    assert photz.gridz[0] == 0
    assert photz.gridz[5] == 1

    # Instantiate a source (Id, gridz)
    src = lp.onesource(101, photz.gridz)
    # read the source, change Id, attribute mag/err, ...
    cont = 0
    # put a crazy magnitude to test rm_band
    maglib_bc03 = [51.9275, 51.2836, 50.39, 49.9945, 49.7161, 49.4793, 49.0574, 49.5759, 55.004, 56.0812]
    # Substract 25 to get a mass of 10^10 Msol
    maglib_bc03 = [x - 25 for x in maglib_bc03]
    src.readsource(
        "40",
        maglib_bc03,
        [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        cont,
        0.4,
        "ebv0.5_DM41.6845age3.60203e8",
    )

    src.convertFlux("AB", photz.allFilters)
    print("Done creating source")
    assert src.spec == "40"

    photz.prep_data(src)

    a0 = photz.compute_offsets([])
    src.adapt_mag(a0)
    print("Done with offsets")
    assert len(a0) == 10

    photz.fit_onesource(src)
    print("Done with fit")
    assert src.zmin[0] == 0.40
    print("src.dmmin[0]: ", src.dmmin[0])
    assert src.dmmin[0] == pytest.approx(1.0000e10, 1.0e4)
    assert src.imasmin[0] == 2

    photz.uncertainties_onesource(src)
    print("Done with uncertainties")
    max_position = np.argmax(src.pdfmap[11].vPDF)
    assert max_position == 2
    assert len(src.pdfmap[11].xaxis) == 6
    assert np.isclose(src.zs, 0.4, atol=1e-02)
    photz.physpara_onesource(src)
    print("Done with physical parameters")
    assert src.consiz == 0.4

    print(src.results)
    print("synth mag rescaled to observations: ", src.magm)
    print("src.chimin[0]: ", src.chimin[0])
    print("synth mag from template libs: ", photz.fullLib[src.indmin[0]].mag)
    print("mag pseudo-observéés: ", maglib_bc03)

    # Low chi2 since no noise and predicted mag in input
    assert src.chimin[0] < 1.0e-2
    # From the .phys
    # mass 0.0168893 ->log10(0.0168893)+10->8.227
    # sfr  6.50848e-11 -> log10(6.50848e-11)+10 -> -0.1865
    # ssfr -8.4135
    assert np.isclose(src.results["MASS_BEST"], 8.227, atol=1e-02)
    assert np.isclose(src.results["SFR_BEST"], -0.1865, atol=1e-02)
    assert np.isclose(src.results["SSFR_BEST"], -8.4135, atol=1e-02)
    assert np.isclose(src.results["EBV_BEST"], 0.5, atol=1e-02)
    assert np.isclose(src.results["AGE_BEST"], 8.55654, atol=1e-02)
    assert np.isclose(src.results["EXTLAW_BEST"], 1, atol=1e-02)
    assert src.massmed[1] < 8.227 < src.massmed[2]
    assert src.agemed[1] < 8.55654 < src.agemed[2]
    assert src.SFRmed[1] < -0.1865 < src.SFRmed[2]
    # compare the predicted magnitudes
    assert np.allclose(src.magm, maglib_bc03, atol=1e-3)


def test_fit_fir():
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    print(test_dir)

    # Read the config file.
    config_file = os.path.expandvars("$LEPHAREDIR/examples/COSMOS.para")
    config = lp.read_config(config_file)
    print(config_file)
    lsst_files = [
        "lsst/total_u.pb",
        "lsst/total_g.pb",
        "lsst/total_r.pb",
        "lsst/total_i.pb",
        "lsst/total_z.pb",
        "lsst/total_y.pb",
    ]
    spitzer_files = ["spitzer/irac_3.pb", "spitzer/irac_4.pb"]
    herschel_files = ["herschel/PACS_100.pb", "herschel/PACS_160.pb"]
    fltstr = ",".join(lsst_files + spitzer_files + herschel_files)
    config.update(
        {
            "VERBOSE": "NO",
            "FILTER_LIST": fltstr,
            "FILTER_FILE": "filters_fir",
            "FILTER_CALIB": "0,0,0,0,0,0,1,1,1,1",
            "ZPHOTLIB": "BC03_IR",
            "AUTO_ADAPT": "NO",
            "Z_STEP": "0.2,0,1",
            "ZFIX": "YES",
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
            "FIR_LIB": "DALE_FIR",
            "FIR_LMIN": "4.0",
            "FIR_CONT": "-1",
            "FIR_SCALE": "-1",
            "FIR_FREESCALE": "YES",
            "FIR_SUBSTELLAR": "YES",
        }
    )
    keymap = lp.all_types_to_keymap(config)

    # Filters
    filterlib = lp.Filter(config_keymap=keymap)
    filterlib.run()

    ##
    # Download missing data
    lp.data_retrieval.get_auxiliary_data(
        keymap=lp.default_cosmos_config.copy(),
        additional_files=[
            "sed/GAL/BC03_CHAB/bc2003_lr_m62_chab_tau1_dust00.ised_ASCII",
            "sed/GAL/BC03_CHAB/bc2003_lr_m62_chab_tau15_dust00.ised_ASCII",
            "sed/GAL/DALE/dale_1.sed",
            "sed/GAL/DALE/dale_34.sed",
        ],
    )

    ##
    # List of templates and ages

    # Generate a test BC03 template list
    pathlist_fir = os.path.expandvars("$LEPHAREWORK/testFIR_MOD.list")
    print("FIR file ", pathlist_fir)
    lines = [
        "DALE/dale_1.sed   LW          14.3404\n",
        "DALE/dale_34.sed   LW         9.8393\n",
    ]
    with open(os.path.expandvars(pathlist_fir), "w") as file:
        file.writelines(lines)
    # Generate a test BC03 template list
    pathlist_bc03 = os.path.expandvars("$LEPHAREWORK/testBC03_MOD.list")
    print("BC03 file ", pathlist_bc03)
    lines = [
        "BC03_CHAB/bc2003_lr_m62_chab_tau1_dust00.ised_ASCII  BC03\n",
        "BC03_CHAB/bc2003_lr_m62_chab_tau15_dust00.ised_ASCII  BC03\n",
    ]
    with open(os.path.expandvars(pathlist_bc03), "w") as file:
        file.writelines(lines)
    # Generate a test BC03 age list
    pathage_bc03 = os.path.expandvars("$LEPHAREWORK/agesBC03.txt")
    print("BC03 ages ", pathage_bc03)
    lines = ["0.360203008\n", "2\n", "7\n"]
    with open(os.path.expandvars(pathage_bc03), "w") as file:
        file.writelines(lines)

    sedlib = lp.Sedtolib(config_keymap=keymap)
    sedlib.run(
        typ="GAL",
        gal_lib="LIB_BC03",
        gal_sed=pathlist_bc03,
        sel_age=pathage_bc03,
    )

    maglib = lp.MagGal(config_keymap=keymap)
    maglib.run(
        typ="GAL",
        lib_ascii="YES",
        gal_lib_in="LIB_BC03",
        gal_lib_out="BC03_IR",
        eb_v="0.,0.25,0.5",
        mod_extinc="0,12",
        extinc_law="SB_calzetti.dat",
        em_lines="PHYS",
        em_dispersion="1.",
    )

    ##
    # Dale

    # Build the Dale library
    sedlib = lp.Sedtolib(config_keymap=keymap)
    sedlib.run(
        typ="GAL",
        gal_lib="LIB_FIR",
        gal_sed=pathlist_fir,
    )

    maglib = lp.MagGal(config_keymap=keymap)
    maglib.run(
        typ="GAL",
        lib_ascii="YES",
        gal_lib_in="LIB_FIR",
        gal_lib_out="DALE_FIR",
        mod_extinc="0,0",
    )

    print("Done creating libraries")
    assert os.path.exists(os.path.expandvars("$LEPHAREWORK/lib_mag/BC03_IR.doc"))
    assert os.path.exists(os.path.expandvars("$LEPHAREWORK/lib_mag/DALE_FIR.doc"))

    # Instantiate the photoz object
    photz = lp.PhotoZ(keymap)
    print("Done instantiate photoz")
    assert len(photz.gridz) == 6
    assert photz.gridz[0] == 0
    assert photz.gridz[5] == 1

    # Instantiate a source (Id, gridz)
    src = lp.onesource(101, photz.gridz)
    # read the source, change Id, attribute mag/err, ...
    cont = 0
    maglib_bc03 = [51.9275, 51.2836, 50.39, 49.9945, 49.7161, 49.4793, 49.0574, 49.5759, 55.004, 56.0812]
    # Substract 25 to get a mass of 10^10 Msol
    maglib_bc03 = [x - 25 for x in maglib_bc03]
    # LDUST predcited for this src would be 9.85586393 from previous test.
    # So, include the dale template Dale 34 corresponding to this LIR
    # Predicted mag in irac 3, 4, PACS 100, 160: 24.6431 21.9891 17.0063 16.2928
    maglib_dale = [1000, 1000, 1000, 1000, 1000, 1000, 24.6431, 21.9891, 17.0063, 16.2928]
    # Convert magnitudes to flux
    flux_bc03 = [10 ** (-0.4 * (mag + 48.6)) for mag in maglib_bc03]
    flux_dale = [10 ** (-0.4 * (mag + 48.6)) for mag in maglib_dale]
    # Sum the flux values stellar+dust
    flux_src = [flux_bc03[i] + flux_dale[i] for i in range(len(flux_bc03))]
    # Associate uncertainties
    eflux_src = [0.02 * flux_src[i] / 1.048 for i in range(len(flux_src))]
    src.readsource(
        "40",
        flux_src,
        eflux_src,
        cont,
        0.4,
        "ebv0.5_DM41.6845age3.60203e8",
    )

    photz.prep_data(src)

    photz.fit_onesource(src)
    print("Done with fit")
    assert src.zmin[0] == 0.40
    assert src.dmmin[0] == pytest.approx(1.0000e10, 1.0e4)
    assert src.imasmin[0] == 2

    photz.uncertainties_onesource(src)
    print("Done with uncertainties")
    max_position = np.argmax(src.pdfmap[11].vPDF)
    assert max_position == 2
    assert len(src.pdfmap[11].xaxis) == 6
    assert np.isclose(src.zs, 0.4, atol=1e-02)

    photz.physpara_onesource(src)
    print("Done with physical parameters")
    print(src.results)
    assert src.consiz == 0.4
    assert src.zminIR == pytest.approx(0.4, 1.0e4)
    # compare the predicted magnitudes from the BC03 fit and the input ones
    assert np.allclose(src.magm, maglib_bc03, atol=1e-3)
    # compare the original FIR flux and the one obtained after substracted the BC03 flux
    # only for the filters used in FIR
    assert np.allclose(src.abIR[6:10], flux_dale[6:10], atol=1e-30)
    # Rescaling of the FIR library (none)
    assert src.dmminIR == pytest.approx(1.0000, 1.0e4)
    # Low chi2 since no noise and predicted mag in input
    assert src.chiminIR < 1.0e-4
    assert src.chimin[0] < 1.0e-2
    # From the .phys
    # mass 0.0168893 ->log10(0.0168893)+10->8.227
    # sfr  6.50848e-11 -> log10(6.50848e-11)+10 -> -0.1865
    # ssfr -8.4135
    assert np.isclose(src.results["MASS_BEST"], 8.227, atol=1e-02)
    assert np.isclose(src.results["SFR_BEST"], -0.1865, atol=1e-02)
    assert np.isclose(src.results["SSFR_BEST"], -8.4135, atol=1e-02)
    assert np.isclose(src.results["EBV_BEST"], 0.5, atol=1e-02)
    assert np.isclose(src.results["AGE_BEST"], 8.55654, atol=1e-02)
    assert np.isclose(src.results["EXTLAW_BEST"], 1, atol=1e-02)
    # The LTIR from .phys is slightly different from the one announced by Dale
    assert np.isclose(src.results["LUM_TIR_BEST"], 9.74753, atol=1e-02)
    assert np.isclose(src.results["EXTLAW_BEST"], 1, atol=1e-02)
    assert src.massmed[1] < 8.227 < src.massmed[2]
    assert src.agemed[1] < 8.55654 < src.agemed[2]
    assert src.SFRmed[1] < -0.1865 < src.SFRmed[2]
    # compare the predicted magnitudes
    assert np.isclose(src.LIRmed[0], 9.74753, atol=1e-01)
    assert src.LIRmed[1] < 9.74753 < src.LIRmed[2]

    minl = 3000.0
    maxl = 200000.0
    gal1 = photz.besttemplate_onesource(src, 0, minl, maxl)
    gal2 = photz.besttemplate_onesource(src, 1, minl, maxl)
    fir = photz.besttemplate_onesource(src, 2, minl, maxl)
    qso = photz.besttemplate_onesource(src, 3, minl, maxl)
    star = photz.besttemplate_onesource(src, 4, minl, maxl)
    print("Done with templates generation")
    assert all(2999.0 < x < 200001.0 for x in gal1[0])
    assert all(10.0 < x < 1000.1 for x in gal1[1])
    assert len(gal2[0]) == 0
    assert all(2999.0 < x < 200001.0 for x in fir[0])
    assert all(10.0 < x < 1000.1 for x in fir[1])
    assert all(2999.0 < x < 200001.0 for x in qso[0])
    assert all(10.0 < x < 1000.1 for x in qso[1])
    assert all(2999.0 < x < 200001.0 for x in star[0])
    assert all(10.0 < x < 1000.1 for x in star[1])
    print("Done with the spectra")
