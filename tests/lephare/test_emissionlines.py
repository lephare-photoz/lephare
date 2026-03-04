import os

import lephare as lp
import numpy as np
import pytest


def test_emlines():
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")

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
    fltstr = ",".join(lsst_files)
    config.update(
        {
            "VERBOSE": "NO",
            "FILTER_LIST": fltstr,
            "FILTER_FILE": "filters_lsst",
            "FILTER_CALIB": "0,0,0,0,0,0",
            "ZPHOTLIB": "BC03_EL",
            "AUTO_ADAPT": "NO",
            "Z_STEP": "0.2,0,1",
            "ZFIX": "YES",
            "EB_V": "0.,0.1,0.2,0.3",
            "MOD_EXTINC": "0,100",
            "ADD_EMLINES": "0,100",
            "EM_LINES": "PHYS",
            "EM_DISPERSION": "0.5,1.,2",
            "ERR_SCALE": " 0.0",
            "ERR_FACTOR": " 1.",
            "Z_INTERP": "NO",
            "MAG_ABS": "-24,-5",
            "MAG_REF": "2",
        }
    )
    keymap = lp.all_types_to_keymap(config)

    # Filters
    filterlib = lp.Filter(config_keymap=keymap)
    filterlib.run()

    # Generate a test BC03 template list
    pathlist_bc03 = os.path.expandvars("$LEPHAREWORK/testEL_MOD.list")
    print("BC03 file ", pathlist_bc03)
    lines = [
        "BC03_CHAB/bc2003_lr_m62_chab_tau15_dust00.ised_ASCII  BC03\n",
    ]
    with open(os.path.expandvars(pathlist_bc03), "w") as file:
        file.writelines(lines)
    # Generate a test BC03 age list
    pathage_bc03 = os.path.expandvars("$LEPHAREWORK/agesBC03.txt")
    print("BC03 ages ", pathage_bc03)
    lines = ["2\n"]
    with open(os.path.expandvars(pathage_bc03), "w") as file:
        file.writelines(lines)

    # Download large files if needed
    lp.data_retrieval.get_auxiliary_data(
        keymap=lp.default_cosmos_config.copy(),
        additional_files=[
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
        gal_lib_out="BC03_EL",
        eb_v="0.,0.5",
        mod_extinc="0,12",
        extinc_law="SB_calzetti.dat",
    )

    print("Done creating libraries")
    assert os.path.exists(os.path.expandvars("$LEPHAREWORK/lib_mag/BC03_EL.doc"))

    # Instantiate the photoz object
    photz = lp.PhotoZ(keymap)
    print("Done instantiate photoz")
    assert len(photz.gridz) == 6
    assert photz.gridz[0] == 0
    assert photz.gridz[5] == 1

    # Instantiate a source (Id, gridz)
    src = lp.onesource(101, photz.gridz)
    # read the source, change Id, attribute mag/err, ...
    maglib_bc03 = [48.2199, 48.0719, 47.5242, 47.4207, 47.3486, 47.297]
    # flux from the .dat
    em_bc03 = [4.04849e-42, 1.52993e-40, 1.91321e-40, 4.96433e-40, 4.48886e-40, 3.64144e-40]
    # Substract 25 to get a mass of 10^10 Msol
    flux_bc03 = [10 ** (-0.4 * (x - 25 + 48.6)) for x in maglib_bc03]
    em_bc03 = [x * 1.0e10 for x in em_bc03]
    # Sum continuum and emmision lines
    sum_em = [x + y for x, y in zip(flux_bc03, em_bc03)]
    esum_em = [x * 0.01 for x in sum_em]
    src.readsource("40", sum_em, esum_em, 0, 0.4, "bid")
    print("Done creating source")

    photz.prep_data(src)
    a0 = photz.compute_offsets([])
    src.adapt_mag(a0)
    print("Done with offsets")
    assert len(a0) == 6

    photz.fit_onesource(src)
    print("Done with fit")
    assert src.zmin[0] == 0.40
    print("src.dmmin[0]: ", src.dmmin[0])
    assert src.dmmin[0] == pytest.approx(1.0000e10, 1.0e4)
    assert src.imasmin[0] == 1

    photz.uncertainties_onesource(src)
    print("Done with uncertainties")
    max_position = np.argmax(src.pdfmap[11].vPDF)
    assert max_position == 2
    assert len(src.pdfmap[11].xaxis) == 6
    assert src.zs == pytest.approx(0.4, abs=1e-02)
    photz.physpara_onesource(src)
    print("Done with physical parameters")
    assert src.consiz == 0.4

    # Low chi2 since no noise and predicted mag in input
    assert src.chimin[0] < 1.0e-2
    # From the .phys
    # mass 0.077602 ->log10(0.077602)+10->8.8898
    # sfr  5.83449e-11 -> log10(5.83449e-11)+10 -> -0.234
    assert src.results["MASS_BEST"] == pytest.approx(8.8898, abs=1e-02)
    assert src.results["SFR_BEST"] == pytest.approx(-0.234, abs=1e-02)
    assert src.results["EBV_BEST"] == pytest.approx(0.0, abs=1e-02)

    # Test the value of the reconstructed emission line
    print(src.results_emission_lines)
    # qi: 722.312 * mass
    # include the distance 10x10pc2 / dl^2 in pc
    # Hbeta expected = 722.312 * 1.e10 * 100. / 4718018.41 * 1.e-12   * 4.780e-13;
    hb_expected = 7.3180116e-17
    fhb = src.results_emission_lines["EM_FLUX_HB"]
    fha = src.results_emission_lines["EM_FLUX_HA"]
    ha_hb_ratio = fha / fhb
    assert ha_hb_ratio == pytest.approx(2.87, abs=1e-02)
    assert fhb == pytest.approx(hb_expected, abs=1e-18)
