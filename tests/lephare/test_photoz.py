import os

import lephare as lp
import numpy as np
from astropy.table import Table


def test_photoz(test_data_dir: str):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    lp.prepare(config)
    photz = lp._lephare.PhotoZ(config)
    # Check the redshift zero case works
    for i, sed in enumerate(photz.fullLib):
        # Check some spectra data is available
        # assert np.array(sed.get_data_vector(0.0, 1000000.0, True, 0.0)).shape != (2, 0)
        # remove for now to use other tests
        # this should eventually test the python code for getting the spectra
        if not np.isclose(photz.zLib[i], 0.0):
            # Check redshifted sources have same model as previous sed
            assert sed.nummod == photz.fullLib[i - 1].nummod


def test_reddening(test_data_dir: str):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    config["Z_STEP"] = "1.,0.,2."  # Fake star SED gets redshifted out of B band at low z
    # keymap=lp.all_types_to_keymap(config)
    config["APPLY_MW_EXTINCTION"] = "GALAMETZ"
    config["EXT_MW_CURVE"] = "LMC_Fitzpatrick.dat"
    config["MW_REFERENCE_MODEL"] = "sed/STAR/PICKLES/b5i.sed"
    mw_ebv_test_file = os.path.join(test_data_dir, "examples/mw_ebv.dat")
    config["MW_EBV_FILE"] = mw_ebv_test_file
    # the traditional file must be written later
    config["CAT_IN"] = os.path.join(test_data_dir, "examples/COSMOS_first100specz_reduced.in")

    lp.prepare(config)

    # The reddening calculator
    albd_lib = lp.compute_model_reddening(config)
    assert albd_lib.shape == (19, 2)
    # test impact of ebv on a source fit
    # Read the test input catalogue
    input_file = os.path.join(test_data_dir, "examples/COSMOS_first100specz.fits")
    input = Table.read(input_file)
    ebv_test = np.linspace(0.0, 0.3, len(input))
    out = Table()
    out["id"] = input[input.colnames[0]]
    out["ebv"] = ebv_test
    out.write(mw_ebv_test_file, format="ascii.no_header", overwrite=True)

    test_string = "te s"  # Test with spaces
    input["string_input"][0] = test_string
    # Make a reduced column set for the minimal test
    reduced_cols = []
    for c in input.colnames:
        if not c.startswith("f"):
            reduced_cols.append(c)
        elif "IB527" in c:
            reduced_cols.append(c)
        elif "IB679" in c:
            reduced_cols.append(c)

    output, photozlist = lp.process(
        config, input[reduced_cols], write_outputs=False, mw_ebv=[0.1] * len(input)
    )
    assert np.isclose(np.sum(albd_lib), 105.84491378690295)

    # Check it can read the ebv values from the file and apply them to the sources
    input[reduced_cols].write(
        os.path.join(test_data_dir, "examples/COSMOS_first100specz_reduced.in"),
        format="ascii.no_header",
        overwrite=True,
    )
    photz = lp.PhotoZ(lp.all_types_to_keymap(config))
    sources = photz.read_photoz_sources()
    photz.read_mw_ebv(sources)
    # This is not working at the moment as the sources are not being read in
    # # with the ebv values, but this should be tested eventually
    for n, s in enumerate(sources):
        assert np.isclose(s.mw_ebv, ebv_test[n])

    # This is getting the ebv values from the file.
    output, _ = lp.process(config, input[reduced_cols], write_outputs=False)
    assert np.isclose(np.sum(output["Z_BEST"]), 84.0)
    # Now test with SMC Prevot curve
    config["EXT_MW_CURVE"] = "MW_seaton.dat"
    albd_lib = lp.compute_model_reddening(config)
    # Check AUTO_ADAPT can run with ebv
    config["AUTO_ADAPT"] = "YES"

    output, photozlist = lp.process(
        config, input[reduced_cols], write_outputs=False, mw_ebv=[0.1] * len(input)
    )

    # Test the band pass correction
    bpc = lp.compute_band_pass_correction(config)
    assert np.isclose(np.sum(bpc), 12.141620831207208)


def test_build_output_tables_fills_missing_values_with_nan(test_data_dir, tmp_path, monkeypatch):
    """Output keys that index past the end of a source attribute become NaN.

    alloutputkeys.txt maps each output name onto an attribute of the source
    object, sometimes with an index. If that index (or dict key) is absent the
    column is filled with NaN rather than raising.
    """
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")

    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    lp.prepare(config)
    input_table = Table.read(os.path.join(test_data_dir, "examples/COSMOS_first100specz.fits"))
    reduced_cols = [c for c in input_table.colnames if not c.startswith("f") or "IB527" in c or "IB679" in c]
    _, srclist = lp.process(config, input_table[reduced_cols][:5], write_outputs=False)

    # A key table pointing at one valid value, one out-of-range index and one
    # absent dict key.
    key_table = tmp_path / "alloutputkeys.txt"
    key_table.write_text(
        "Z_BEST\tfloat\tzgmin[0]\n"
        "OUT_OF_RANGE\tfloat\tzgmin[9]\n"
        'MISSING_KEY\tfloat\tresults["NOT_A_REAL_KEY"]\n'
    )
    para_out = tmp_path / "custom_output.para"
    para_out.write_text("Z_BEST\nOUT_OF_RANGE\nMISSING_KEY\n")

    monkeypatch.setenv("LEPHAREDIR", str(tmp_path))
    photz = lp.PhotoZ(lp.all_types_to_keymap(config))
    table = photz.build_output_tables(srclist, para_out=str(para_out))

    assert set(table.colnames) >= {"Z_BEST", "OUT_OF_RANGE", "MISSING_KEY"}
    # The valid index produced real redshifts
    assert np.all(np.isfinite(table["Z_BEST"]))
    # Both unavailable values fell back to NaN
    assert np.all(np.isnan(table["OUT_OF_RANGE"]))
    assert np.all(np.isnan(table["MISSING_KEY"]))
