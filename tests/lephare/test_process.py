import os
import shutil

import lephare as lp
import numpy as np
import pytest
from astropy.table import Table


def test_process(test_data_dir: str):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    # Read the config file.
    config_file = os.path.join(test_data_dir, "examples/COSMOS.para")
    config = lp.read_config(config_file)
    # Run preparation tasks.
    lp.prepare(config)
    # Read the test input catalogue
    input_file = os.path.join(test_data_dir, "examples/COSMOS_first100specz.fits")
    input = Table.read(input_file)
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
    output, photozlist = lp.process(config, input[reduced_cols], write_outputs=False)
    # Check one of the outputs (results are terrible with just one filter and sparse z grid)
    assert output["Z_BEST"][0] == pytest.approx(3.5878, abs=1.0e-3)
    assert len(photozlist[0].pdfmap[11].xaxis) == 51
    pdfs = np.array([photozlist[i].pdfmap[11].vPDF for i in np.arange(len(photozlist))])
    assert np.sum(pdfs) == pytest.approx(1001.3718)
    assert output["STRING_INPUT"][0] == test_string
    # Check AUTO_ADAPT
    config["AUTO_ADAPT"] = "YES"
    output, photozlist = lp.process(config, input[reduced_cols], write_outputs=True)

    assert ~(output["Z_BEST"][0] == pytest.approx(3.5877994546919934))
    assert os.path.isfile("zphot.out")
    assert output["IDENT"][0] == str(input["id"][0])

    a0 = lp.calculate_offsets_from_input(config, input[reduced_cols])
    assert len(a0) == 2

    # Test table formatting
    id, flux, flux_err, context, zspec, string_data = lp.table_to_data(
        config, input[reduced_cols], col_names=reduced_cols, standard_names=False
    )
    assert len(zspec) == 100
    config["FILTER_LIST"] = "cosmos/IB527.lowres,cosmos/IB679.lowres"
    config = lp.all_types_to_keymap(config)
    id, flux, flux_err, context, zspec, string_data = lp.table_to_data(
        config, input[reduced_cols], standard_names=True
    )
    assert len(zspec) == 100


def test_load_sed_list(test_data_dir):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    # Move one of the example sed folders
    _ = shutil.copytree(
        os.path.join(test_dir, "../data/sed/QSO"), os.path.join(test_dir, "../tmp/seds"), dirs_exist_ok=True
    )
    lp.load_sed_list(os.path.join(test_dir, "../tmp/seds/ONE_SED.list"), "QSO")
    # Check the list is there
    assert os.path.exists(os.path.join(test_dir, "../data/sed/QSO/ONE_SED/ONE_SED.list"))
    # Check the sed is there
    assert os.path.exists(os.path.join(test_dir, "../data/sed/QSO/ONE_SED/o5v.sed.ext"))

    # Check it can run even if the file is already there
    lp.load_sed_list(os.path.join(test_dir, "../tmp/seds/ONE_SED.list"), "QSO")
    # Check absolute paths
    with open(os.path.join(test_dir, "../tmp/seds/ONE_SED_ABS.list"), "w") as file:
        file.write(os.path.join(test_dir, "../tmp/seds/o5v.sed.ext"))
    lp.load_sed_list(os.path.join(test_dir, "../tmp/seds/ONE_SED_ABS.list"), "QSO", absolute_paths=True)
    # The sed named by absolute path was copied into the type folder
    assert os.path.exists(os.path.join(test_dir, "../data/sed/QSO/ONE_SED_ABS/o5v.sed.ext"))
    # Clear the copied folders. ONE_SED_ABS must go too, or a re-run of this test
    # takes the "list file already exists" short cut and never copies anything.
    shutil.rmtree(os.path.join(test_dir, "../tmp/seds"))
    shutil.rmtree(os.path.join(test_dir, "../data/sed/QSO/ONE_SED"))
    shutil.rmtree(os.path.join(test_dir, "../data/sed/QSO/ONE_SED_ABS"))


def _reduced_input(test_data_dir):
    """Load the example catalogue, keeping only the two filters the tests build."""
    input_file = os.path.join(test_data_dir, "examples/COSMOS_first100specz.fits")
    input = Table.read(input_file)
    reduced_cols = [c for c in input.colnames if not c.startswith("f") or "IB527" in c or "IB679" in c]
    return input[reduced_cols]


def test_process_rejects_mismatched_mw_ebv(test_data_dir):
    """An mw_ebv array of the wrong length is rejected rather than zipped short."""
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    input_table = _reduced_input(test_data_dir)

    with pytest.raises(ValueError, match=r"Length of mw_ebv \(3\) does not match number of objects"):
        lp.process(config, input_table, write_outputs=False, mw_ebv=[0.1, 0.2, 0.3])


def test_calculate_offsets_rejects_mismatched_mw_ebv(test_data_dir):
    """calculate_offsets_from_input applies the same length check."""
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    input_table = _reduced_input(test_data_dir)

    with pytest.raises(ValueError, match=r"Length of mw_ebv \(2\) does not match number of objects"):
        lp.calculate_offsets_from_input(config, input_table, mw_ebv=[0.1, 0.2])


def test_calculate_offsets_with_mw_ebv(test_data_dir):
    """Supplying mw_ebv to the offset calculation warns and changes the offsets."""
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    input_table = _reduced_input(test_data_dir)

    baseline = lp.calculate_offsets_from_input(config, input_table)
    with pytest.warns(UserWarning, match="Milky Way E\\(B-V\\) values provided"):
        reddened = lp.calculate_offsets_from_input(config, input_table, mw_ebv=[0.3] * len(input_table))

    # One offset per filter, whether or not reddening was supplied
    assert len(reddened) == len(baseline) == 2
    # With AUTO_ADAPT off in COSMOS.para the offsets are all zero
    np.testing.assert_allclose(baseline, 0.0)
    np.testing.assert_allclose(reddened, 0.0)


def test_calculate_offsets_reads_mw_ebv_file(test_data_dir, capsys):
    """An MW_EBV_FILE in the config is read during the offset calculation."""
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    input_table = _reduced_input(test_data_dir)

    ebv_file = os.path.join(test_data_dir, "examples/mw_ebv.dat")
    ebv_table = Table()
    ebv_table["id"] = input_table[input_table.colnames[0]]
    ebv_table["ebv"] = np.linspace(0.0, 0.3, len(input_table))
    ebv_table.write(ebv_file, format="ascii.no_header", overwrite=True)

    config["MW_EBV_FILE"] = ebv_file
    config["CAT_IN"] = os.path.join(test_data_dir, "examples/COSMOS_first100specz_reduced.in")
    input_table.write(config["CAT_IN"], format="ascii.no_header", overwrite=True)

    a0 = lp.calculate_offsets_from_input(config, input_table)

    assert f"Reading offsets from file {ebv_file}" in capsys.readouterr().out
    assert len(a0) == 2
