import os
import shutil

import lephare as lp
import numpy as np
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
    assert np.isclose(output["Z_BEST"][0], 3.5877994546919934)
    assert len(photozlist[0].pdfmap[11].xaxis) == 51
    pdfs = np.array([photozlist[i].pdfmap[11].vPDF for i in np.arange(len(photozlist))])
    assert np.isclose(np.sum(pdfs), 1001.2774052829275)
    assert output["STRING_INPUT"][0] == test_string
    # Check AUTO_ADAPT
    config["AUTO_ADAPT"] = "YES"
    output, photozlist = lp.process(config, input[reduced_cols], write_outputs=True)

    assert ~np.isclose(output["Z_BEST"][0], 3.5877994546919934)
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
    # Clear the copied folders
    shutil.rmtree(os.path.join(test_dir, "../tmp/seds"))
    shutil.rmtree(os.path.join(test_dir, "../data/sed/QSO/ONE_SED"))
