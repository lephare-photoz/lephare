import os

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
    # Make a reduced column set for the minimal test
    reduced_cols = []
    for c in input.colnames:
        if not c.startswith("f"):
            reduced_cols.append(c)
        elif "IB527" in c:
            reduced_cols.append(c)
        elif "IB679" in c:
            reduced_cols.append(c)
    output, pdfs, zgrid = lp.process(config, input[reduced_cols])
    # Check one of the outputs (results are terrible with just one filter and sparse z grid)
    assert np.isclose(output["Z_BEST"][0], 3.58577857627795)
    assert len(zgrid) == 51
    assert np.isclose(np.sum(pdfs), 1001.3400039697148)
