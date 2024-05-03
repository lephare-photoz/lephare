import os

import lephare as lp
import numpy as np
from astropy.table import Table


def test_process(test_data_dir: str):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    # Read the config file.
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    # Run preparation tasks.
    lp.prepare(config)
    # Test on standard table format
    input = Table.read("../data/examples/COSMOS_first100specz.fits")
    output, pdfs = lp.process(input, config)
    # Check one fo the outputs
    assert np.isclose(output["ZBEST"][0], 4.6945872848190575)
