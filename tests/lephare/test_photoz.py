import os

import lephare as lp
import numpy as np


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
