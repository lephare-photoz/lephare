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


def test_reddening(test_data_dir: str):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    # keymap=lp.all_types_to_keymap(config)
    lp.prepare(config)
    # Test the reddening computation function
    xtest = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    ytest = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    assert len(lp.multiply_on_grids(xtest, ytest, xtest, ytest, xtest)[0]) == len(xtest)

    # The reddening calculator
    albd_lib = lp.compute_model_reddening(config)
    assert albd_lib.shape == (307, 2)
    # test impact of ebv on a source fit
