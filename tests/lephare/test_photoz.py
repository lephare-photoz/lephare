import os

import lephare as lp
import numpy as np
from astropy.table import Table

out_z_true = np.array(
    [
        3.57892305,
        3.18375599,
        0.7,
        3.85,
        2.09518357,
        0.80034386,
        1.48439925,
        2.2939652,
        1.31775377,
        0.86908708,
        3.15323574,
        2.98449005,
        0.98326861,
        0.55153602,
        1.18570296,
        3.53229024,
        2.70674634,
        2.83820855,
        2.94946781,
        3.15146252,
        1.5013684,
        2.84600841,
        1.70175661,
        1.1829986,
        1.32712949,
        1.70385503,
        1.1149645,
        0.98851633,
        3.24705196,
        2.9,
        3.17973438,
        0.81103435,
        1.2,
        2.71450225,
        3.4,
        1.3242741,
        0.98035073,
        3.18077313,
        1.05401207,
        1.49872757,
        1.06706739,
        0.98592478,
        3.1,
        0.85546416,
        1.21506092,
        0.84430386,
        3.24247515,
        1.59507083,
        1.13636467,
        1.18293656,
        1.48957004,
        1.21998967,
        3.24113914,
        1.1697655,
        1.50354007,
        3.35078231,
        2.09561018,
        1.06851229,
        3.4,
        1.16582576,
        3.14287899,
        1.073411,
        3.35438822,
        0.85920482,
        1.32161593,
        2.67604136,
        1.19677899,
        3.2467587,
        1.21231392,
        0.98022891,
        1.49795109,
        3.56991808,
        1.58425944,
        1.31882606,
        1.3859304,
        1.19218187,
        2.93987793,
        1.12442176,
        3.2439558,
        1.51679845,
        2.67859993,
        1.70380004,
        1.03814283,
        1.60779959,
        2.84701408,
        0.82986633,
        1.07926986,
        1.0,
        1.18625949,
        1.33592357,
        0.97763736,
        1.35481174,
        1.03893961,
        1.59652143,
        0.8,
        0.83029193,
        1.47577102,
        3.18796741,
        0.85742579,
        2.83919201,
    ]
)


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

    # Check it can run without an x axis
    assert len(lp.multiply_on_grids(xtest, ytest, xtest, ytest)[0]) == len(xtest)

    # The reddening calculator
    albd_lib = lp.compute_model_reddening(config)
    assert albd_lib.shape == (307, 2)
    # test impact of ebv on a source fit
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
    output, photozlist = lp.process(
        config, input[reduced_cols], write_outputs=False, reddening=albd_lib, ebv=[0.1] * len(input)
    )
    np.testing.assert_allclose(np.array(output["Z_BEST"]), out_z_true, rtol=1e-7, atol=0)
