import os

import lephare as lp


def test_photoz(test_data_dir: str):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    lp.prepare(config)
    photz = lp._lephare.PhotoZ(config)
    # Check the redshift zero case works
    print(photz.fullLib[0].data().shape)
    assert photz.fullLib[0].data().shape == (2, 4860)
    # Check the non redshift zero case works
    assert photz.fullLib[1].data().shape == (2, 2492)
