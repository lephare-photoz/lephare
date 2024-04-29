import os

import lephare as lp


def test_prepare(test_data_dir):
    print(os.environ["LEPHAREDIR"])
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))

    lp.prepare(config)
