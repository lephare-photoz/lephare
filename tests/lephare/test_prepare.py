import os
import shutil

import lephare as lp


def test_prepare(test_data_dir):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    lp.prepare(config)
    # Check the default config is consistent
    default_config = lp.default_cosmos_config
    assert config["Z_STEP"] != default_config["Z_STEP"]
    assert default_config["FILTER_REP"] == str(os.path.join(test_data_dir, "filt"))
    # Check it made the galaxy binary file
    assert os.path.exists(os.path.join(test_dir, "../tmp/lib_mag/CE_COSMOS.bin"))


def test_load_sed_list(test_data_dir):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    # Move one of the example sed folders
    _ = shutil.copytree(os.path.join(test_dir, "../data/sed/QSO"), os.path.join(test_dir, "../tmp/seds"))
    lp.load_sed_list(os.path.join(test_dir, "../tmp/seds/ONE_SED.list"), "QSO")
    # Check the list is there
    assert os.path.exists(os.path.join(test_dir, "../data/sed/QSO/ONE_SED/ONE_SED.list"))
    # Check the sed is there
    assert os.path.exists(os.path.join(test_dir, "../data/sed/QSO/ONE_SED/o5v.sed.ext"))
    # Clear the copied folder
    shutil.rmtree(os.path.join(test_dir, "../tmp/seds"))
