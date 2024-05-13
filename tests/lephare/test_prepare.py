import os
import shutil

import lephare as lp


def test_prepare(test_data_dir):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    lp.prepare(config)


def test_load_sed_list(test_data_dir):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    # Move one of the example sed folders
    _ = shutil.copytree(os.path.join(test_dir, "../data/sed/QSO"), os.path.join(test_dir, "../tmp/tmp"))
    lp.load_sed_list(os.path.join(test_dir, "../tmp/tmp/ONE_SED.list"), "QSO")
    # Check the list is there
    assert os.path.exists(os.path.join(test_dir, "../data/sed/QSO/ONE_SED/ONE_SED.list"))
    # Check the sed is there
    assert os.path.exists(os.path.join(test_dir, "../data/sed/QSO/ONE_SED/o5v.sed.ext"))
