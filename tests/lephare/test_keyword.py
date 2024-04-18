import os
import tempfile

import numpy as np
from lephare import keyword, read_command, read_config

TESTDIR = os.path.abspath(os.path.dirname(__file__))
LEPHAREDIR = os.path.join(TESTDIR, "..")


def test_keyword_basic():
    # test that default is taken when value is empty
    k = keyword("k1", "")
    assert k.split_string("DEFAULT", 1)[0] == "DEFAULT"
    # check basic assignment
    k = keyword("k2", "1")
    assert k.name == "k2"
    assert k.value == "1"
    # test mismatch in size
    assert np.array_equal(k.split_string("0", 2), ["1", "1"])
    k = keyword("k3", "1,1")
    assert np.array_equal(k.split_string("0", 3), ["0", "0", "0"])


def test_keyword_array():
    k = keyword("k", "1,1")
    assert np.array_equal(k.split_string("0", 2), ["1", "1"])
    assert np.array_equal(k.split_int("0", 2), [1, 1])
    assert np.array_equal(k.split_long("0", 2), [1, 1])
    assert np.array_equal(k.split_bool("0", 2), [True, True])
    assert np.array_equal(k.split_double("0", 2), [1.0, 1.0])


def test_keyword_bool():
    k = keyword("k", "1")
    assert k.split_bool("0", 1)[0]
    k = keyword("k", "YES")
    assert k.split_bool("0", 1)[0]
    k = keyword("k", "yes")
    assert k.split_bool("0", 1)[0]
    k = keyword("k", "True")
    assert k.split_bool("0", 1)[0]
    k = keyword("k", "true")
    assert k.split_bool("0", 1)[0]
    k = keyword("k", "ANYTHING ELSE")
    assert not k.split_bool("0", 1)[0]


def test_keyword_expand_path():
    lepharedir = os.environ["LEPHAREDIR"]
    home = os.environ["HOME"]
    k = keyword("k", "$LEPHAREDIR/unittest")
    assert k.value == os.path.join(lepharedir, "unittest")
    k = keyword("k", "$HOME/unittest")
    assert k.value == os.path.join(home, "unittest")


def test_keyword_read_config(test_data_dir):
    config_file = os.path.join(test_data_dir, "examples/COSMOS.para")
    km = read_config(config_file)
    for k in km:
        assert k[0] != "#"
        assert k[0] != " "
        assert k[0] != "\t"


def test_keyword_read_command():
    km = read_command(["exec", "-c", "config", "--arg1", "val1", "-arg2", "val2"])
    for k in km:
        if k == "c":
            assert km[k].value == "config"
        if k == "arg1":
            assert km[k].value == "val1"
        if k == "arg2":
            assert km[k].value == "val2"


def test_keyword_misformed_kwd():
    # normal
    k = keyword("k", "1,2,3")
    assert np.array_equal(k.split_string("0", -1), ["1", "2", "3"])
    # space
    k = keyword("k", "1,2, 3")
    assert np.array_equal(k.split_string("0", -1), ["1", "2", "3"])
    # tab
    k = keyword("k", "1,2,\t3")
    assert np.array_equal(k.split_string("0", -1), ["1", "2", "3"])

    with tempfile.TemporaryDirectory() as dir_name:
        filename = os.path.join(dir_name, "tmp")
        with open(filename, "w") as f:
            f.writelines(
                [
                    "FILTER_CALIB    0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,\t0",
                    "\nFILTER_FILE     filter_cosmos",
                ]
            )
        km = read_config(filename)
        assert np.array_equal(km["FILTER_CALIB"].value, "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0")
        assert np.array_equal(km["FILTER_FILE"].value, "filter_cosmos")
