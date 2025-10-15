import os

import lephare as lp
from astropy.table import Table
from lephare.prepare import all_types_to_keymap


def test_version():
    """Check we have a version."""
    assert isinstance(lp.__version__, str)


def test_prepare(test_data_dir: str):
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
    gal_mag = Table.read(os.path.join(test_dir, "../tmp/lib_mag/CE_COSMOS.dat"), format="ascii")
    assert len(gal_mag.colnames) == 14  # Contains additional EM_DISPERSION column


def test_all_types_to_keymap():
    for in_dict in [{"key": "map"}, {"key": lp.keyword("key", "map")}]:
        out_dict = all_types_to_keymap(in_dict)
        for i, o in zip(in_dict, out_dict):
            assert i == o
            assert out_dict[o].__class__ == lp.keyword
            assert out_dict[o].name == "key"
            assert out_dict[o].value == "map"


def test_config_formatting():
    """Some simple tests of configs and switching between formats"""
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    # Make a string dict config and keymap
    config = lp.default_cosmos_config.copy()
    keymap = all_types_to_keymap(config)
    # Check they are consistent
    assert config["FILTER_FILE"] == keymap["FILTER_FILE"].value
    # Check the type of the keymap is a lephare.keyword
    assert type(keymap["Z_STEP"]) == lp.keyword
    lp.write_para_config(config, os.path.join(os.environ["LEPHAREWORK"], "test.para"))
    # Check it made the file
    assert os.path.exists(os.path.join(os.environ["LEPHAREWORK"], "test.para"))
    os.remove(os.path.join(os.environ["LEPHAREWORK"], "test.para"))
    # Make the file again from the keymap instead
    lp.write_para_config(keymap, os.path.join(os.environ["LEPHAREWORK"], "test.para"))
    # Check it also made the file from the keymap
    assert os.path.exists(os.path.join(os.environ["LEPHAREWORK"], "test.para"))
    os.remove(os.path.join(os.environ["LEPHAREWORK"], "test.para"))
    # Check it converts correctly back to a string dict
    assert lp.string_dict_to_keymap(config)["FILTER_FILE"].value == "filter_cosmos"
    # Check the reverse conversion
    assert lp.keymap_to_string_dict(keymap)["FILTER_FILE"] == "filter_cosmos"
