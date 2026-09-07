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
    assert len(gal_mag.colnames) == 12  # Contains additional EM_DISPERSION column


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


def test_prepare_skips_object_types(test_data_dir, set_env_vars):
    """Passing None for an object type skips its sedtolib/mag_gal stage entirely."""
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    work = os.environ["LEPHAREWORK"]

    # read_config yields keyword objects, so take .value for the file names
    keymap = all_types_to_keymap(config)

    # Only run the galaxy stage, leaving gal_config at its default
    star_yaml = os.path.join(work, "lib_bin", f"{keymap['STAR_LIB'].value}_star_config.yaml")
    if os.path.exists(star_yaml):
        os.remove(star_yaml)

    lp.prepare(config, star_config=None, qso_config=None)

    # The galaxy stage ran and wrote its config alongside the library
    assert os.path.exists(os.path.join(work, "lib_bin", f"{keymap['GAL_LIB'].value}_gal_config.yaml"))
    # The skipped types wrote nothing
    assert not os.path.exists(star_yaml)
    # The filter stage always runs, regardless of which types are selected
    assert os.path.exists(os.path.join(work, "filt", f"{keymap['FILTER_FILE'].value}.dat"))


def test_prepare_overrides_config_per_type(test_data_dir, set_env_vars):
    """A per-type config overrides the base config for that type only."""
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    work = os.environ["LEPHAREWORK"]

    lp.prepare(config, star_config=None, gal_config={"GAL_LIB_OUT": "OVERRIDDEN_GAL"}, qso_config=None)

    # The override is reflected in the name of the written magnitude library
    assert os.path.exists(os.path.join(work, "lib_mag", "OVERRIDDEN_GAL_gal_config.yaml"))
    written = lp.read_yaml_config(os.path.join(work, "lib_mag", "OVERRIDDEN_GAL_gal_config.yaml"))
    assert written["GAL_LIB_OUT"].value == "OVERRIDDEN_GAL"


def test_overwrite_config():
    """overwrite_config merges the second config over the first."""
    base = all_types_to_keymap({"A": "1", "B": "2"})
    override = all_types_to_keymap({"B": "changed", "C": "new"})

    merged = lp.overwrite_config(base, override)
    assert merged["A"].value == "1"
    assert merged["B"].value == "changed"
    assert merged["C"].value == "new"


def test_overwrite_config_with_none():
    """A None override leaves the base config untouched."""
    base = all_types_to_keymap({"A": "1"})
    assert lp.overwrite_config(base, None) is base


def test_yaml_config_round_trip(tmp_path):
    """A keymap survives a write/read cycle through yaml."""
    keymap = all_types_to_keymap({"Z_STEP": "0.04,0,6", "FILTER_FILE": "filter_cosmos"})
    yaml_path = tmp_path / "config.yaml"

    lp.write_yaml_config(keymap, str(yaml_path))
    assert yaml_path.exists()
    # The file is stamped with the version that wrote it
    assert lp.__version__ in yaml_path.read_text().splitlines()[0]

    restored = lp.read_yaml_config(str(yaml_path))
    assert set(restored) == set(keymap)
    for key in keymap:
        assert type(restored[key]) == lp.keyword
        assert restored[key].value == keymap[key].value


def test_write_yaml_config_accepts_string_dict(tmp_path):
    """write_yaml_config converts a plain string dict before writing."""
    yaml_path = tmp_path / "strings.yaml"
    lp.write_yaml_config({"Z_STEP": "0.04,0,6"}, str(yaml_path))

    restored = lp.read_yaml_config(str(yaml_path))
    assert restored["Z_STEP"].value == "0.04,0,6"


def test_write_para_config_round_trip(tmp_path):
    """A .para file written from a keymap can be read back by a Runner."""
    para_path = tmp_path / "config.para"
    lp.write_para_config({"CAT_IN": "somefile.in"}, str(para_path))

    runner = lp.Runner(config_keys={"CAT_IN": "help"}, config_file=str(para_path))
    assert runner.keymap["CAT_IN"].value == "somefile.in"
