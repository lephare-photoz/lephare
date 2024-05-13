import os

import pytest
from lephare import GalMag, MagSvc, QSOMag, StarMag


def test_magsvc(test_data_dir):
    # Test MagSvc class
    config_file_path = os.path.join(test_data_dir, "examples/COSMOS.para")
    keymap = MagSvc.from_config(objtype="GAL", config_file=config_file_path)
    assert type(keymap) == GalMag

    keymap = MagSvc.from_config(objtype="Star", config_file=config_file_path)
    assert type(keymap) == StarMag

    keymap = MagSvc.from_config(objtype="qso", config_file=config_file_path)
    assert type(keymap) == QSOMag


def test_magsvc_unknown_objtype(test_data_dir):
    # Test MagSvc class with unknown objtype
    config_file_path = os.path.join(test_data_dir, "examples/COSMOS.para")

    with pytest.raises(ValueError) as excinfo:
        _ = MagSvc.from_config(objtype="UNKNOWN", config_file=config_file_path)
        assert excinfo.value == "Unknown objtype: UNKNOWN"
