import os
from unittest.mock import ANY, MagicMock, patch

import pytest
from lephare import Zphota

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


# Test init


def test_no_config():
    with pytest.raises(SystemExit):
        _ = Zphota()


def test_with_config_file():
    """Use a config file to instantiate the Zphota class. Expect the keymap to be
    populated."""
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    zp = Zphota(config_file=config_file_path)
    assert len(zp.keymap)
    assert zp.keymap["CAT_IN"].value == "bidon"


def test_zphota_config_file_not_found():
    """Use a non-existent config file to instantiate the Zphota class. Expect a
    RuntimeError."""
    config_file_path = os.path.join(TESTDATADIR, "foo/bar.para")

    with pytest.raises(RuntimeError) as excinfo:
        _ = Zphota(config_file=config_file_path)
        assert excinfo.value == f"File {config_file_path} not found"


def test_with_command_line_config(monkeypatch):
    """Use command line arguments to instantiate the Zphota class. Expect the keymap
    to be populated."""
    monkeypatch.setattr("sys.argv", ["zphota.py", "--CAT_IN", "cat_foo", "--Z_INTERP", "YES"])
    zp = Zphota()
    assert zp.keymap["CAT_IN"].value == "cat_foo"
    assert zp.keymap["Z_INTERP"].value == "YES"


def test_with_command_line_unrecogized_config(monkeypatch):
    """Use command line arguments to instantiate the Zphota class, but include
    an unrecognized argument."""
    with pytest.raises(SystemExit):
        monkeypatch.setattr("sys.argv", ["zphota.py", "--FOO", "bar"])
        _ = Zphota()


def test_override_config_file_with_command_line(monkeypatch):
    """Use the command line to feed in a config file, but override one of the values
    with a command line value.

     This is the expected behavior, but only if we pass the config file in on the
    command line. Otherwise, as we test in the last part of this, we do not override
    the value."""
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    monkeypatch.setattr("sys.argv", ["zphota.py", "--config", config_file_path, "--CAT_IN", "cat_foo"])
    zp = Zphota()
    assert zp.keymap["CAT_IN"].value == "cat_foo"
    assert zp.keymap["Z_INTERP"].value == "YES"


# Test run


def run_configured_zp(zp, config_file_path):
    # Create mock onesource objects
    mock_sources = [
        MagicMock(spec=["spec", "zs"]),
        MagicMock(spec=["spec", "zs"]),
        MagicMock(spec=["spec", "zs"]),
    ]

    # Set the attributes
    mock_sources[0].spec = "spec1"
    mock_sources[0].zs = 0.1
    mock_sources[1].spec = "spec2"
    mock_sources[1].zs = 0.2
    mock_sources[2].spec = "spec3"
    mock_sources[2].zs = 0.3

    # Run zphota with cpp functions mocked
    with patch(
        "lephare._lephare.PhotoZ.read_photoz_sources", return_value=mock_sources
    ) as mock_read_photoz_sources:
        with patch("lephare._lephare.PhotoZ.run_photoz") as mock_run_photoz:
            with patch("lephare._lephare.PhotoZ.write_outputs") as mock_write_outputs:
                mock_run_photoz.return_value = None
                mock_write_outputs.return_value = None

                zp.run()

                mock_read_photoz_sources.assert_called_once()
                mock_run_photoz.assert_called_once()
                mock_run_photoz.assert_called_with(mock_sources, ANY, ANY)
                mock_write_outputs.assert_called_once()
                mock_write_outputs.assert_called_with(mock_sources, ANY)

                assert os.path.realpath(zp.keymap["c"].value) == os.path.realpath(config_file_path)


def test_run_zp_config_file(monkeypatch):
    """Use a config file to run the Zphota class (ie, running from a python session).

    Note that we mock the cpp functions that are called by the run method."""
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    zp = Zphota(config_file=config_file_path)

    run_configured_zp(zp, config_file_path)


def test_run_zp_command_line(monkeypatch):
    """Use the command line to run the Zphota class."""
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    monkeypatch.setattr("sys.argv", ["zphota.py", "--config", config_file_path, "--timer"])
    zp = Zphota()

    run_configured_zp(zp, config_file_path)


def test_run_zp_extra_command_line_flags(monkeypatch):
    """Run the zphota class with timer and verbose flags."""
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    monkeypatch.setattr("sys.argv", ["zphota.py", "--config", config_file_path, "--timer", "--verbose"])
    zp = Zphota(config_file=config_file_path)

    run_configured_zp(zp, config_file_path)


def test_run_zp_auto_adapt_yes(monkeypatch):
    """Run the zphota class with AUTO_ADAPT toggled on."""
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    monkeypatch.setattr("sys.argv", ["zphota.py", "--config", config_file_path, "--AUTO_ADAPT", "YES"])
    zp = Zphota()

    run_configured_zp(zp, config_file_path)


def test_run_zp_auto_adapt_no(monkeypatch):
    """Run the zphota class with AUTO_ADAPT toggled off."""
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    monkeypatch.setattr("sys.argv", ["zphota.py", "--config", config_file_path, "--AUTO_ADAPT", "NO"])
    zp = Zphota()

    run_configured_zp(zp, config_file_path)
