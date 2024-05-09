import os

import lephare as lp
import pytest
from lephare._lephare import keyword

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


def test_runner_base():
    """This is the most rudimentary test of the Runner class. We have to provide
    both a list of config_keys and a path to a config file. We want to check only
    that the class can be instantiated and that the resulting keymap isn't empty."""

    test_keys = ["key1", "key2", "key3"]
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    runner = lp.Runner(config_keys=test_keys, config_file=config_file_path)
    assert len(runner.keymap)


def test_runner_no_config_keys():
    """Expect that a RuntimeError is raised when no config_keys or
    config_keys=None is passed to Runner."""
    with pytest.raises(RuntimeError) as excinfo:
        lp.Runner()
        assert excinfo.value == "Runner is a base class and cannot be initialized"

    with pytest.raises(RuntimeError) as excinfo:
        lp.Runner(config_keys=None)
        assert excinfo.value == "Runner is a base class and cannot be initialized"


def test_runner_config_keymap_updates():
    """This test checks that the keymap is updated when a config_keymap is passed
    to the Runner class."""
    test_keys = ["key1", "key2", "key3"]
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    config_keymap = {
        "STAR_SED": keyword("STAR_SED", "foo"),
        "LIMITS_MAPP_CUT": keyword("LIMITS_MAPP_CUT", "42"),
    }
    runner = lp.Runner(config_keys=test_keys, config_file=config_file_path, config_keymap=config_keymap)

    resulting_keymap = runner.keymap
    assert resulting_keymap["STAR_SED"].value == "foo"
    assert resulting_keymap["LIMITS_MAPP_CUT"].value == "42"


def test_runner_cannot_run():
    """Since the runner class is an abstract class, we expect that it cannot call
    the run method directly."""

    test_keys = ["key1", "key2", "key3"]
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    runner = lp.Runner(config_keys=test_keys, config_file=config_file_path)

    with pytest.raises(Exception) as excinfo:
        runner.run()
        assert excinfo.value == "runner.py is an abstract class"


def test_runner_verbosity():
    """Check to make sure that verbosity is set correctly via the config_keymap"""
    test_keys = ["key1", "key2", "key3"]
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    config_keymap = {}
    runner = lp.Runner(config_keys=test_keys, config_file=config_file_path, config_keymap=config_keymap)

    assert not runner.verbose

    config_keymap = {"VERBOSE": keyword("VERBOSE", "NO")}
    runner = lp.Runner(config_keys=test_keys, config_file=config_file_path, config_keymap=config_keymap)

    assert not runner.verbose

    config_keymap = {"VERBOSE": keyword("VERBOSE", "YES")}
    runner = lp.Runner(config_keys=test_keys, config_file=config_file_path, config_keymap=config_keymap)

    assert runner.verbose


def test_runner_config_file_not_found():
    """Pass in a config file that does not exist and expect a RuntimeError."""

    test_keys = ["key1", "key2", "key3"]
    config_file_path = os.path.join(TESTDATADIR, "foo/bar.para")

    with pytest.raises(RuntimeError) as excinfo:
        _ = lp.Runner(config_keys=test_keys, config_file=config_file_path)
        assert excinfo.value == f"File {config_file_path} not found"


def test_command_line_argument_parsing_basic(monkeypatch):
    """Check to make sure that command line arguments are parsed correctly."""
    test_keys = ["key1", "key2", "key3"]
    monkeypatch.setattr("sys.argv", ["runner.py", "--key1", "foo", "--key2", "42"])
    runner = lp.Runner(config_keys=test_keys)

    assert runner.keymap["key1"].value == "foo"
    assert runner.keymap["key2"].value == "42"
    assert runner.keymap["key3"].value == ""
    assert len(runner.keymap) == 3
    assert runner.args.config == ""
    assert runner.args.timer is False
    assert runner.args.verbose is False


def test_command_line_argument_parsing_with_known_args(monkeypatch):
    """Check to make sure that command line arguments are parsed correctly when
    including known arguments."""
    test_keys = ["key1", "key2", "key3"]
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    monkeypatch.setattr(
        "sys.argv", ["runner.py", "--key1", "foo", "--key2", "42", "--config", config_file_path, "--timer"]
    )
    runner = lp.Runner(config_keys=test_keys)

    assert runner.keymap["key1"].value == "foo"
    assert runner.keymap["key2"].value == "42"
    assert runner.keymap["key3"].value == ""
    assert runner.keymap["QSO_FSCALE"].value == "1."
    assert len(runner.keymap) > 3
    assert runner.args.config == config_file_path
    assert runner.args.timer is True
    assert runner.args.verbose is False
    with pytest.raises(AttributeError) as excinfo:
        _ = runner.args.typ
        assert excinfo.value == "'Runner' object has no attribute 'typ'"


def test_command_line_argument_parsing_with_subclass(monkeypatch):
    """Want to cover the case where the `typ` argument is passed to the runner."""

    class Sedtolib(lp.Runner):
        @property
        def __class__(self):
            return type("Sedtolib", (object,), {})

    test_keys = ["key1", "key2", "key3"]
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    monkeypatch.setattr(
        "sys.argv",
        [
            "runner.py",
            "--typ",
            "BAR",
            "--key1",
            "foo",
            "--key2",
            "42",
            "--config",
            config_file_path,
            "--timer",
        ],
    )
    runner = Sedtolib(config_keys=test_keys)

    assert runner.keymap["key1"].value == "foo"
    assert runner.keymap["key2"].value == "42"
    assert runner.keymap["key3"].value == ""
    assert runner.keymap["QSO_FSCALE"].value == "1."
    assert len(runner.keymap) > 3
    assert runner.args.config == config_file_path
    assert runner.args.timer is True
    assert runner.args.verbose is False
    assert runner.args.typ == "BAR"
    assert runner.typ == "BAR"
