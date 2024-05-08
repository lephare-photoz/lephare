import os

import lephare as lp
import pytest

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
    config_keymap = {"STAR_SED": "foo", "LIMITS_MAPP_CUT": 42}
    runner = lp.Runner(config_keys=test_keys, config_file=config_file_path, config_keymap=config_keymap)

    resulting_keymap = runner.keymap
    assert resulting_keymap["STAR_SED"] == "foo"
    assert resulting_keymap["LIMITS_MAPP_CUT"] == 42


def test_runner_cannot_run():
    """Since the runner class is an abstract class, we expect that it cannot call
    the run method directly."""

    test_keys = ["key1", "key2", "key3"]
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    runner = lp.Runner(config_keys=test_keys, config_file=config_file_path)

    with pytest.raises(Exception) as excinfo:
        runner.run()
        assert excinfo.value == "runner.py is an abstract class"
