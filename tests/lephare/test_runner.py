import os

import lephare as lp
import pytest
from lephare._lephare import keyword

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


def test_runner_base():
    """This is the most rudimentary test of the Runner class. We have to provide
    both a list of config_keys and a path to a config file. We want to check only
    that the class can be instantiated and that the resulting keymap is not empty."""

    test_keys = {
        "CAT_IN": "help1",
    }
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    runner = lp.Runner(config_keys=test_keys, config_file=config_file_path)
    assert len(runner.keymap)
    runner = lp.Runner(config_keys=test_keys, config_keymap={"CAT_IN": keyword("CAT_IN", "dummy")})
    assert len(runner.keymap)
    runner = lp.Runner(config_keys=test_keys, CAT_IN="dummy")
    assert len(runner.keymap)
    runner = lp.Runner(config_keys=test_keys, cat_in="dummy")
    assert len(runner.keymap)

    for verbose_val in ["yes", "YES", True, 1]:
        r = lp.runner.Runner(config_keys={}, verbose=verbose_val)
        assert r.verbose
        assert r.keymap["VERBOSE"].value == "YES"

    for verbose_val in ["no", "NO", False, 0, "any"]:
        r = lp.runner.Runner(config_keys={}, verbose=verbose_val)
        assert not r.verbose
        assert r.keymap["VERBOSE"].value == "NO"


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
    test_keys = {"STAR_SED": "help1", "LIMITS_MAPP_CUT": "help2"}
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

    test_keys = {"key1": "help1", "key2": "help2", "key3": "help3"}
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    runner = lp.Runner(config_keys=test_keys, config_file=config_file_path)

    with pytest.raises(Exception) as excinfo:
        runner.run()
        assert excinfo.value == "runner.py is an abstract class"


def test_runner_verbosity():
    """Check to make sure that verbosity is set correctly via the config_keymap"""
    test_keys = {"key1": "help2", "key2": "help2", "key3": "help3"}
    config_file_path = os.path.join(TESTDATADIR, "examples/COSMOS.para")
    config_keymap = {}
    runner = lp.Runner(config_keys=test_keys, config_file=config_file_path, config_keymap=config_keymap)
    # VERBOSE is set to NO in COSMOS.para
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


def test_kwargs_arguments():
    test_keys = {"key1": "help1", "key2": "help2", "key3": "help3"}
    runner = lp.Runner(config_keys=test_keys)
    assert [k in runner.keymap for k in test_keys]

    runner = lp.Runner(config_keys={"A": "help"}, a="dummy2")
    assert runner.keymap["A"].value == "dummy2"

    with pytest.raises(RuntimeError) as excinfo:
        runner = lp.Runner(config_keys={"key1": "help"}, key2="unauthorized key")
        assert excinfo.value == f"key2 is not a recognized argument of {runner.__class__.__name__}."


def test_runner_timer_reports_elapsed(capsys):
    """With the timer on, end() prints how long the run took."""
    runner = lp.Runner(config_keys={"A": "help"})
    runner.timer = True
    runner.run()
    runner.end()

    out = capsys.readouterr().out
    assert "execution time:" in out


def test_runner_end_is_silent_without_timer(capsys):
    """Without the timer, end() prints nothing (and needs no start time)."""
    runner = lp.Runner(config_keys={"A": "help"})
    runner.run()
    runner.end()

    assert capsys.readouterr().out == ""


def test_runner_run_updates_keymap():
    """run() only copies through kwargs that are recognised config keys."""
    runner = lp.Runner(config_keys={"A": "help", "B": "help"})
    runner.run(a=1, b="two", unknown="ignored")

    assert runner.keymap["A"].value == "1"
    assert runner.keymap["B"].value == "two"
    assert "UNKNOWN" not in runner.keymap


def test_runner_run_verbose_kwarg():
    """VERBOSE passed to run() flips verbosity and is not stored as a config key."""
    runner = lp.Runner(config_keys={"A": "help"})
    assert not runner.verbose

    runner.run(VERBOSE="YES")
    assert runner.verbose

    runner.run(VERBOSE="NO")
    assert not runner.verbose


def test_runner_run_typ_kwarg():
    """typ is upper-cased by run(), but only for runners that declare it."""
    runner = lp.Runner(config_keys={"typ": "help", "A": "help"})
    runner.run(typ="gal")
    assert runner.typ == "GAL"

    # A runner without a "typ" config key leaves the attribute alone
    other = lp.Runner(config_keys={"A": "help"})
    other.run(typ="gal")
    assert other.typ is None


def test_runner_typ_key_is_accepted():
    """TYP is allowed as a kwarg when the runner declares a lower-case "typ" key."""
    runner = lp.Runner(config_keys={"typ": "help"})
    runner.validate_config_dict({"TYP": "GAL"}, no_raise=False)

    # ...but not when the runner has no typ key at all
    with pytest.raises(RuntimeError, match="TYP is not a recognized argument"):
        lp.Runner(config_keys={"A": "help"}).validate_config_dict({"TYP": "GAL"}, no_raise=False)


def test_runner_config_file_skips_comments_and_short_lines(tmp_path):
    """Comments, blank lines and one-token lines are ignored when reading a .para."""
    config_file = tmp_path / "sparse.para"
    config_file.write_text(
        "# a comment\n"
        "\n"
        "   \n"
        "LONELY\n"  # only one token, so no value
        "A first_value trailing comment text\n"
        "IGNORED not_a_config_key\n"
    )
    runner = lp.Runner(config_keys={"A": "help", "B": "help"}, config_file=str(config_file))

    assert runner.keymap["A"].value == "first_value"
    # Keys not present in the file still get an empty keyword
    assert runner.keymap["B"].value == ""
    # Keys in the file but not in config_keys are dropped
    assert "IGNORED" not in runner.keymap
    assert "LONELY" not in runner.keymap
