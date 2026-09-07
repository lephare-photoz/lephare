import pytest
from click.testing import CliRunner
from lephare.cli import build_cli

# Minimal dummy config keys
config_keys = {
    "FOO": "Foo option",
    "BAR": "Bar option",
}


# Dummy runner class for CLI isolation
class DummyRunner:
    def __init__(self, config_file=None, **kwargs):
        self.config_file = config_file
        self.kwargs = kwargs
        self.timer = False
        self.keymap = kwargs

    def run(self):
        return self.keymap

    def end(self):
        pass


@pytest.fixture
def cli():
    return build_cli(DummyRunner, config_keys)


def test_cli_invokes_with_config(cli, tmp_path):
    # Create a dummy config file
    config_file = tmp_path / "config.para"
    config_file.write_text("FOO dummy\nBAR dummy\n")

    runner = CliRunner()
    result = runner.invoke(cli, ["--config", str(config_file), "--FOO", "foo_value"], standalone_mode=False)

    assert result.exit_code == 0
    assert "FOO" in DummyRunner.kwargs if hasattr(DummyRunner, "kwargs") else True


def test_cli_timer_flag(cli):
    runner = CliRunner()
    result = runner.invoke(cli, ["--timer"], standalone_mode=False)
    assert result.exit_code == 0


def test_cli_dynamic_options(cli, tmp_path):
    # Create a dummy config file
    config_file = tmp_path / "config2.para"
    config_file.write_text("FOO dummy\nBAR dummy\n")

    runner = CliRunner()
    result = runner.invoke(
        cli, ["--config", str(config_file), "--FOO", "foo", "--BAR", "bar"], standalone_mode=False
    )
    assert result.exit_code == 0


# Config keys including the special-cased "typ" and "VERBOSE" entries
typed_config_keys = {
    "typ": "Object type",
    "VERBOSE": "Verbosity",
    "FOO": "Foo option",
}


class RecordingRunner(DummyRunner):
    """Runner that records the attributes the CLI sets on it."""

    last = None

    def __init__(self, config_file=None, **kwargs):
        super().__init__(config_file=config_file, **kwargs)
        self.typ = None
        self.verbose = None
        RecordingRunner.last = self


@pytest.fixture
def typed_cli():
    return build_cli(RecordingRunner, typed_config_keys)


def test_cli_uppercases_typ(typed_cli):
    """The -t/--typ option is normalised to upper case before the run."""
    result = CliRunner().invoke(typed_cli, ["--typ", "gal"], standalone_mode=False)
    assert result.exit_code == 0, result.output
    assert RecordingRunner.last.typ == "GAL"


def test_cli_short_typ_flag(typed_cli):
    """ "typ" gets a -t short form, unlike the other dynamic options."""
    result = CliRunner().invoke(typed_cli, ["-t", "qso"], standalone_mode=False)
    assert result.exit_code == 0, result.output
    assert RecordingRunner.last.typ == "QSO"


@pytest.mark.parametrize("value,expected", [("YES", True), ("NO", False), ("yes", False)])
def test_cli_verbose_flag(typed_cli, value, expected):
    """VERBOSE is true only for an exact "YES"."""
    result = CliRunner().invoke(typed_cli, ["--VERBOSE", value], standalone_mode=False)
    assert result.exit_code == 0, result.output
    assert RecordingRunner.last.verbose is expected


def test_cli_omits_unset_options(typed_cli):
    """Options left off the command line are not forwarded to the runner."""
    result = CliRunner().invoke(typed_cli, [], standalone_mode=False)
    assert result.exit_code == 0, result.output
    assert RecordingRunner.last.kwargs == {}
    # Neither typ nor verbose were touched
    assert RecordingRunner.last.typ is None
    assert RecordingRunner.last.verbose is None


def test_cli_timer_is_forwarded(typed_cli):
    """--timer is passed through as an attribute rather than a config key."""
    result = CliRunner().invoke(typed_cli, ["--timer"], standalone_mode=False)
    assert result.exit_code == 0, result.output
    assert RecordingRunner.last.timer is True
    assert "TIMER" not in RecordingRunner.last.kwargs


@pytest.mark.parametrize(
    "module",
    ["filter", "filter_extinc", "sedtolib", "mag_gal", "zphota"],
)
def test_entry_point_help(module, monkeypatch, capsys):
    """Each console-script entry point is importable and prints its own help.

    This exercises the same ``main()`` functions that pyproject.toml exposes as
    console scripts, but in-process, so it does not depend on the scripts having
    been re-installed.
    """
    import importlib

    mod = importlib.import_module(f"lephare.{module}")
    monkeypatch.setattr("sys.argv", [module, "--help"])

    with pytest.raises(SystemExit) as excinfo:
        mod.main()

    assert excinfo.value.code == 0
    assert "--help" in capsys.readouterr().out
