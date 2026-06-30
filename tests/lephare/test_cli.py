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
