# cli.py (or at bottom of filter.py)
import click


def build_cli(runner_class, config_keys):
    @click.command(context_settings={"help_option_names": ["-h", "--help"]})
    @click.option(
        "-c",
        "--config",
        type=click.Path(exists=True),
        help="Path to ASCII config file",
    )
    @click.option("--timer", is_flag=True, help="Switch timer on")
    def cli(config, timer, **kwargs):
        normalized = {k.upper(): v for k, v in kwargs.items() if v is not None}
        runner = runner_class(
            config_file=config,
            **normalized,
        )
        runner.timer = timer
        if "TYP" in normalized:
            runner.typ = normalized["TYP"]
        if "VERBOSE" in normalized:
            runner.verbose = normalized["VERBOSE"] == "YES"
        runner.run()
        runner.end()

    # dynamic options
    for key, help_text in config_keys.items():
        if key == "typ":
            cli = click.option(
                "-t",
                "--typ",
                type=str,
                help=help_text,
            )(cli)
        else:
            cli = click.option(
                f"--{key}",
                type=str,
                help=help_text,
            )(cli)

    return cli
