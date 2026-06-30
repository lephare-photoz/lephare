import subprocess


def test_command_line_tasks():
    """Simply check that they are there and return basic outputs"""
    for task in ["filter", "mag_gal", "sedtolib", "zphota", "filter_extinc"]:
        result = subprocess.run(
            [task, "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert result.stdout.startswith("User defined")
