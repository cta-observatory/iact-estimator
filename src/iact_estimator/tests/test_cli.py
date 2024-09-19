"""Tests for the Command-line interface."""


def test_cli_version(script_runner):
    from .. import __version__

    result = script_runner.run(["iact-estimator", "--version"])
    assert result.returncode == 0
    assert result.stderr == f"{__version__}\n"


def test_cli_config(script_runner, tmp_path):
    result = script_runner.run(["iact-estimator", "config", "--to", tmp_path])
    assert result.returncode == 0
    assert result.stderr == "You're halfway there! Have fun!\n"

    assert (tmp_path / "config.yml").is_file()


def test_cli_run(script_runner, tmp_path):
    script_runner.run(["iact-estimator", "config", "--to", tmp_path])
    result = script_runner.run(
        [
            "iact-estimator",
            "run",
            "--config",
            tmp_path / "config.yml",
            "--source-name",
            "Crab Nebula",
            "--output-path",
            tmp_path,
        ],
    )
    assert result.returncode == 0
