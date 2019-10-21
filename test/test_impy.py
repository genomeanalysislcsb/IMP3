import click
import os
from click.testing import CliRunner

runner = CliRunner()
from impy import *

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def test_configuration_path_not_exists():
    result = runner.invoke(cli, ['-c', 'doesnotexists', 'init'])
    assert result.exit_code != 0


def test_configuration_path_isdir():
    result = runner.invoke(cli, ['-c', '.', 'init'])
    assert result.exit_code != 0


def test_configuration_path_wrong_file():
    result = runner.invoke(cli, ['-c', os.path.join(BASE_DIR, 'MT.R1.small.fq'), 'init'])
    assert result.exit_code != 0


def test_screen_path_not_exists():
    result = runner.invoke(cli, ['init', '--screen', 'doesnotexists'])
    assert result.exit_code != 0


def test_screen_path_isdir():
    result = runner.invoke(cli, ['init', '--screen', '.'])
    assert result.exit_code != 0
