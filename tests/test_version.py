import pytest
from typer.testing import CliRunner

from tiny.cli import app


@pytest.mark.quick
def test_version_flag_prints_version():
    runner = CliRunner()
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert "0.2.0" in result.stdout


@pytest.mark.quick
def test_tiny_module_exposes_version():
    import tiny

    assert hasattr(tiny, "__version__")
    assert tiny.__version__ == "0.2.0"
