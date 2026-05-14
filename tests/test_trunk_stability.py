"""The architectural contract: the trunk must remain operable even when leaf
dependencies are missing or broken. These tests enforce the principle.
"""

import subprocess
import sys
import textwrap

import pytest


@pytest.mark.quick
def test_tiny_imports_cleanly():
    """Bare `import tiny` succeeds with no side effects."""
    result = subprocess.run(
        [sys.executable, "-c", "import tiny; print('OK')"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


@pytest.mark.quick
def test_tiny_help_does_not_import_pysam():
    """Running `tiny --help` must not trigger pysam import (architectural contract)."""
    code = textwrap.dedent("""
        import sys
        from typer.testing import CliRunner
        from tiny.cli import app
        runner = CliRunner()
        result = runner.invoke(app, ['--help'])
        assert result.exit_code == 0, result.stdout
        assert 'pysam' not in sys.modules, 'pysam was imported during --help!'
        print('OK')
    """)
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


@pytest.mark.quick
def test_tiny_analyze_works_without_pysam():
    """`tiny analyze ATCG` runs even when pysam is unavailable."""
    code = textwrap.dedent("""
        import sys
        sys.modules['pysam'] = None  # simulate pysam being uninstallable
        from typer.testing import CliRunner
        from tiny.cli import app
        runner = CliRunner()
        result = runner.invoke(app, ['analyze', 'ATCG'])
        assert result.exit_code == 0, result.stdout
        print('OK')
    """)
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout
