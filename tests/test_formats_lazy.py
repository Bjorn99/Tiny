import subprocess
import sys
import textwrap
from pathlib import Path

import pytest


@pytest.mark.quick
def test_importing_formats_does_not_import_pysam():
    """The trunk must not pull in pysam just by being imported."""
    code = textwrap.dedent("""
        import sys
        import tiny.core.formats  # noqa: F401
        assert "pysam" not in sys.modules, f"pysam was imported! Modules: {[m for m in sys.modules if 'pysam' in m]}"
        print("OK")
    """)
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


@pytest.mark.quick
def test_read_sam_bam_raises_optional_dependency_error_when_pysam_missing(monkeypatch):
    """Calling read_sam_bam without pysam installed raises a typed, helpful error."""
    monkeypatch.setitem(sys.modules, "pysam", None)
    from tiny.core.errors import OptionalDependencyError
    from tiny.core.formats import FormatHandler

    with pytest.raises(OptionalDependencyError) as exc:
        FormatHandler.read_sam_bam(Path("/nonexistent.bam"))
    assert "pysam" in str(exc.value)
    assert "poetry install -E sam" in str(exc.value)
