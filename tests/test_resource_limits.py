import pytest

from tiny.core.errors import ResourceLimitError
from tiny.core.sequence import DNASequence


@pytest.mark.quick
def test_sequence_within_limit_accepted():
    DNASequence("A" * 10000)  # at default limit


@pytest.mark.quick
def test_sequence_above_limit_rejected():
    with pytest.raises(ResourceLimitError) as exc:
        DNASequence("A" * 10001)
    assert "10001" in str(exc.value)
    assert "10000" in str(exc.value)


@pytest.mark.quick
def test_limit_overridable_via_env(monkeypatch):
    monkeypatch.setenv("TINY_MAX_SEQUENCE", "20000")
    DNASequence("A" * 20000)  # should not raise
