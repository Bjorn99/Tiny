import pytest

from tiny.core.errors import (
    InvalidSequenceError,
    OptionalDependencyError,
    ResourceLimitError,
    TinyError,
    UnsupportedFormatError,
)


@pytest.mark.quick
def test_all_inherit_from_tiny_error():
    for cls in (
        InvalidSequenceError,
        UnsupportedFormatError,
        OptionalDependencyError,
        ResourceLimitError,
    ):
        assert issubclass(cls, TinyError)


@pytest.mark.quick
def test_optional_dependency_error_carries_install_hint():
    err = OptionalDependencyError("pysam", "poetry install -E sam")
    assert "pysam" in str(err)
    assert "poetry install -E sam" in str(err)


@pytest.mark.quick
def test_resource_limit_error_carries_actual_and_limit():
    err = ResourceLimitError("sequence length", actual=20000, limit=10000)
    msg = str(err)
    assert "20000" in msg
    assert "10000" in msg
