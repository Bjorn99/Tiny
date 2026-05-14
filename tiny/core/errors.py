"""Typed exception hierarchy for Tiny.

Trunk and algorithm modules raise these; they never sys.exit or print.
Only tiny/cli.py translates exceptions into user-facing output.
"""


class TinyError(Exception):
    """Base class for all Tiny-raised errors."""


class InvalidSequenceError(TinyError):
    """A sequence string contained invalid characters or violated invariants."""


class UnsupportedFormatError(TinyError):
    """A file extension or format identifier is not handled by Tiny."""


class OptionalDependencyError(TinyError):
    """An optional dependency is missing. Carries an install hint."""

    def __init__(self, dependency: str, install_hint: str):
        super().__init__(
            f"Optional dependency '{dependency}' is not installed. " f"Install with: {install_hint}"
        )
        self.dependency = dependency
        self.install_hint = install_hint


class ResourceLimitError(TinyError):
    """An input exceeds a configured size/length limit."""

    def __init__(self, kind: str, actual: int, limit: int):
        super().__init__(
            f"{kind} exceeds limit: got {actual}, max allowed {limit}. "
            f"Override with environment variables if needed."
        )
        self.kind = kind
        self.actual = actual
        self.limit = limit
