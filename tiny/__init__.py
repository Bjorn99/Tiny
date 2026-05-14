"""Tiny — a focused bioinformatics CLI for classical DNA sequence algorithms."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("tiny")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"
