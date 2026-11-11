"""Compatibility utilities for Python 3.9 through 3.13+."""

from dataclasses import dataclass
import sys


def slotted_dataclass(*args, **kwargs):
    """Dataclass decorator supporting slots=True on Python 3.10+ while falling back cleanly on 3.9."""
    if sys.version_info < (3, 10):
        kwargs.pop("slots", None)
    return dataclass(*args, **kwargs)

