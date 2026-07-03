"""Shared helpers for the Python command-line entry points."""

import argparse
import gzip
import sys
from typing import TextIO

from motif_bridge.core import Background, _background_values


def background_prob(value: str) -> Background:
    """Parse and validate a scalar or comma-separated background vector."""
    values = []
    for part in value.split(","):
        try:
            values.append(float(part))
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"Invalid -b value: {value}") from exc

    background: Background = values if len(values) > 1 else values[0]
    try:
        _background_values(background, len(values))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return background


def open_input(path: str) -> TextIO:
    """Open plain, gzip, or stdin input."""
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")
