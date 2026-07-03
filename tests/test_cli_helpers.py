import argparse
from io import StringIO

import pytest

from python_scripts._cli import background_prob, open_input


def test_background_prob_parses_scalar() -> None:
    assert background_prob("0.25") == pytest.approx(0.25)


def test_background_prob_parses_vector() -> None:
    assert background_prob("0.29,0.21,0.21,0.29") == pytest.approx([0.29, 0.21, 0.21, 0.29])


def test_background_prob_reuses_core_validation() -> None:
    with pytest.raises(argparse.ArgumentTypeError, match="sum to 1.0"):
        background_prob("0.5,0.5,0.5,0.5")


def test_open_input_returns_stdin(monkeypatch: pytest.MonkeyPatch) -> None:
    stdin = StringIO("motif\n")
    monkeypatch.setattr("sys.stdin", stdin)

    assert open_input("-") is stdin
