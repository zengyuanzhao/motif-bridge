from io import StringIO

import pytest

from motif_bridge.io import is_logodds_row, logodds_to_prob, read_meme, write_meme


def test_logodds_to_prob_normalizes_rows():
    row = logodds_to_prob([2.0, 0.0, -1.0, -2.0], pseudocount=0.01, background=0.25)

    assert sum(row) == pytest.approx(1.0)
    assert row[0] > row[1] > row[2] > row[3]


def test_is_logodds_row_respects_explicit_input_format():
    assert is_logodds_row([0.2, 0.3, 0.3, 0.2], "probability") is False
    assert is_logodds_row([0.2, 0.3, 0.3, 0.2], "logodds") is True
    assert is_logodds_row([2.0, 0.0, -1.0, -2.0], "auto") is True
    assert is_logodds_row([0.25, 0.25, 0.25, 0.25], "auto") is False


def test_read_meme_warns_on_alength_conflict_and_keeps_alphabet_width(capsys):
    uniform_row = " ".join(["0.05"] * 20)
    mixed_row = " ".join(["0.10"] + ["0.05"] * 17 + ["0.00", "0.05"])
    source = f"""MEME version 4

ALPHABET= PROTEIN

MOTIF P0001 protein_motif

letter-probability matrix: alength= 4 w= 2 nsites= 20 E= 0
  {uniform_row}
  {mixed_row}
//
"""

    motifs = list(read_meme(StringIO(source)))

    captured = capsys.readouterr()
    assert "alength=4 conflicts with alphabet PROTEIN" in captured.err
    assert len(motifs) == 1
    assert motifs[0].alphabet == "PROTEIN"
    assert len(motifs[0].matrix) == 2
    assert all(len(row) == 20 for row in motifs[0].matrix)


def test_write_meme_omits_strands_for_protein():
    source = (
        "ALPHABET= PROTEIN\nMOTIF P1\nletter-probability matrix:\n" + " ".join(["0.05"] * 20) + "\n"
    )
    motifs = list(read_meme(StringIO(source)))
    output = StringIO()

    write_meme(motifs, output)

    rendered = output.getvalue()
    assert "ALPHABET= PROTEIN" in rendered
    assert "strands: + -" not in rendered
    assert "A 0.05" in rendered
