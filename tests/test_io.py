from io import StringIO

import pytest

from motif_bridge.io import (
    is_logodds_row,
    logodds_to_prob,
    read_homer,
    read_meme,
    write_homer,
    write_meme,
)


def test_logodds_to_prob_normalizes_rows():
    row = logodds_to_prob([2.0, 0.0, -1.0, -2.0], pseudocount=0.01, background=0.25)

    assert sum(row) == pytest.approx(1.0)
    assert row[0] > row[1] > row[2] > row[3]


def test_is_logodds_row_respects_explicit_input_format():
    assert is_logodds_row([0.2, 0.3, 0.3, 0.2], "probability") is False
    assert is_logodds_row([0.2, 0.3, 0.3, 0.2], "logodds") is True
    assert is_logodds_row([2.0, 0.0, -1.0, -2.0], "auto") is True
    assert is_logodds_row([0.25, 0.25, 0.25, 0.25], "auto") is False


def test_is_logodds_row_warns_on_nonnegative_auto_gray_zone(capsys):
    assert is_logodds_row([0.30, 0.30, 0.30, 0.15], "auto") is True

    captured = capsys.readouterr()
    assert "auto-detection boundary" in captured.err


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


def test_write_meme_metadata_overrides_and_renormalizes_rows():
    motif = read_homer(
        StringIO(">M1\tdesc\t1\t0\t0\t0\n0.20\t0.20\t0.20\t0.20\n"),
        input_format="probability",
    )
    output = StringIO()

    write_meme(motif, output, nsites=4000, evalue=0.00001, renormalize=True)

    rendered = output.getvalue()
    assert "nsites= 4000 E= 0.000010" in rendered
    assert "  0.250000  0.250000  0.250000  0.250000" in rendered


def test_write_homer_warns_on_zero_threshold(capsys):
    motif = read_meme(StringIO("MOTIF WEAK\nletter-probability matrix:\n0.25 0.25 0.25 0.25\n"))
    output = StringIO()

    write_homer(motif, output)

    captured = capsys.readouterr()
    assert "clipped to 0" in captured.err
    assert output.getvalue().startswith(">WEAK\tWEAK\t0.000000\t")


def test_write_homer_can_keep_source_threshold():
    motifs = read_homer(
        StringIO(">KEEP\tdesc\t9.5\t0\t0\t0\n0.25\t0.25\t0.25\t0.25\n"),
        input_format="probability",
    )
    output = StringIO()

    write_homer(motifs, output, keep_threshold=True)

    assert output.getvalue().splitlines()[0].split("\t")[2] == "9.500000"
