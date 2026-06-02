import math

import pytest

from motif_bridge.core import Motif


def test_calculate_ic_dna_conserved_and_uniform_rows():
    motif = Motif("m1", "dna", [[1.0, 0.0, 0.0, 0.0], [0.25, 0.25, 0.25, 0.25]])

    ic = motif.calculate_ic()

    assert ic[0] == pytest.approx(2.0)
    assert ic[1] == pytest.approx(0.0)
    assert motif.total_ic() == pytest.approx(2.0)


def test_calculate_ic_protein_uses_twenty_letter_maximum():
    conserved = [1.0] + [0.0] * 19
    uniform = [0.05] * 20
    motif = Motif("p1", "protein", [conserved, uniform], alphabet="PROTEIN")

    ic = motif.calculate_ic()

    assert ic[0] == pytest.approx(math.log2(20))
    assert ic[1] == pytest.approx(0.0)


def test_reverse_complement_reverses_rows_and_swaps_columns():
    motif = Motif(
        "m1",
        "dna",
        [
            [0.1, 0.2, 0.3, 0.4],
            [0.5, 0.6, 0.7, 0.8],
        ],
    )

    motif.reverse_complement()

    assert motif.id == "m1_RC"
    assert motif.matrix == [
        [0.8, 0.7, 0.6, 0.5],
        [0.4, 0.3, 0.2, 0.1],
    ]


def test_reverse_complement_rejects_protein():
    motif = Motif("p1", "protein", [[1.0] + [0.0] * 19], alphabet="PROTEIN")

    with pytest.raises(ValueError, match="Reverse complement not supported"):
        motif.reverse_complement()


def test_trim_edges_removes_low_ic_flanks():
    motif = Motif(
        "m1",
        "dna",
        [
            [0.25, 0.25, 0.25, 0.25],
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.25, 0.25, 0.25, 0.25],
        ],
    )

    motif.trim_edges(0.5)

    assert motif.matrix == [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
    ]


def test_calculate_score_uses_column_specific_background():
    motif = Motif(
        "m1",
        "dna",
        [
            [0.1, 0.2, 0.3, 0.4],
            [0.6, 0.2, 0.1, 0.1],
        ],
    )

    score = motif.calculate_score([0.30, 0.20, 0.20, 0.30], 0.0)

    assert score == pytest.approx(math.log2(0.4 / 0.30) + math.log2(0.6 / 0.30))


def test_calculate_score_rejects_background_width_mismatch():
    motif = Motif("m1", "dna", [[0.25, 0.25, 0.25, 0.25]])

    with pytest.raises(ValueError, match="background length"):
        motif.calculate_score([0.25, 0.25], 0.0)
