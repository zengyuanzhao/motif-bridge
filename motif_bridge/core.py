import math
import sys
from typing import Any, Dict, List, Optional, Sequence, Union

Background = Union[float, Sequence[float]]


def _background_values(background: Background, width: int) -> List[float]:
    if isinstance(background, (list, tuple)):
        values = [float(v) for v in background]
        if len(values) != width:
            raise ValueError(f"background length {len(values)} does not match row width {width}")
    else:
        values = [float(background)] * width
    if any(v <= 0 for v in values):
        raise ValueError("background values must be > 0")
    return values


def _renormalized(row: List[float]) -> List[float]:
    total = sum(row)
    return [v / total for v in row] if total > 0 else row[:]


class Motif:
    def __init__(
        self,
        id: str,
        description: str,
        matrix: List[List[float]],
        alphabet: str = "ACGT",
        threshold: Optional[float] = None,
    ):
        self.id = id
        self.description = description
        self.matrix = matrix
        self.alphabet = alphabet
        self.threshold = threshold

    def calculate_ic(self) -> List[float]:
        """Calculate Information Content for each position."""
        ic_list = []
        is_protein = self.alphabet == "PROTEIN"
        max_ic = math.log2(20) if is_protein else 2.0

        for row in self.matrix:
            h = 0.0
            for p in row:
                if p > 0:
                    h -= p * math.log2(p)
            ic = max_ic - h
            ic_list.append(max(0.0, ic))
        return ic_list

    def total_ic(self) -> float:
        """Return the total information content of the motif."""
        return sum(self.calculate_ic())

    def trim_edges(self, threshold: float) -> None:
        """Trim positions from both edges that have Information Content < threshold."""
        ic_list = self.calculate_ic()

        start = 0
        while start < len(ic_list) and ic_list[start] < threshold:
            start += 1

        end = len(ic_list)
        while end > start and ic_list[end - 1] < threshold:
            end -= 1

        if start < end:
            self.matrix = self.matrix[start:end]
        else:
            sys.stderr.write(
                f"Warning: motif '{self.id}' trimmed to empty matrix (IC threshold={threshold})\n"
            )
            self.matrix = []

    def reverse_complement(self) -> None:
        """Reverse complement the motif matrix in-place. Only valid for DNA/RNA."""
        if self.alphabet not in ("ACGT", "ACGU"):
            raise ValueError(f"Reverse complement not supported for alphabet: {self.alphabet}")
        if not self.matrix:
            return
        expected = len(self.matrix[0])
        if expected != 4:
            raise ValueError(f"Reverse complement requires 4 columns, got {expected}")

        rev_matrix = self.matrix[::-1]
        rc_matrix = []
        for row in rev_matrix:
            # Swap A(0)<->T/U(3) and C(1)<->G(2)
            rc_row = [row[3], row[2], row[1], row[0]]
            rc_matrix.append(rc_row)

        self.matrix = rc_matrix
        self.id = self.id + "_RC"

    def to_dict(self) -> Dict[str, Any]:
        """Serialize motif to dictionary."""
        data = {
            "id": self.id,
            "description": self.description,
            "matrix": [row[:] for row in self.matrix],
            "alphabet": self.alphabet,
        }
        if self.threshold is not None:
            data["threshold"] = self.threshold
        return data

    def calculate_score(
        self, bg: Background, t_offset: float, renormalize: bool = False
    ) -> float:
        """Calculate HOMER log-odds threshold from a probability matrix."""
        score = 0.0
        for row in self.matrix:
            if not row:
                continue
            values = _renormalized(row) if renormalize else row
            backgrounds = _background_values(bg, len(values))
            best_idx = max(range(len(values)), key=values.__getitem__)
            max_p = values[best_idx]
            if max_p > 0:
                score += math.log2(max_p / backgrounds[best_idx])
        score -= t_offset
        return max(score, 0.0)
