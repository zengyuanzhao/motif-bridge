import math
import sys
from typing import Any, Dict, List


class Motif:
    def __init__(
        self, id: str, description: str, matrix: List[List[float]], alphabet: str = "ACGT"
    ):
        self.id = id
        self.description = description
        self.matrix = matrix
        self.alphabet = alphabet

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
        return {
            "id": self.id,
            "description": self.description,
            "matrix": [row[:] for row in self.matrix],
            "alphabet": self.alphabet,
        }

    def calculate_score(self, bg: float, t_offset: float) -> float:
        """Calculate HOMER log-odds threshold from a probability matrix."""
        score = 0.0
        for row in self.matrix:
            max_p = max(row) if row else 0.0
            if max_p > 0:
                score += math.log2(max_p / bg)
        score -= t_offset
        return max(score, 0.0)
