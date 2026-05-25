import json
import sys
import math
from typing import List, Iterator, Iterable, TextIO

from .core import Motif

ALPHABETS = {
    "ACGT": "ACGT",
    "ACGU": "ACGU",
    "PROTEIN": "ACDEFGHIKLMNPQRSTVWY",
}

def _alphabet_letters(alphabet: str) -> str:
    return ALPHABETS.get(alphabet, alphabet)

def _format_background_line(alphabet: str) -> str:
    letters = _alphabet_letters(alphabet)
    if not letters:
        return ""
    freq = 1.0 / len(letters)
    fmt = "{:.2f}" if len(letters) in (4, 20) else "{:.6f}"
    parts = [f"{letter} {fmt.format(freq)}" for letter in letters]
    return " ".join(parts)

def _json_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=True)

def logodds_to_prob(
    row: List[float], pseudocount: float = 0.01, background: float = 0.25
) -> List[float]:
    """Convert a HOMER log-odds row to a probability row."""
    raw = [2**v * background for v in row]
    total = sum(raw) + pseudocount * len(raw)
    return [(v + pseudocount) / total for v in raw]

def is_logodds_row(row: List[float], input_format: str) -> bool:
    """Determine if row is log-odds based on input_format setting."""
    if input_format == "logodds":
        return True
    if input_format == "probability":
        return False
    s = sum(row)
    return not (0.98 <= s <= 1.02)

def read_meme(fh: TextIO, alphabet: str = "ACGT") -> Iterator[Motif]:
    """Parse MEME format from a file-like object and yield Motif instances."""
    in_motif = False
    in_matrix = False
    motif_id = ""
    description = ""
    matrix: List[List[float]] = []
    expected_cols = len(_alphabet_letters(alphabet))

    for raw_line in fh:
        line = raw_line.rstrip("\n")
        stripped = line.strip()

        if stripped.startswith("MOTIF"):
            if len(stripped) > 5 and not stripped[5].isspace():
                continue
            if in_motif and matrix:
                yield Motif(motif_id, description, matrix, alphabet)

            parts = stripped.split()
            mid = parts[1] if len(parts) > 1 else ""
            if not mid:
                sys.stderr.write(f"Warning: skipping malformed MOTIF line without ID: {stripped}\n")
                in_motif = False
                in_matrix = False
                matrix = []
                continue
            original_name = " ".join(parts[2:]) if len(parts) > 2 else (mid or "motif")

            motif_id = mid
            description = original_name
            matrix = []
            in_motif = True
            in_matrix = False
            continue

        if not in_motif:
            continue

        if stripped.startswith("URL"):
            continue
        if stripped.startswith("//"):
            if matrix:
                yield Motif(motif_id, description, matrix, alphabet)
            matrix = []
            in_motif = False
            in_matrix = False
            continue

        if stripped.startswith("letter-probability matrix:"):
            in_matrix = True
            if "alength=" in stripped:
                try:
                    alength = int(stripped.split("alength=")[1].split()[0])
                    expected_cols = alength
                except (ValueError, IndexError):
                    pass
            continue

        if in_matrix and stripped and (
            stripped[0].isdigit() or stripped.startswith(".") or stripped.startswith("-")
        ):
            tokens = stripped.split()
            try:
                row = [float(t) for t in tokens]
                if len(row) == expected_cols:
                    if any(v < 0 for v in row):
                        sys.stderr.write(
                            f"Warning: negative value in matrix row (expected probabilities): {stripped}\n"
                        )
                    matrix.append(row)
                elif row:
                    sys.stderr.write(
                        f"Warning: skipping malformed matrix row "
                        f"(expected {expected_cols} cols, got {len(row)}): {stripped}\n"
                    )
            except ValueError:
                pass

    if in_motif and matrix:
        yield Motif(motif_id, description, matrix, alphabet)

def read_homer(
    fh: TextIO,
    pseudocount: float = 0.01,
    input_format: str = "auto",
    alphabet: str = "ACGT",
    background: float = 0.25,
) -> Iterator[Motif]:
    """Parse HOMER format from a file-like object and yield Motif instances."""
    in_motif = False
    motif_id = ""
    description = ""
    matrix: List[List[float]] = []

    expected_cols = len(_alphabet_letters(alphabet))

    for raw_line in fh:
        line = raw_line.rstrip("\n")
        stripped = line.strip()

        if not stripped:
            continue

        if stripped.startswith(">"):
            if in_motif and matrix:
                yield Motif(motif_id, description, matrix, alphabet)

            parts = stripped[1:].split("\t")
            mid = parts[0] if len(parts) > 0 else "motif"
            desc = parts[1] if len(parts) > 1 else mid

            motif_id = mid
            description = desc
            matrix = []
            in_motif = True
            continue

        if not in_motif:
            continue

        tokens = stripped.split()
        try:
            row = [float(t) for t in tokens]
            if not row:
                continue
            if len(row) != expected_cols:
                sys.stderr.write(
                    f"Warning: skipping malformed row (expected {expected_cols} cols, got {len(row)}): {stripped}\n"
                )
                continue
            if is_logodds_row(row, input_format):
                row = logodds_to_prob(row, pseudocount, background)
            matrix.append(row)
        except ValueError:
            pass

    if in_motif and matrix:
        yield Motif(motif_id, description, matrix, alphabet)

def read_json(
    fh: TextIO, pseudocount: float = 0.01, input_format: str = "auto", background: float = 0.25
) -> Iterator[Motif]:
    """Parse JSON format from a file-like object and yield Motif instances."""
    data = json.load(fh)
    for motif in data.get("motifs", []):
        mid = motif.get("id", "motif")
        desc = motif.get("description", mid)
        matrix = motif.get("matrix", [])
        alphabet = motif.get("alphabet", "ACGT")

        if not matrix:
            continue

        expected_cols = len(_alphabet_letters(alphabet))
        processed_matrix = []
        for row in matrix:
            if len(row) != expected_cols:
                sys.stderr.write(
                    f"Warning: skipping malformed matrix row (expected {expected_cols} cols, got {len(row)})\n"
                )
                continue
            if is_logodds_row(row, input_format):
                row = logodds_to_prob(row, pseudocount, background)
            processed_matrix.append(row)

        if processed_matrix:
            yield Motif(mid, desc, processed_matrix, alphabet)

def _calculate_homer_score(matrix: List[List[float]], background: float, threshold_offset: float) -> float:
    score = 0.0
    for row in matrix:
        max_p = max(row) if row else 0.0
        if max_p > 0:
            score += math.log2(max_p / background)
    score -= threshold_offset
    return max(score, 0.0)

def write_homer(motifs: Iterable[Motif], fh: TextIO, background: float = 0.25, threshold_offset: float = 4.0) -> None:
    """Write an iterable of Motif objects to a file-like object in HOMER format."""
    for m in motifs:
        score = _calculate_homer_score(m.matrix, background, threshold_offset)
        fh.write(f">{m.id}\t{m.description}\t{score:.6f}\t0\t0\t0\n")
        for row in m.matrix:
            fh.write("\t".join(f"{v:.6f}" for v in row) + "\n")

def write_meme(motifs: Iterable[Motif], fh: TextIO) -> None:
    """Write an iterable of Motif objects to a file-like object in MEME format."""
    header_printed = False
    header_alphabet = ""
    
    for m in motifs:
        if not header_printed:
            header_alphabet = m.alphabet
            fh.write("MEME version 4\n\n")
            fh.write(f"ALPHABET= {m.alphabet}\n\n")
            fh.write("strands: + -\n\n")
            fh.write("Background letter frequencies\n")
            fh.write(_format_background_line(m.alphabet) + "\n\n")
            header_printed = True
        elif m.alphabet != header_alphabet:
            sys.stderr.write(
                f"Warning: skipping motif '{m.id}' with alphabet {m.alphabet} "
                f"(header uses {header_alphabet})\n"
            )
            continue

        width = len(m.matrix)
        expected_cols = len(ALPHABETS.get(m.alphabet, m.alphabet))
        fh.write(f"MOTIF {m.id} {m.description}\n\n")
        fh.write(f"letter-probability matrix: alength= {expected_cols} w= {width} nsites= 20 E= 0\n")
        for row in m.matrix:
            fh.write("  " + "  ".join(f"{v:.6f}" for v in row) + "\n")
        fh.write("\n")

def write_json(motifs: Iterable[Motif], fh: TextIO) -> None:
    """Write an iterable of Motif objects to a file-like object in JSON format."""
    motifs_list = list(motifs)
    fh.write("{\n")
    fh.write('  "version": "1.0",\n')
    fh.write('  "source": "meme",\n')
    fh.write('  "motifs": [\n')
    for mi, motif in enumerate(motifs_list):
        fh.write("    {\n")
        fh.write(f'      "id": {_json_string(motif.id)},\n')
        fh.write(f'      "description": {_json_string(motif.description)},\n')
        if motif.alphabet:
            fh.write(f'      "alphabet": {_json_string(motif.alphabet)},\n')
        fh.write('      "matrix": [\n')
        for ri, row in enumerate(motif.matrix):
            values = ", ".join(f"{v:.6f}" for v in row)
            suffix = "," if ri + 1 < len(motif.matrix) else ""
            fh.write(f"        [{values}]{suffix}\n")
        fh.write("      ]\n")
        motif_suffix = "," if mi + 1 < len(motifs_list) else ""
        fh.write(f"    }}{motif_suffix}\n")
    fh.write("  ]\n")
    fh.write("}\n")
