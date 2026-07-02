import json
import sys
from typing import Iterable, Iterator, List, Optional, TextIO

from .core import Background, Motif, _background_values, _renormalized

ALPHABETS = {
    "ACGT": "ACGT",
    "ACGU": "ACGU",
    "PROTEIN": "ACDEFGHIKLMNPQRSTVWY",
}
PROBABILITY_SUM_TOLERANCE = 1e-3


def _alphabet_letters(alphabet: str) -> str:
    return ALPHABETS.get(alphabet, alphabet)


def _format_background_line(alphabet: str, background: Optional[Background] = None) -> str:
    letters = _alphabet_letters(alphabet)
    if not letters:
        return ""
    if isinstance(background, (list, tuple)) and len(background) > 1:
        values = _background_values(background, len(letters))
        parts = [f"{letter} {value:.6f}" for letter, value in zip(letters, values)]
        return " ".join(parts)
    freq = 1.0 / len(letters)
    fmt = "{:.2f}" if len(letters) in (4, 20) else "{:.6f}"
    parts = [f"{letter} {fmt.format(freq)}" for letter in letters]
    return " ".join(parts)


def _json_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def _without_warning_prefix(message: str) -> str:
    prefix = "Warning: "
    return message[len(prefix) :] if message.startswith(prefix) else message


def _validate_probability_row(row: List[float], context: str) -> None:
    if any(v < 0 or v > 1 for v in row):
        raise ValueError(f"{context} values must be in [0, 1]")
    total = sum(row)
    if abs(total - 1.0) > PROBABILITY_SUM_TOLERANCE:
        raise ValueError(f"{context} must sum to 1.0, got {total:.6f}")


def logodds_to_prob(
    row: List[float], pseudocount: float = 0.01, background: Background = 0.25
) -> List[float]:
    """Convert a HOMER log-odds row to a probability row."""
    bg_values = _background_values(background, len(row))

    raw = [2**v * bg_values[i] for i, v in enumerate(row)]
    total = sum(raw) + pseudocount * len(raw)
    return [(v + pseudocount) / total for v in raw]


def is_logodds_row(row: List[float], input_format: str, strict: bool = False) -> bool:
    """Determine if row is log-odds based on input_format setting."""
    if input_format == "logodds":
        return True
    if input_format == "probability":
        return False
    s = sum(row)
    is_probability = 0.98 <= s <= 1.02
    if not is_probability and all(v >= 0 for v in row) and 0.90 <= s <= 1.10:
        message = (
            f"row sum {s:.4f} is close to the auto-detection boundary; "
            "use --input-format logodds or --input-format probability for reproducibility"
        )
        if strict:
            raise ValueError(message)
        sys.stderr.write(f"Warning: {message}\n")
    return not is_probability


def read_meme(
    fh: TextIO, alphabet_override: Optional[str] = None, strict: bool = False
) -> Iterator[Motif]:
    """Parse MEME format from a file-like object and yield Motif instances."""
    in_motif = False
    in_matrix = False
    motif_id = ""
    description = ""
    nsites: Optional[int] = None
    evalue: Optional[float] = None
    matrix: List[List[float]] = []
    alphabet = alphabet_override or "ACGT"
    expected_cols = len(_alphabet_letters(alphabet))

    for raw_line in fh:
        line = raw_line.rstrip("\n")
        stripped = line.strip()

        if stripped.startswith("ALPHABET="):
            if alphabet_override is None:
                detected = stripped.split("=", 1)[1].strip().split()[0]
                if detected:
                    alphabet = detected
                    expected_cols = len(_alphabet_letters(alphabet))
            continue

        if stripped.startswith("MOTIF"):
            if len(stripped) > 5 and not stripped[5].isspace():
                in_matrix = False
                continue
            if in_motif and matrix:
                yield Motif(motif_id, description, matrix, alphabet, nsites=nsites, evalue=evalue)

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
            nsites = None
            evalue = None
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
                yield Motif(motif_id, description, matrix, alphabet, nsites=nsites, evalue=evalue)
            matrix = []
            in_motif = False
            in_matrix = False
            continue

        if stripped.startswith("letter-probability matrix:"):
            in_matrix = True
            if "nsites=" in stripped:
                try:
                    nsites = int(stripped.split("nsites=")[1].split()[0])
                except (ValueError, IndexError):
                    nsites = None
            if "E=" in stripped:
                try:
                    evalue = float(stripped.split("E=")[1].split()[0])
                except (ValueError, IndexError):
                    evalue = None
            if "alength=" in stripped:
                try:
                    alength_token = stripped.split("alength=")[1].split()[0]
                except IndexError:
                    pass
                else:
                    try:
                        alength = int(alength_token)
                    except ValueError:
                        pass
                    else:
                        if alength != expected_cols:
                            message = (
                                f"alength={alength} conflicts with alphabet {alphabet} "
                                f"(expected {expected_cols} cols); using alphabet-derived width"
                            )
                            if strict:
                                raise ValueError(message)
                            sys.stderr.write(f"Warning: {message}\n")
                        else:
                            expected_cols = alength
            continue

        if (
            in_matrix
            and stripped
            and (stripped[0].isdigit() or stripped.startswith(".") or stripped.startswith("-"))
        ):
            tokens = stripped.split()
            try:
                row = [float(t) for t in tokens]
                if len(row) == expected_cols:
                    if strict:
                        _validate_probability_row(row, f"MEME matrix row for motif '{motif_id}'")
                    elif any(v < 0 for v in row):
                        msg = (
                            "Warning: negative value in matrix row "
                            f"(expected probabilities): {stripped}\n"
                        )
                        sys.stderr.write(msg)
                    matrix.append(row)
                elif row:
                    message = (
                        "Warning: skipping malformed matrix row "
                        f"(expected {expected_cols} cols, got {len(row)}): {stripped}\n"
                    )
                    if strict:
                        raise ValueError(_without_warning_prefix(message).rstrip())
                    sys.stderr.write(message)
            except ValueError:
                if strict:
                    raise
                pass

    if in_motif and matrix:
        yield Motif(motif_id, description, matrix, alphabet, nsites=nsites, evalue=evalue)


def read_homer(
    fh: TextIO,
    pseudocount: float = 0.01,
    input_format: str = "auto",
    alphabet: str = "ACGT",
    background: Background = 0.25,
    strict: bool = False,
) -> Iterator[Motif]:
    """Parse HOMER format from a file-like object and yield Motif instances."""
    in_motif = False
    motif_id = ""
    description = ""
    threshold: Optional[float] = None
    matrix: List[List[float]] = []

    expected_cols = len(_alphabet_letters(alphabet))

    for raw_line in fh:
        line = raw_line.rstrip("\n")
        stripped = line.strip()

        if not stripped:
            continue

        if stripped.startswith(">"):
            if in_motif and matrix:
                yield Motif(motif_id, description, matrix, alphabet, threshold=threshold)

            parts = stripped[1:].split("\t")
            mid = parts[0] if len(parts) > 0 else "motif"
            desc = parts[1] if len(parts) > 1 else mid
            source_threshold = None
            if len(parts) > 2:
                try:
                    source_threshold = float(parts[2])
                except ValueError:
                    source_threshold = None

            motif_id = mid
            description = desc
            threshold = source_threshold
            matrix = []
            in_motif = True
            continue

        if not in_motif:
            continue

        tokens = stripped.split()
        try:
            row = [float(t) for t in tokens]
        except ValueError:
            continue
        if not row:
            continue
        if len(row) != expected_cols:
            message = (
                "Warning: skipping malformed row "
                f"(expected {expected_cols} cols, got {len(row)}): {stripped}\n"
            )
            if strict:
                raise ValueError(_without_warning_prefix(message).rstrip())
            sys.stderr.write(message)
            continue
        if is_logodds_row(row, input_format, strict=strict):
            row = logodds_to_prob(row, pseudocount, background)
        elif strict:
            _validate_probability_row(row, f"HOMER probability row for motif '{motif_id}'")
        matrix.append(row)

    if in_motif and matrix:
        yield Motif(motif_id, description, matrix, alphabet, threshold=threshold)


def read_json(
    fh: TextIO,
    pseudocount: float = 0.01,
    input_format: str = "auto",
    background: Background = 0.25,
    strict: bool = False,
) -> Iterator[Motif]:
    """Parse JSON format from a file-like object and yield Motif instances."""
    data = json.load(fh)
    for motif in data.get("motifs", []):
        mid = motif.get("id", "motif")
        desc = motif.get("description", mid)
        matrix = motif.get("matrix", [])
        alphabet = motif.get("alphabet", "ACGT")
        threshold = motif.get("threshold")
        if threshold is not None:
            try:
                threshold = float(threshold)
            except (TypeError, ValueError):
                threshold = None
        nsites = motif.get("nsites")
        if nsites is not None:
            try:
                nsites = int(nsites)
            except (TypeError, ValueError):
                nsites = None
        evalue = motif.get("evalue")
        if evalue is not None:
            try:
                evalue = float(evalue)
            except (TypeError, ValueError):
                evalue = None

        if not matrix:
            continue

        expected_cols = len(_alphabet_letters(alphabet))
        processed_matrix = []
        for row in matrix:
            if not isinstance(row, list):
                message = "Warning: skipping malformed matrix row (expected array)\n"
                if strict:
                    raise ValueError(_without_warning_prefix(message).rstrip())
                sys.stderr.write(message)
                continue
            try:
                row_values = [float(v) for v in row]
            except (TypeError, ValueError):
                message = "Warning: skipping malformed matrix row (expected numeric values)\n"
                if strict:
                    raise ValueError(_without_warning_prefix(message).rstrip())
                sys.stderr.write(message)
                continue
            if len(row_values) != expected_cols:
                message = (
                    "Warning: skipping malformed matrix row "
                    f"(expected {expected_cols} cols, got {len(row_values)})\n"
                )
                if strict:
                    raise ValueError(_without_warning_prefix(message).rstrip())
                sys.stderr.write(message)
                continue
            if is_logodds_row(row_values, input_format, strict=strict):
                row_values = logodds_to_prob(row_values, pseudocount, background)
            elif strict:
                _validate_probability_row(row_values, f"JSON probability row for motif '{mid}'")
            processed_matrix.append(row_values)

        if processed_matrix:
            yield Motif(
                mid,
                desc,
                processed_matrix,
                alphabet,
                threshold=threshold,
                nsites=nsites,
                evalue=evalue,
            )


def write_homer(
    motifs: Iterable[Motif],
    fh: TextIO,
    background: Background = 0.25,
    threshold_offset: float = 4.0,
    keep_threshold: bool = False,
    renormalize: bool = False,
) -> None:
    """Write an iterable of Motif objects to a file-like object in HOMER format."""
    for m in motifs:
        if keep_threshold and m.threshold is not None:
            score = m.threshold
        else:
            score = m.calculate_score(background, threshold_offset, renormalize=renormalize)
            if score == 0.0:
                sys.stderr.write(
                    f"Warning: HOMER threshold for motif '{m.id}' was clipped to 0 "
                    f"(t_offset={threshold_offset}); scanning may be very permissive\n"
                )
        fh.write(f">{m.id}\t{m.description}\t{score:.6f}\t0\t0\t0\n")
        for row in m.matrix:
            values = _renormalized(row) if renormalize else row
            fh.write("\t".join(f"{v:.6f}" for v in values) + "\n")


def write_meme(
    motifs: Iterable[Motif],
    fh: TextIO,
    nsites: Optional[int] = None,
    evalue: Optional[float] = None,
    renormalize: bool = False,
    background: Optional[Background] = None,
    strict: bool = False,
) -> None:
    """Write an iterable of Motif objects to a file-like object in MEME format."""
    header_printed = False
    header_alphabet = ""
    for m in motifs:
        if not header_printed:
            header_alphabet = m.alphabet
            fh.write("MEME version 4\n\n")
            fh.write(f"ALPHABET= {m.alphabet}\n\n")
            if m.alphabet in ("ACGT", "ACGU"):
                fh.write("strands: + -\n\n")
            fh.write("Background letter frequencies\n")
            fh.write(_format_background_line(m.alphabet, background) + "\n\n")
            header_printed = True
        elif m.alphabet != header_alphabet:
            sys.stderr.write(
                f"Warning: skipping motif '{m.id}' with alphabet {m.alphabet} "
                f"(header uses {header_alphabet})\n"
            )
            continue

        width = len(m.matrix)
        expected_cols = len(_alphabet_letters(m.alphabet))
        nsites_value = nsites if nsites is not None else (m.nsites if m.nsites is not None else 20)
        motif_evalue = evalue if evalue is not None else m.evalue
        evalue_value = "0" if motif_evalue is None else f"{motif_evalue:.6f}"
        fh.write(f"MOTIF {m.id} {m.description}\n\n")
        fh.write(
            "letter-probability matrix: "
            f"alength= {expected_cols} w= {width} nsites= {nsites_value} E= {evalue_value}\n"
        )
        for row in m.matrix:
            values = _renormalized(row) if renormalize else row
            if strict:
                _validate_probability_row(values, f"MEME output row for motif '{m.id}'")
            fh.write("  " + "  ".join(f"{v:.6f}" for v in values) + "\n")
        fh.write("\n")


def write_json(motifs: Iterable[Motif], fh: TextIO) -> None:
    """Write an iterable of Motif objects to a file-like object in JSON format."""
    motifs_list = list(motifs)
    fh.write("{\n")
    fh.write('  "version": "1.0",\n')
    fh.write('  "format": "motif-bridge-json",\n')
    fh.write('  "motifs": [\n')
    for mi, motif in enumerate(motifs_list):
        fh.write("    {\n")
        fh.write(f'      "id": {_json_string(motif.id)},\n')
        fh.write(f'      "description": {_json_string(motif.description)},\n')
        if motif.alphabet:
            fh.write(f'      "alphabet": {_json_string(motif.alphabet)},\n')
        if motif.threshold is not None:
            fh.write(f'      "threshold": {motif.threshold:.6f},\n')
        if motif.nsites is not None:
            fh.write(f'      "nsites": {motif.nsites},\n')
        if motif.evalue is not None:
            fh.write(f'      "evalue": {motif.evalue:.6f},\n')
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
