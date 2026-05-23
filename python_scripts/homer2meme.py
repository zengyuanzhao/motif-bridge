#!/usr/bin/env python3
"""
homer2meme.py - Convert HOMER motif format to MEME motif format

Usage:
    python homer2meme.py -i motifs.homer > motifs.meme
    python homer2meme.py -i motifs.homer.gz > motifs.meme
    python homer2meme.py -i motifs.json -f json > motifs.meme
    python homer2meme.py -i motifs.homer --input-format logodds > motifs.meme
    cat motifs.homer | python homer2meme.py -i -

Requires: Python 3.8+, no external dependencies
"""

import argparse
import gzip
import json
import sys
from typing import List

try:
    from motif_bridge.core import Motif
except ImportError:
    import os

    # Add project root to sys.path
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if root not in sys.path:
        sys.path.insert(0, root)
    from motif_bridge.core import Motif


def positive_float(value: str) -> float:
    try:
        v = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid float value: {value}") from exc
    if v <= 0:
        raise argparse.ArgumentTypeError("-a must be > 0.")
    return v


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert HOMER motif format to MEME format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s -i results/motifs.homer > raw/motifs.meme\n"
            "  %(prog)s -i results/motifs.homer.gz > raw/motifs.meme\n"
            '  %(prog)s -i results/motifs.homer -e "CTCF/Jaspar"\n'
            "  %(prog)s -i motifs.json -f json > motifs.meme\n"
            "  %(prog)s -i motifs.homer --input-format logodds\n"
            "  cat motifs.homer | %(prog)s -i -\n"
        ),
    )
    parser.add_argument(
        "-i",
        metavar="<file>",
        required=True,
        help="Input HOMER format file (or '-' for stdin, supports .gz)",
    )
    parser.add_argument(
        "-e",
        metavar="<string>",
        default="",
        help="Extract only specified motif by id or description",
    )
    parser.add_argument(
        "-a",
        metavar="<float>",
        type=positive_float,
        default=0.01,
        help="Pseudocount for log-odds -> probability (default: 0.01)",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=["homer", "json"],
        default="homer",
        help="Input format: homer (default) or json",
    )
    parser.add_argument(
        "--input-format",
        choices=["auto", "logodds", "probability"],
        default="auto",
        help="Matrix type: auto (default), logodds, or probability",
    )
    parser.add_argument(
        "--rc",
        action="store_true",
        help="Output the reverse complement of the motif (DNA/RNA only).",
    )
    parser.add_argument(
        "--trim-edges",
        metavar="<float>",
        type=float,
        default=0.0,
        help="Trim edges with information content below threshold.",
    )
    parser.add_argument(
        "--min-ic",
        metavar="<float>",
        type=float,
        default=0.0,
        help="Filter out motifs with total information content below threshold.",
    )
    return parser.parse_args()


def open_input(path: str):
    """Open plain, gzip, or stdin input."""
    if path == "-":
        return sys.stdin
    elif path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    else:
        return open(path, "r", encoding="utf-8")


def logodds_to_prob(row: List[float], pseudocount: float = 0.01) -> List[float]:
    """
    Convert a HOMER log-odds row to a probability row.
    HOMER stores log2(p / 0.25). Reverse: p = 2^x * 0.25, then normalize.
    """
    background = 0.25
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


def print_meme_header(alphabet: str = "ACGT") -> None:
    print("MEME version 4")
    print()
    print(f"ALPHABET= {alphabet}")
    print()
    print("strands: + -")
    print()
    print("Background letter frequencies")
    print("A 0.25 C 0.25 G 0.25 T 0.25")
    print()


def process_and_output_meme(
    m: Motif, do_rc: bool, trim_edges: float, min_ic: float, header_printed: bool
) -> bool:
    if do_rc:
        m.reverse_complement()
    if trim_edges > 0:
        m.trim_edges(trim_edges)
    if not m.matrix:
        return header_printed
    if min_ic > 0 and m.total_ic() < min_ic:
        return header_printed

    if not header_printed:
        print_meme_header(m.alphabet)
        header_printed = True
    m.print_meme_motif()
    return header_printed


def parse_and_convert_homer(
    fh,
    extract: str,
    pseudocount: float,
    input_format: str,
    do_rc: bool,
    trim_edges: float,
    min_ic: float,
) -> None:
    header_printed = False

    in_motif = False
    motif_id = ""
    description = ""
    matrix: List[List[float]] = []

    for raw_line in fh:
        line = raw_line.rstrip("\n")
        stripped = line.strip()

        if not stripped:
            continue

        if stripped.startswith(">"):
            if in_motif and matrix:
                m = Motif(motif_id, description, matrix, "ACGT")
                header_printed = process_and_output_meme(
                    m, do_rc, trim_edges, min_ic, header_printed
                )

            parts = stripped[1:].split("\t")
            mid = parts[0] if len(parts) > 0 else "motif"
            desc = parts[1] if len(parts) > 1 else mid

            if extract and mid != extract and desc != extract:
                in_motif = False
                matrix = []
                continue

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
            if len(row) != 4:
                sys.stderr.write(
                    f"Warning: skipping malformed row "
                    f"(expected 4 cols, got {len(row)}): {stripped}\n"
                )
                continue
            if is_logodds_row(row, input_format):
                row = logodds_to_prob(row, pseudocount)
            matrix.append(row)
        except ValueError:
            pass

    if in_motif and matrix:
        m = Motif(motif_id, description, matrix, "ACGT")
        header_printed = process_and_output_meme(m, do_rc, trim_edges, min_ic, header_printed)


def parse_and_convert_json(
    fh, extract: str, pseudocount: float, do_rc: bool, trim_edges: float, min_ic: float
) -> None:
    data = json.load(fh)
    header_printed = False

    for motif in data.get("motifs", []):
        mid = motif.get("id", "motif")
        desc = motif.get("description", mid)
        matrix = motif.get("matrix", [])
        alphabet = motif.get("alphabet", "ACGT")

        if extract and mid != extract and desc != extract:
            continue

        if not matrix:
            continue

        processed_matrix = []
        for row in matrix:
            if len(row) != 4:
                sys.stderr.write(
                    f"Warning: skipping malformed matrix row (expected 4 cols, got {len(row)})\n"
                )
                continue
            if is_logodds_row(row, "auto"):
                row = logodds_to_prob(row, pseudocount)
            processed_matrix.append(row)

        if processed_matrix:
            m = Motif(mid, desc, processed_matrix, alphabet)
            header_printed = process_and_output_meme(m, do_rc, trim_edges, min_ic, header_printed)


def main() -> None:
    args = parse_args()
    fh = None
    try:
        fh = open_input(args.i)
        if args.format == "json":
            parse_and_convert_json(
                fh,
                extract=args.e,
                pseudocount=args.a,
                do_rc=args.rc,
                trim_edges=args.trim_edges,
                min_ic=args.min_ic,
            )
        else:
            parse_and_convert_homer(
                fh,
                extract=args.e,
                pseudocount=args.a,
                input_format=args.input_format,
                do_rc=args.rc,
                trim_edges=args.trim_edges,
                min_ic=args.min_ic,
            )
    except FileNotFoundError:
        sys.exit(f"Error: Cannot open file: {args.i}")
    except BrokenPipeError:
        pass
    finally:
        if fh is not None and args.i != "-":
            try:
                fh.close()
            except Exception:
                pass


if __name__ == "__main__":
    main()
