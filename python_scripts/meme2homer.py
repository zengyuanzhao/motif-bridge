#!/usr/bin/env python3
"""
meme2homer.py - Convert MEME motif format to HOMER motif format

Usage:
    python meme2homer.py -i motifs.meme -j JASPAR2026 > motifs.homer
    python meme2homer.py -i motifs.meme.gz -e MA0021.1
    python meme2homer.py -i motifs.meme -f json > motifs.json
    python meme2homer.py -i motifs.meme --alphabet ACGU > motifs.homer
    cat motifs.meme | python meme2homer.py -i -

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

ALPHABETS = {
    "ACGT": "ACGT",
    "ACGU": "ACGU",
    "PROTEIN": "ACDEFGHIKLMNPQRSTVWY",
}


def bg_prob(value: str) -> float:
    try:
        v = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid float value: {value}") from exc
    if not (0 < v <= 1):
        raise argparse.ArgumentTypeError("-b must be in (0, 1].")
    return v


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert MEME format to HOMER motif format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s -i raw/motifs.meme -j JASPAR2026 > results/motifs.homer\n"
            "  %(prog)s -i raw/motifs.meme.gz -e MA0021.1\n"
            "  %(prog)s -i motifs.meme -b 0.25 -t 6\n"
            "  %(prog)s -i motifs.meme -f json > motifs.json\n"
            "  %(prog)s -i motifs.meme --alphabet ACGU\n"
            "  cat motifs.meme | %(prog)s -i -\n"
        ),
    )
    parser.add_argument(
        "-i",
        metavar="<file>",
        required=True,
        help="Input MEME format file (or '-' for stdin, supports .gz)",
    )
    parser.add_argument("-j", metavar="<string>", default="NA", help="Database name (default: NA)")
    parser.add_argument(
        "-k",
        metavar="<string>",
        default="",
        help="Override motif name (default: name from MEME file)",
    )
    parser.add_argument(
        "-e", metavar="<string>", default="", help="Extract only specified motif by id or name"
    )
    parser.add_argument(
        "-b",
        metavar="<float>",
        type=bg_prob,
        default=0.25,
        help="Background nucleotide probability (default: 0.25)",
    )
    parser.add_argument(
        "-t",
        metavar="<float>",
        type=float,
        default=4.0,
        help="Threshold offset in log2 bits (default: 4.0)",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=["homer", "json"],
        default="homer",
        help="Output format: homer (default) or json",
    )
    parser.add_argument(
        "--alphabet",
        choices=["ACGT", "ACGU", "PROTEIN"],
        default=None,
        help="Alphabet: ACGT (DNA, default), ACGU (RNA), or PROTEIN",
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


def print_motif_homer(
    motif_id: str,
    description: str,
    matrix: List[List[float]],
    background: float,
    threshold_offset: float,
) -> None:
    """Print a single motif in HOMER format (6 tab-separated header fields)."""
    m = Motif(motif_id, description, matrix)
    score = m.calculate_score(background, threshold_offset)
    print(f">{motif_id}\t{description}\t{score:.6f}\t0\t0\t0")
    for row in matrix:
        print("\t".join(f"{v:.6f}" for v in row))


def process_and_output_motif(
    m: Motif,
    background: float,
    threshold_offset: float,
    output_format: str,
    do_rc: bool,
    trim_edges: float,
    min_ic: float,
    motifs_list: List[dict],
) -> None:
    if do_rc:
        m.reverse_complement()
    if trim_edges > 0:
        m.trim_edges(trim_edges)
    if not m.matrix:
        return
    if min_ic > 0 and m.total_ic() < min_ic:
        return

    if output_format == "json":
        motifs_list.append(m.to_dict())
    else:
        print_motif_homer(m.id, m.description, m.matrix, background, threshold_offset)


def parse_and_convert(
    fh,
    db: str,
    motif_name: str,
    extract: str,
    background: float,
    threshold_offset: float,
    output_format: str = "homer",
    alphabet: str = "ACGT",
    do_rc: bool = False,
    trim_edges: float = 0.0,
    min_ic: float = 0.0,
) -> List[dict]:
    in_motif = False
    in_matrix = False
    motif_id = ""
    description = ""
    matrix: List[List[float]] = []
    motifs: List[dict] = []
    expected_cols = len(ALPHABETS.get(alphabet, alphabet))

    for raw_line in fh:
        line = raw_line.rstrip("\n")
        stripped = line.strip()

        if stripped.startswith("MOTIF"):
            if in_motif and matrix:
                m = Motif(motif_id, description, matrix, alphabet)
                process_and_output_motif(
                    m,
                    background,
                    threshold_offset,
                    output_format,
                    do_rc,
                    trim_edges,
                    min_ic,
                    motifs,
                )

            parts = stripped.split()
            mid = parts[1] if len(parts) > 1 else ""
            if not mid:
                sys.stderr.write(f"Warning: skipping malformed MOTIF line without ID: {stripped}\n")
                in_motif = False
                in_matrix = False
                matrix = []
                continue
            original_name = " ".join(parts[2:]) if len(parts) > 2 else (mid or "motif")

            if extract and mid != extract and original_name != extract:
                in_motif = False
                in_matrix = False
                matrix = []
                continue

            motif_id = mid
            description = f"{motif_name or original_name}/{db}"
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
                m = Motif(motif_id, description, matrix, alphabet)
                process_and_output_motif(
                    m,
                    background,
                    threshold_offset,
                    output_format,
                    do_rc,
                    trim_edges,
                    min_ic,
                    motifs,
                )
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

        if in_matrix and stripped and (stripped[0].isdigit() or stripped.startswith(".")):
            tokens = stripped.split()
            try:
                row = [float(t) for t in tokens]
                if len(row) == expected_cols:
                    matrix.append(row)
                elif row:
                    sys.stderr.write(
                        f"Warning: skipping malformed matrix row "
                        f"(expected {expected_cols} cols, got {len(row)}): {stripped}\n"
                    )
            except ValueError:
                pass

    if in_motif and matrix:
        m = Motif(motif_id, description, matrix, alphabet)
        process_and_output_motif(
            m, background, threshold_offset, output_format, do_rc, trim_edges, min_ic, motifs
        )

    return motifs


def main() -> None:
    args = parse_args()
    fh = None
    try:
        fh = open_input(args.i)
        motifs = parse_and_convert(
            fh,
            db=args.j,
            motif_name=args.k,
            extract=args.e,
            background=args.b,
            threshold_offset=args.t,
            output_format=args.format,
            alphabet=args.alphabet or "ACGT",
            do_rc=args.rc,
            trim_edges=args.trim_edges,
            min_ic=args.min_ic,
        )
        if args.format == "json":
            json.dump({"version": "1.0", "source": "meme", "motifs": motifs}, sys.stdout, indent=2)
            print()
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
