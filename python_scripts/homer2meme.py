#!/usr/bin/env python3
"""
homer2meme.py - Convert HOMER motif format to MEME motif format

Usage:
    python homer2meme.py -i motifs.homer > motifs.meme
    python homer2meme.py -i motifs.homer.gz > motifs.meme
    python homer2meme.py -i motifs.json -f json > motifs.meme
    python homer2meme.py -i motifs.homer --input-format logodds > motifs.meme
    python homer2meme.py -i motifs.homer --input-format logodds -b 0.2 > motifs.meme
    python homer2meme.py -i motifs.homer --alphabet ACGU > motifs.meme
    cat motifs.homer | python homer2meme.py -i -

Requires: Python 3.8+, no external dependencies
"""

import argparse
import gzip
import sys

from motif_bridge import __version__
from motif_bridge.io import read_homer, read_json, write_meme


def positive_float(value: str) -> float:
    try:
        v = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid float value: {value}") from exc
    if v <= 0:
        raise argparse.ArgumentTypeError("-a must be > 0.")
    return v


def background_prob(value: str):
    parts = value.split(",")
    values = []
    for part in parts:
        try:
            v = float(part)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"Invalid -b value: {value}") from exc
        if not (0 < v <= 1):
            raise argparse.ArgumentTypeError("-b values must be in (0, 1].")
        values.append(v)
    if len(values) > 1 and abs(sum(values) - 1.0) > 1e-3:
        raise argparse.ArgumentTypeError("-b vector must sum to 1.0.")
    return values if len(values) > 1 else values[0]


def positive_int(value: str) -> int:
    try:
        v = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid integer value: {value}") from exc
    if v <= 0:
        raise argparse.ArgumentTypeError("--nsites must be > 0.")
    return v


def nonnegative_float(value: str) -> float:
    try:
        v = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid float value: {value}") from exc
    if v < 0:
        raise argparse.ArgumentTypeError("--evalue must be >= 0.")
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
            "  %(prog)s -i motifs.homer --input-format logodds -b 0.2\n"
            "  %(prog)s -i motifs.homer --input-format logodds -b 0.29,0.21,0.21,0.29\n"
            "  %(prog)s -i motifs.homer --alphabet ACGU\n"
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
        "-b",
        "--background",
        metavar="<float[,float...]>",
        type=background_prob,
        default=0.25,
        help="Background probability, either scalar or comma-separated vector (default: 0.25)",
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
        "--alphabet",
        choices=["ACGT", "ACGU", "PROTEIN"],
        default="ACGT",
        help="Alphabet for HOMER input: ACGT (default), ACGU, or PROTEIN",
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
    parser.add_argument(
        "--nsites",
        metavar="<int>",
        type=positive_int,
        default=None,
        help="Override MEME nsites metadata in output.",
    )
    parser.add_argument(
        "--evalue",
        metavar="<float>",
        type=nonnegative_float,
        default=None,
        help="Override MEME E metadata in output.",
    )
    parser.add_argument(
        "--renormalize",
        action="store_true",
        help="Renormalize each row before writing MEME output.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Fail on malformed matrix rows, ambiguous auto-detection, or invalid probability rows."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser.parse_args()


def open_input(path: str):
    """Open plain, gzip, or stdin input."""
    if path == "-":
        return sys.stdin
    elif path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    else:
        return open(path, "r", encoding="utf-8")


def process_motifs(motifs, args):
    """Filter and apply operations to motifs."""
    for m in motifs:
        if args.e and m.id != args.e and m.description != args.e:
            continue

        if args.rc:
            try:
                m.reverse_complement()
            except ValueError as exc:
                sys.stderr.write(f"Warning: skipping motif '{m.id}': {exc}\n")
                continue
        if args.trim_edges > 0:
            m.trim_edges(args.trim_edges)
        if not m.matrix:
            continue
        if args.min_ic > 0 and m.total_ic() < args.min_ic:
            continue

        yield m


def main() -> None:
    args = parse_args()
    fh = None
    try:
        fh = open_input(args.i)

        if args.format == "json":
            raw_motifs = read_json(
                fh,
                pseudocount=args.a,
                input_format=args.input_format,
                background=args.background,
                strict=args.strict,
            )
        else:
            raw_motifs = read_homer(
                fh,
                pseudocount=args.a,
                input_format=args.input_format,
                alphabet=args.alphabet,
                background=args.background,
                strict=args.strict,
            )

        processed_motifs = process_motifs(raw_motifs, args)
        write_meme(
            processed_motifs,
            sys.stdout,
            nsites=args.nsites,
            evalue=args.evalue,
            renormalize=args.renormalize,
            background=args.background,
            strict=args.strict,
        )

    except FileNotFoundError:
        sys.exit(f"Error: Cannot open file: {args.i}")
    except ValueError as exc:
        sys.exit(f"Error: {exc}")
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
