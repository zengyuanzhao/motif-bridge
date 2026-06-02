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
import sys

from motif_bridge import __version__
from motif_bridge.io import read_meme, write_homer, write_json


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert MEME format to HOMER motif format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s -i raw/motifs.meme -j JASPAR2026 > results/motifs.homer\n"
            "  %(prog)s -i raw/motifs.meme.gz -e MA0021.1\n"
            "  %(prog)s -i motifs.meme -b 0.25 -t 6\n"
            "  %(prog)s -i motifs.meme -b 0.29,0.21,0.21,0.29\n"
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
        metavar="<float[,float...]>",
        type=background_prob,
        default=0.25,
        help="Background probability, either scalar or comma-separated vector (default: 0.25)",
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
        help="Alphabet override for MEME input: ACGT, ACGU, or PROTEIN",
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
        "--renormalize",
        action="store_true",
        help="Renormalize each row before writing HOMER output.",
    )
    parser.add_argument(
        "--keep-threshold",
        action="store_true",
        help="Keep an existing Motif threshold when present instead of recalculating it.",
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
        original_name = m.description
        if args.e and m.id != args.e and original_name != args.e:
            continue

        m.description = f"{args.k or original_name}/{args.j}"

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
        raw_motifs = read_meme(fh, alphabet_override=args.alphabet)
        processed_motifs = process_motifs(raw_motifs, args)

        if args.format == "json":
            write_json(processed_motifs, sys.stdout)
        else:
            write_homer(
                processed_motifs,
                sys.stdout,
                background=args.b,
                threshold_offset=args.t,
                keep_threshold=args.keep_threshold,
                renormalize=args.renormalize,
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
