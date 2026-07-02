# motif-bridge

[![Bioconda](https://img.shields.io/conda/vn/bioconda/motif-bridge.svg)](https://anaconda.org/bioconda/motif-bridge)
[![Crates.io](https://img.shields.io/crates/v/motif-bridge.svg)](https://crates.io/crates/motif-bridge)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21136623.svg)](https://doi.org/10.5281/zenodo.21136623)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A cross-language toolkit for bridging the **MEME** and **HOMER** motif analysis ecosystems.

MEME and HOMER are widely used platforms for motif analysis, but their motif file formats are not directly interoperable. This project provides bidirectional motif matrix conversion tools implemented in **Perl**, **Python**, and **Rust**, allowing users to move standard motif matrices between the two ecosystems.

The project is also available as a reusable library:
- **Python SDK**: `from motif_bridge.core import Motif`
- **Rust Crate**: `motif_bridge::Motif`

The converters preserve motif matrix content for standard A/C/G/T motifs after numeric formatting. RNA and protein alphabets are supported for parsing, writing, JSON I/O, IC calculation, and filtering operations. Some format-specific metadata, such as HOMER log-p, pseudo-counts, sites, and MEME `nsites` or `E` values, may be regenerated or set to documented defaults during conversion.

---

## Installation

> **Published package status**: Crates.io carries the current Rust release. The Bioconda package installs the Rust binaries from the crates.io source package and may lag while recipe updates are reviewed. For unreleased `main` changes, use the Git-based Cargo install or build from source.

### Bioconda (published releases)

```bash
conda install -c bioconda motif-bridge
```

This installs the precompiled Rust binaries (`meme2homer` and `homer2meme`) directly to your PATH. No compiler or manual build steps required.

The Bioconda recipe packages the Rust implementation only. The Python package and Perl scripts remain available from the source repository.

### Cargo (Rust)

```bash
# Published Rust crate
cargo install motif-bridge

# Latest source version from main
cargo install --git https://github.com/zengyuanzhao/motif-bridge
```

### From source

```bash
git clone https://github.com/zengyuanzhao/motif-bridge
cd motif-bridge/rust_scripts
cargo build --release
# Produces: target/release/meme2homer
#           target/release/homer2meme
```

The Perl scripts can be run directly. The Python package can be installed in editable mode:

```bash
pip install -e .
meme2homer-py -i motifs.meme -j JASPAR2026 > motifs.homer
homer2meme-py -i motifs.homer > motifs.meme
```

Alternatively, set `PYTHONPATH` to the repository root and run the scripts directly.

**Requirements:** Python 3.8+, Perl 5 with `IO::Uncompress::Gunzip`, or Rust 1.70+.

---

## Citation

If you use motif-bridge in a workflow or publication, cite the software using [`CITATION.cff`](CITATION.cff). The Zenodo concept DOI for all versions is [`10.5281/zenodo.21136623`](https://doi.org/10.5281/zenodo.21136623); the archived `v0.3.1` release DOI is [`10.5281/zenodo.21136624`](https://doi.org/10.5281/zenodo.21136624).

---

## Why This Exists

| | MEME Suite | HOMER |
|---|---|---|
| **Motif format** | `.meme`, probability matrix with header block | `.motif`, log-odds or probability, tab-separated |
| **Threshold** | E-value / p-value based | Per-motif log-odds score |
| **Typical use** | *de novo* motif discovery, FIMO scanning | ChIP-seq peak annotation, known motif enrichment |
| **Downstream tools** | FIMO, MAST, TOMTOM | `findMotifsGenome.pl`, `annotatePeaks.pl` |

Converting between the two lets you:
- Use MEME-discovered motifs directly in HOMER annotation pipelines
- Scan HOMER motif databases with FIMO/MAST
- Build unified motif databases, such as JASPAR-derived databases for both platforms

---

## Repository Layout

```
motif-bridge/
├── README.md
├── LICENSE
├── motif_bridge/          # Reusable Python library
│   ├── __init__.py
│   ├── core.py
│   └── io.py
├── perl_scripts/          # Perl 5, runs on most servers without compilation
│   ├── meme2homer.pl
│   └── homer2meme.pl
├── python_scripts/        # Python 3.8+ CLI wrappers
│   ├── meme2homer.py
│   └── homer2meme.py
└── rust_scripts/          # Rust implementation, also available as Crate
    ├── Cargo.toml
    ├── Cargo.lock
    └── src/
        ├── lib.rs         # Reusable Rust library (motif_bridge crate)
        └── bin/
            ├── meme2homer.rs
            └── homer2meme.rs
```

All three language implementations share the same main CLI flags and are intended to produce matching output for valid standard MEME/HOMER inputs. Choose the implementation that fits your environment.

---

## Quick Start

### MEME → HOMER

```bash
# Bioconda / Cargo install
meme2homer -i motifs.meme -j JASPAR2026 > motifs.homer

# Perl
perl perl_scripts/meme2homer.pl -i motifs.meme -j JASPAR2026 > motifs.homer

# Python script
python3 python_scripts/meme2homer.py -i motifs.meme -j JASPAR2026 > motifs.homer

# Python package entry point
meme2homer-py -i motifs.meme -j JASPAR2026 > motifs.homer

# Rust, build from source
cd rust_scripts && cargo build --release && cd ..
./rust_scripts/target/release/meme2homer -i motifs.meme -j JASPAR2026 > motifs.homer
```

### HOMER → MEME

```bash
# Bioconda / Cargo install
homer2meme -i motifs.homer > motifs.meme

# Perl
perl perl_scripts/homer2meme.pl -i motifs.homer > motifs.meme

# Python script
python3 python_scripts/homer2meme.py -i motifs.homer > motifs.meme

# Python package entry point
homer2meme-py -i motifs.homer > motifs.meme

# Rust, build from source
cd rust_scripts && cargo build --release && cd ..
./rust_scripts/target/release/homer2meme -i motifs.homer > motifs.meme
```

---

## Shared CLI Flags

### meme2homer (all implementations)

| Flag | Description | Default |
|---|---|---|
| `-i <file>` | Input MEME file (`-` for stdin, `.gz` supported) | *(required)* |
| `-j <string>` | Database name appended to motif description | `NA` |
| `-k <string>` | Override motif name from file | *(from file)* |
| `-e <string>` | Extract only the specified motif by ID or name | *(all motifs)* |
| `-b <float[,float...]>` | Background probability used for HOMER threshold calculation; scalar or per-column vector matching matrix width, values in `(0,1]`, vector sum within `1e-3` of `1.0` | `0.25` |
| `-t <float>` | Threshold offset subtracted from log-odds score (log2 bits) | `4.0` |
| `-f, --format <fmt>` | Output format: `homer` or `json` | `homer` |
| `--alphabet <string>` | Alphabet override (`ACGT`, `ACGU`, or `PROTEIN`) | *(auto from MEME header)* |
| `--rc` | Output the reverse complement of the motif (DNA/RNA only) | *(off)* |
| `--trim-edges <float>` | Trim edges with Information Content below threshold | `0.0` |
| `--min-ic <float>` | Filter out motifs with total Information Content below threshold | `0.0` |
| `--renormalize` | Renormalize rows before writing HOMER output | *(off)* |
| `--keep-threshold` | Keep an existing motif threshold when present instead of recalculating; plain MEME CLI input has no threshold metadata, so this matters only for library/JSON-derived motif objects | *(off)* |
| `--strict` | Fail on malformed matrix rows, values outside `[0, 1]`, row sums outside tolerance, or `alength=`/alphabet conflicts | *(off)* |
| `--version` | Show version and exit | |
| `-h` | Show help | |

### homer2meme (all implementations)

| Flag | Description | Default |
|---|---|---|
| `-i <file>` | Input HOMER motif file (`-` for stdin, `.gz` supported, `.json` supported) | *(required)* |
| `-e <string>` | Extract only the specified motif by ID or description | *(all motifs)* |
| `-a <float>` | Pseudocount for log-odds to probability conversion | `0.01` |
| `-b <float[,float...]>` | Background probability for log-odds conversion; scalar or per-column vector matching matrix width, values in `(0,1]`, vector sum within `1e-3` of `1.0` | `0.25` |
| `-f, --format <fmt>` | Input format: `homer` or `json` | `homer` |
| `--input-format <fmt>` | Matrix type: `auto`, `logodds`, or `probability` | `auto` |
| `--alphabet <string>` | Alphabet (`ACGT`, `ACGU`, or `PROTEIN`) | `ACGT` |
| `--rc` | Output the reverse complement of the motif (DNA/RNA only) | *(off)* |
| `--trim-edges <float>` | Trim edges with Information Content below threshold | `0.0` |
| `--min-ic <float>` | Filter out motifs with total Information Content below threshold | `0.0` |
| `--nsites <int>` | Override MEME `nsites` metadata in output | `20` |
| `--evalue <float>` | Override MEME `E` metadata in output | `0` |
| `--renormalize` | Renormalize rows before writing MEME output | *(off)* |
| `--strict` | Fail on malformed rows, ambiguous `--input-format auto` rows, or invalid probability rows | *(off)* |
| `--version` | Show version and exit | |
| `-h` | Show help | |

Note: `homer2meme` auto-detects log-odds vs probability rows by checking whether row sum is near 1.0 (`[0.98, 1.02]`). This is a practical heuristic and may be ambiguous for edge-case inputs whose rows are nonnegative and sum close to the boundary. When the source format is known, prefer `--input-format logodds` or `--input-format probability` to avoid misclassification.

---

## Alphabet Behavior

`meme2homer` reads the MEME `ALPHABET=` line automatically when `--alphabet` is omitted. Use `--alphabet ACGT`, `--alphabet ACGU`, or `--alphabet PROTEIN` to override the detected alphabet.

`homer2meme` uses `--alphabet` because HOMER motif files do not carry a global alphabet declaration. For RNA output, background frequencies use `U` instead of `T`. For protein output, the MEME writer omits the `strands: + -` line.

Reverse complement is only defined for DNA and RNA motifs (`ACGT` and `ACGU`). Protein motifs are skipped with a warning when `--rc` is requested.

---

## Format Reference

### MEME format (input/output)

```
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF MA0021.1 CTCF

letter-probability matrix: alength= 4 w= 19 nsites= 4000 E= 0
  0.347475  0.154857  0.163364  0.334304
  ...
```

For `ACGU`, background frequencies use `U`. For `PROTEIN`, the `strands: + -` line is omitted.

### HOMER format (input/output)

```
>MA0021.1  CTCF/JASPAR2026  10.234560  0  0  0
0.347475  0.154857  0.163364  0.334304
...
```

HOMER header fields: `>id \t description \t threshold \t log-p \t pseudo \t sites`

The **threshold score** is computed as:

```
threshold = sum_over_positions( log2(max_prob / background) ) - t_offset
threshold = max(threshold, 0)
```

---

## JSON I/O

Both converters support JSON as an intermediate exchange format:

```bash
meme2homer -i motifs.meme -f json > motifs.json
homer2meme -i motifs.json -f json > motifs.meme
```

The JSON document uses top-level `"format": "motif-bridge-json"` and `"version": "1.0"` markers. Each motif stores `id`, `description`, `alphabet`, and `matrix`. Optional metadata fields `threshold`, `nsites`, and `evalue` are included when they are present on the motif object, so library workflows such as `read_homer` -> `write_json` -> `read_json` -> `write_homer(..., keep_threshold=True)` can preserve a parsed HOMER threshold. Unicode motif descriptions are preserved.

For JSON input containing mixed alphabets, `homer2meme` writes one MEME file using the alphabet of the first emitted motif. Later motifs with incompatible alphabets are skipped with warnings.

---

## Motif Operations

### Reverse Complement (`--rc`)

Reverses the matrix order and swaps A↔T, C↔G columns for DNA motifs, or A↔U, C↔G columns for RNA motifs. Automatically appends `_RC` to the motif ID.

```bash
meme2homer -i motifs.meme --rc > motifs_rc.homer
```

### Edge Trimming (`--trim-edges`)

Calculates Information Content (IC) per position and trims non-conserved edges below the threshold (in bits).

```bash
meme2homer -i motifs.meme --trim-edges 0.5 > motifs_trimmed.homer
```

### Information Content Filter (`--min-ic`)

Filters out low-quality motifs whose total IC falls below the threshold.

```bash
meme2homer -i motifs.meme --min-ic 5.0 > motifs_filtered.homer
```

---

## SDK / Library Usage

### Python SDK

```python
import sys
from motif_bridge import Motif, read_meme, write_homer

# Read motifs from a MEME file, yields a generator to save memory
with open("motifs.meme", "r") as f:
    motifs = list(read_meme(f))

for motif in motifs:
    # Filter by Information Content
    if motif.total_ic() > 5.0:
        # Trim low-IC edges
        motif.trim_edges(threshold=0.5)
        # Reverse complement, DNA/RNA only
        motif.reverse_complement()

# Write processed motifs to stdout in HOMER format
write_homer(motifs, sys.stdout, background=0.25)
```

### Rust Crate

Add to your `Cargo.toml` (using git dependency for the latest source version):

```toml
[dependencies]
motif-bridge = { git = "https://github.com/zengyuanzhao/motif-bridge" }
```

```rust
use motif_bridge::io::{read_meme, write_homer};
use std::fs::File;
use std::io::{BufReader, BufWriter};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open("motifs.meme")?;
    let reader = BufReader::new(file);

    // Read motifs from a MEME file
    let mut motifs = read_meme(reader, None)?;

    for motif in motifs.iter_mut() {
        if motif.total_ic() > 5.0 {
            motif.trim_edges(0.5);
            motif.reverse_complement()?;
        }
    }

    // Write to a new HOMER file
    let out_file = File::create("processed.homer")?;
    let mut writer = BufWriter::new(out_file);
    write_homer(&mut writer, &motifs, &[0.25], 4.0, false, false)?;

    Ok(())
}
```

---

## Metadata and Round-trip Limits

The current converters focus on matrix-level conversion. During conversion, some metadata fields are written as defaults:

| Direction | Metadata behavior |
|---|---|
| MEME to HOMER | HOMER `log-p`, `pseudo`, and `sites` are written as `0`; threshold is recalculated from the matrix |
| HOMER to MEME | Plain HOMER text has no MEME metadata, so `nsites` is written as `20` and `E` as `0` unless overridden with `--nsites` / `--evalue`; JSON input can carry `nsites` / `evalue` metadata |
| Both directions | Matrix values are printed to six decimal places |

For this reason, round-trip conversion should be interpreted as matrix-level consistency rather than byte-for-byte recovery of the original file.

### Downstream Interpretation Caveats

The converted files are suitable as matrix-level exchange artifacts, but they should not be treated as statistically identical to the source files in motif-scanning workflows.

- HOMER thresholds are recalculated during `meme2homer` output as `sum(log2(max_prob / background)) - t_offset`, then clipped to be non-negative. A warning is emitted when the calculated threshold clips to `0`, because HOMER scans may become very permissive. Original HOMER detection thresholds are not preserved through plain `homer2meme` -> `meme2homer` text round trips; JSON can carry `threshold`, and library callers can retain a parsed threshold with `write_homer(..., keep_threshold=True)`.
- `homer2meme --input-format auto` classifies rows as probabilities only when the row sum is near 1.0 (`[0.98, 1.02]`). This is practical for standard HOMER known motifs, but rounded, counted, or externally generated matrices can be misclassified. Nonnegative rows close to the auto boundary emit a warning. For reproducible analyses, pass `--input-format probability` or `--input-format logodds` explicitly.
- MEME `nsites` and `E` metadata are regenerated as `nsites= 20` and `E= 0` for plain HOMER input unless overridden. MEME and JSON readers preserve `nsites` / `evalue` metadata when available. Tools such as FIMO or TOMTOM can use `nsites` in pseudocount or statistical calculations, so pass `--nsites` / `--evalue` when source-specific values matter downstream.
- Background defaults are uniform (`0.25`). The `-b` option can be either a scalar or a comma-separated per-column vector, for example `-b 0.29,0.21,0.21,0.29`, used by threshold and log-odds conversions. Valid vectors are also written to the MEME `Background letter frequencies` header during `homer2meme`; scalar values other than the uniform default remain conversion parameters and are not written as normalized MEME background vectors. Retune or validate thresholds when scanner background assumptions matter.
- With `--rc`, motif columns are reverse-complemented back into the same alphabet order, but a non-symmetric `-b` vector is not reverse-complemented. This does not matter for symmetric vectors such as `0.29,0.21,0.21,0.29`, but asymmetric vectors can shift recalculated thresholds slightly.
- Matrix rows are printed to six decimal places and are not renormalized by default. Most motif tools tolerate tiny row-sum drift, but consumers can use `--renormalize` to scale each row before writing or `--strict` to fail on invalid probability rows instead of warning/skipping.

---

## Validation

The project was validated on a server using the following real motif datasets:

| File | Description |
|---|---|
| `homer.known.motifs` | HOMER built-in known motif library (436 motifs) |
| `JASPAR2024_small.meme` | First 200 lines of `JASPAR2024_vertebrates.meme` (12 motifs), used for rapid iteration |
| `JASPAR2024_vertebrates.meme` | Full JASPAR 2024 vertebrate motif database (879 motifs, 469 KB) |

### Test coverage

Run `bash test_motif_bridge.sh` locally. The test suite covers:

| Test stage | Description |
|---|---|
| 1. Cross-language consistency | Diff between Perl, Python, and Rust outputs |
| 2. Cross-language consistency (homer2meme) | Diff between Perl, Python, and Rust outputs |
| 3. Single motif extraction | Test `-e` flag for extracting a named motif |
| 4. stdin pipeline | Test `cat file | tool -i -` |
| 5. gzip input | Test automatic decompression of `.gz` inputs |
| 6. Round-trip | Validate matrix-level consistency for `meme→homer→meme` and `homer→meme→homer` |
| 7. Log-odds conversion | Verify log-odds to probability conversion |
| 8. Format compliance | Validate MEME headers, HOMER column counts, probability row sums |
| 9. JSON I/O | Validate JSON output, Unicode descriptions, JSON input, round-trip, and mixed-alphabet handling |
| 10. Explicit input format | Test `--input-format` flag (`auto`, `logodds`, `probability`) |
| 11. Alphabet support | Test `--alphabet`, MEME `ALPHABET=` auto-detection, RNA headers, and backgrounds |
| 12. Motif operations | Test `--rc`, `--trim-edges`, and `--min-ic` flags |
| 13. MOTIF word boundary | Ensure `MOTIF` lines require a word boundary (`MOTIFY` ignored) |
| 14. Negative matrix warnings | Ensure negative MEME values trigger warnings |
| 15. Version and parser metadata | Test `--version`, Perl version parsing, and `ALPHABET=`/`alength=` conflicts |
| 16. Safety flags and metadata | Test threshold/auto warnings, metadata overrides, renormalization, strict-mode failures, background-vector validation/header output, tied-max scoring parity, and JSON metadata preservation |

### Continuous Integration

GitHub Actions currently checks:

| Area | Coverage |
|---|---|
| Python versions | Python 3.8, 3.10, and 3.12 |
| CLI regression tests | Full `test_motif_bridge.sh` run |
| Python packaging | Editable install, wheel build, wheel install, and entry point checks |
| Python style | `ruff check` and `ruff format --check` on `python_scripts/` and `motif_bridge/` |
| Rust quality | `cargo fmt --check`, `cargo clippy -- -D warnings`, and `cargo test` |
| Perl quality | Syntax checks for both Perl scripts |

### Latest test result

The exact number of shell-level checks may vary depending on whether Rust is available in the local environment. In CI, all regression, packaging, linting, and formatting checks are expected to pass before release.

### Performance snapshot (large-file benchmark, 879 motifs)

| Implementation | Time | Relative speed |
|---|---|---|
| Perl | 125 ms | Baseline |
| Python | 100 ms | 1.25× |
| Rust | **23 ms** | **5.4×** |

---

## Choosing a Language Implementation

| Scenario | Recommended |
|---|---|
| Quickest install, conda environment | `conda install -c bioconda motif-bridge` |
| Server without compiler, Perl available | `perl_scripts/` |
| Conda / Python environment | `python_scripts/` or `meme2homer-py` / `homer2meme-py` |
| Large-scale batch processing | `rust_scripts/` or Cargo-installed binaries |
| `.gz` compressed input | All three implementations |
| Fastest runtime on tested server data | Rust implementation |

---
