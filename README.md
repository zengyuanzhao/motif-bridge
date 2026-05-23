# motif-bridge

[![Bioconda](https://img.shields.io/conda/vn/bioconda/motif-bridge.svg)](https://anaconda.org/bioconda/motif-bridge)
[![Crates.io](https://img.shields.io/crates/v/motif-bridge.svg)](https://crates.io/crates/motif-bridge)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A cross-language toolkit for bridging the **MEME** and **HOMER** motif analysis ecosystems.

MEME and HOMER are widely used platforms for motif analysis, but their motif file formats are not directly interoperable. This project provides bidirectional motif matrix conversion tools implemented in **Perl**, **Python**, and **Rust**, allowing users to move standard DNA motif matrices between the two ecosystems.

The project is also available as a reusable library:
- **Python SDK**: `from motif_bridge.core import Motif`
- **Rust Crate**: `motif_bridge::Motif`

The converters preserve motif matrix content for standard A/C/G/T motifs after numeric formatting. Some format-specific metadata, such as HOMER log-p, pseudo-counts, sites, and MEME `nsites` or `E` values, may be regenerated or set to documented defaults during conversion.

---

## Installation

> **⚠️ Note on version 0.2.0**: The new features described in this document (Python SDK, Rust Crate library, JSON I/O, Motif operations like `--rc` and `--trim-edges`) are currently in the `main` branch. The published versions on Bioconda and Crates.io are currently at `v0.1.0`. To use the `v0.2.0` features, please install from source.

### Bioconda (v0.1.0 only)

```bash
conda install -c bioconda motif-bridge
```

This installs the precompiled Rust binaries (`meme2homer` and `homer2meme`) directly to your PATH. No compiler or manual build steps required.

### Cargo (Rust)

```bash
cargo install motif-bridge
```

### From source

```bash
git clone https://github.com/zengyuanzhao/motif-bridge
cd motif-bridge/rust_scripts
cargo build --release
# Produces: target/release/meme2homer
#           target/release/homer2meme
```

For the Perl and Python implementations, no installation is needed. The scripts can be run directly.

**Requirements:** Python 3.8+, Perl 5 with `IO::Uncompress::Gunzip`, or Rust 1.70+.

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
- Build unified motif databases (e.g., JASPAR to both platforms)

---

## Repository Layout

```
motif-bridge/
├── README.md
├── LICENSE
├── perl_scripts/          # Perl 5, runs on most servers without compilation
│   ├── meme2homer.pl
│   └── homer2meme.pl
├── python_scripts/        # Python 3.8+, also available as SDK
│   ├── motif_bridge/      # Reusable Python library
│   │   ├── __init__.py
│   │   └── core.py
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

All three language implementations share the same main CLI flags and are intended to produce matching output for valid standard DNA MEME/HOMER inputs. Choose the implementation that fits your environment.

---

## Quick Start

### MEME → HOMER

```bash
# Bioconda / Cargo install
meme2homer -i motifs.meme -j JASPAR2026 > motifs.homer

# Perl
perl perl_scripts/meme2homer.pl -i motifs.meme -j JASPAR2026 > motifs.homer

# Python
python3 python_scripts/meme2homer.py -i motifs.meme -j JASPAR2026 > motifs.homer

# Rust (build from source)
cd rust_scripts && cargo build --release && cd ..
./rust_scripts/target/release/meme2homer -i motifs.meme -j JASPAR2026 > motifs.homer
```

### HOMER → MEME

```bash
# Bioconda / Cargo install
homer2meme -i motifs.homer > motifs.meme

# Perl
perl perl_scripts/homer2meme.pl -i motifs.homer > motifs.meme

# Python
python3 python_scripts/homer2meme.py -i motifs.homer > motifs.meme

# Rust (build from source)
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
| `-b <float>` | Background nucleotide probability | `0.25` |
| `-t <float>` | Threshold offset subtracted from log-odds score (log2 bits) | `4.0` |
| `-f, --format <fmt>` | Output format: `homer` or `json` | `homer` |
| `--alphabet <string>` | Alphabet (`ACGT`, `ACGU`, or `PROTEIN`) | `ACGT` |
| `--rc` | Output the reverse complement of the motif (DNA/RNA only) | *(off)* |
| `--trim-edges <float>` | Trim edges with Information Content below threshold | `0.0` |
| `--min-ic <float>` | Filter out motifs with total Information Content below threshold | `0.0` |
| `-h` | Show help | |

### homer2meme (all implementations)

| Flag | Description | Default |
|---|---|---|
| `-i <file>` | Input HOMER motif file (`-` for stdin, `.gz` supported, `.json` supported) | *(required)* |
| `-e <string>` | Extract only the specified motif by ID or description | *(all motifs)* |
| `-a <float>` | Pseudocount for log-odds to probability conversion | `0.01` |
| `-f, --format <fmt>` | Input format: `homer` or `json` | `homer` |
| `--input-format <fmt>` | Matrix type: `auto`, `logodds`, or `probability` | `auto` |
| `--rc` | Output the reverse complement of the motif (DNA/RNA only) | *(off)* |
| `--trim-edges <float>` | Trim edges with Information Content below threshold | `0.0` |
| `--min-ic <float>` | Filter out motifs with total Information Content below threshold | `0.0` |
| `-h` | Show help | |

Note: `homer2meme` auto-detects log-odds vs probability rows by checking whether row sum is near 1.0 (`[0.98, 1.02]`). This is a practical heuristic and may be ambiguous for edge-case inputs whose log-odds rows also sum near 1.0.

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

## Motif Operations

### Reverse Complement (`--rc`)

Reverses the matrix order and swaps A↔T, C↔G columns. Automatically appends `_RC` to the motif ID.

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
from motif_bridge.core import Motif

# Create a motif from a probability matrix
motif = Motif(
    id="MA0021.1",
    description="CTCF/JASPAR2026",
    matrix=[[0.35, 0.15, 0.16, 0.34], [0.10, 0.40, 0.40, 0.10]],
    alphabet="ACGT"
)

# Calculate Information Content
print(motif.total_ic())

# Reverse complement
motif.reverse_complement()

# Trim low-IC edges
motif.trim_edges(threshold=0.5)

# Serialize to dict
data = motif.to_dict()
```

### Rust Crate

Add to your `Cargo.toml` (using git dependency since `v0.2.0` is not yet published):

```toml
[dependencies]
motif-bridge = { git = "https://github.com/zengyuanzhao/motif-bridge" }
```

```rust
use motif_bridge::Motif;

let mut motif = Motif::new(
    "MA0021.1".to_string(),
    "CTCF/JASPAR2026".to_string(),
    vec![
        vec![0.35, 0.15, 0.16, 0.34],
        vec![0.10, 0.40, 0.40, 0.10],
    ],
    "ACGT".to_string(),
);

let ic = motif.total_ic();
motif.reverse_complement().unwrap();
motif.trim_edges(0.5);
motif.print_homer(0.25, 4.0);
```

---

## Metadata and Round-trip Limits

The current converters focus on standard DNA motif matrices. During conversion, some metadata fields are written as defaults:

| Direction | Metadata behavior |
|---|---|
| MEME to HOMER | HOMER `log-p`, `pseudo`, and `sites` are written as `0`; threshold is recalculated from the matrix |
| HOMER to MEME | MEME `nsites` is written as `20`; `E` is written as `0` |
| Both directions | Matrix values are printed to six decimal places |

For this reason, round-trip conversion should be interpreted as matrix-level consistency rather than byte-for-byte recovery of the original file.

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
| 2. Single motif extraction | Test `-e` flag for extracting a named motif |
| 3. stdin pipeline | Test `cat file | tool -i -` |
| 4. gzip input | Test automatic decompression of `.gz` inputs |
| 5. Round-trip | Validate matrix-level consistency for `meme→homer→meme` and `homer→meme→homer` |
| 6. Log-odds conversion | Verify log-odds to probability conversion |
| 7. Format compliance | Validate MEME headers, HOMER column counts, probability row sums |
| 8. JSON I/O | Validate JSON output structure and round-trip |
| 9. Explicit input format | Test `--input-format` flag (auto/logodds/probability) |
| 10. Alphabet support | Test `--alphabet` flag for RNA/Protein motifs |
| 11. Motif operations | Test `--rc`, `--trim-edges`, and `--min-ic` flags |

### Latest test result

| Metric | Value |
|---|---|
| Total checks | 52 |
| Passed | **52** ✅ |
| Failed | 0 |
| Skipped | 0 |

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
| Quickest install (conda environment) | `conda install -c bioconda motif-bridge` |
| Server without compiler, Perl available | `perl_scripts/` |
| Conda / Python environment | `python_scripts/` |
| Large-scale batch processing | `rust_scripts/` |
| `.gz` compressed input | All three implementations |
| Fastest runtime on tested server data | `rust_scripts/` |

---
