# motif-bridge

A cross-language toolkit for bridging the **MEME** and **HOMER** motif analysis ecosystems.

MEME and HOMER are the two dominant platforms for ChIP-seq motif analysis, but their motif file formats are fundamentally incompatible. This project provides bidirectional, lossless conversion tools implemented in **Perl**, **Python**, and **Rust**, allowing seamless interoperability between the two ecosystems.

---

## Why This Exists

| | MEME Suite | HOMER |
|---|---|---|
| **Motif format** | `.meme` — probability matrix with header block | `.motif` — log-odds or probability, tab-separated |
| **Threshold** | E-value / p-value based | Per-motif log-odds score |
| **Typical use** | *de novo* motif discovery, FIMO scanning | ChIP-seq peak annotation, known motif enrichment |
| **Downstream tools** | FIMO, MAST, TOMTOM | `findMotifsGenome.pl`, `annotatePeaks.pl` |

Converting between the two lets you:
- Use MEME-discovered motifs directly in HOMER annotation pipelines
- Scan HOMER motif databases with FIMO/MAST
- Build unified motif databases (e.g., JASPAR → both platforms)

---

## Repository Layout

```
motif-bridge/
├── README.md
├── LICENSE
├── perl_scripts/          # Perl 5 — runs on any server, no compilation
│   ├── meme2homer.pl
│   ├── homer2meme.pl
│   └── README.md
├── python_scripts/        # Python 3.6+ — zero external dependencies
│   ├── meme2homer.py
│   ├── homer2meme.py
│   └── README.md
└── rust_scripts/          # Rust — highest throughput, Cargo project
    ├── Cargo.toml
    ├── Cargo.lock
    ├── src/bin/
    │   ├── meme2homer.rs
    │   └── homer2meme.rs
    └── README.md
```

All three language implementations share **identical CLI flags** and produce matching output for valid MEME/HOMER inputs — choose the one that fits your environment.

---

## Quick Start

### MEME → HOMER

```bash
# Perl
perl perl_scripts/meme2homer.pl -i motifs.meme -j JASPAR2026 > motifs.homer

# Python
python3 python_scripts/meme2homer.py -i motifs.meme -j JASPAR2026 > motifs.homer

# Rust (build first)
cd rust_scripts && cargo build --release && cd ..
./rust_scripts/target/release/meme2homer -i motifs.meme -j JASPAR2026 > motifs.homer
```

### HOMER → MEME

```bash
# Perl
perl perl_scripts/homer2meme.pl -i motifs.homer > motifs.meme

# Python
python3 python_scripts/homer2meme.py -i motifs.homer > motifs.meme

# Rust (build first)
cd rust_scripts && cargo build --release && cd ..
./rust_scripts/target/release/homer2meme -i motifs.homer > motifs.meme
```

### Building the Rust binaries

```bash
cd rust_scripts
cargo build --release
# Produces: rust_scripts/target/release/meme2homer
#           rust_scripts/target/release/homer2meme

# Or install to PATH:
cargo install --path .
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
| `-h` | Show help | — |

### homer2meme (all implementations)

| Flag | Description | Default |
|---|---|---|
| `-i <file>` | Input HOMER motif file (`-` for stdin, `.gz` supported) | *(required)* |
| `-e <string>` | Extract only the specified motif by ID or description | *(all motifs)* |
| `-a <float>` | Pseudocount for log-odds → probability conversion | `0.01` |
| `-h` | Show help | — |

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
threshold = sum_over_positions( log2(max_prob / background) ) − t_offset
threshold = max(threshold, 0)
```

---

## Validation

The project was validated on a server using the following real motif datasets:

| File | Description |
|---|---|
| `homer.known.motifs` | HOMER built-in known motif library (436 motifs) |
| `JASPAR2024_small.meme` | First 200 lines of `JASPAR2024_vertebrates.meme` (12 motifs), used for rapid iteration |
| `JASPAR2024_vertebrates.meme` | Full JASPAR 2024 vertebrate motif database (879 motifs, 469 KB) |

### Test coverage

| Test stage | Description |
|---|---|
| 0. Environment check | Verify data files, Perl, Python, and Rust availability |
| 1. meme2homer | MEME → HOMER conversion on small file (12 motifs) |
| 2. homer2meme | HOMER → MEME conversion on known motif library (436 motifs) |
| 3. Output consistency | Cross-language diff between Perl, Python, and Rust outputs |
| 4. Round-trip | Validate `meme→homer→meme` and `homer→meme→homer` losslessness |
| 5. Single motif extraction | Test `-e` flag for extracting a named motif |
| 6. stdin pipeline | Test `cat file \| tool -i -` |
| 7. gzip input | Test automatic decompression of `.gz` inputs |
| 8. Large-file performance | Benchmark all three implementations on 879-motif dataset |
| 9. Format compliance | Validate MEME headers, HOMER column counts, probability row sums |

### Latest server-side result (2026-04-03)

| Metric | Value |
|---|---|
| Total checks | 44 |
| Passed | **44** ✅ |
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
| Server without compiler, Perl available | `perl_scripts/` |
| Conda / Python environment | `python_scripts/` |
| Large-scale batch processing | `rust_scripts/` |
| `.gz` compressed input | All three implementations |
| Fastest runtime on tested server data | `rust_scripts/` |

---

## License

MIT
