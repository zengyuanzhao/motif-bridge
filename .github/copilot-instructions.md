# Copilot instructions for motif-bridge

## Build, test, and lint commands

### Prerequisites
- Python 3.8+, Perl 5 (with `IO::Uncompress::Gunzip`), Rust 1.70+.

### Build
- Rust release binaries:
  - `cd rust_scripts && cargo build --release`
  - Outputs: `rust_scripts/target/release/meme2homer` and `rust_scripts/target/release/homer2meme`

### Test
- Full regression suite (cross-language consistency + fixtures + round-trip + format checks):
  - `bash test_motif_bridge.sh`

- Single regression check (example: one meme2homer expected-output test):
  - `python3 python_scripts/meme2homer.py -i fixtures/test.meme -j JASPAR2026 > /tmp/m2h.out && diff -u /tmp/m2h.out fixtures/expected_meme2homer.homer`

### Lint / format checks
- Rust:
  - `cd rust_scripts && cargo fmt --check`
  - `cd rust_scripts && cargo clippy -- -D warnings`
- Python:
  - `ruff check python_scripts/`
  - `ruff format --check python_scripts/`
- Perl:
  - No automated linter is configured; follow existing script style.

## High-level architecture

`motif-bridge` implements the same two CLI tools in three languages:
- `meme2homer` (MEME -> HOMER)
- `homer2meme` (HOMER -> MEME)

Implementation layout:
- Python: `python_scripts/meme2homer.py`, `python_scripts/homer2meme.py`
- Perl: `perl_scripts/meme2homer.pl`, `perl_scripts/homer2meme.pl`
- Rust: `rust_scripts/src/bin/meme2homer.rs`, `rust_scripts/src/bin/homer2meme.rs`

Core design across languages:
- The CLIs are kept behaviorally aligned (same core flags and conversion behavior).
- Both directions support plain files, `.gz`, and stdin (`-i -`).
- `meme2homer` parses MEME motifs and emits HOMER text (or JSON with `-f json`).
- `homer2meme` parses HOMER (or JSON with `-f json`) and emits MEME text.
- `test_motif_bridge.sh` is the parity gate: it compares Python/Perl/Rust outputs, validates fixture expectations, checks round-trip matrix consistency, and verifies format compliance.

## Key repository conventions

- Treat Python implementation as the reference when changing conversion logic; then port the same behavior to Perl and Rust.
- Keep output parity across languages for standard DNA motif inputs; changes to conversion behavior should be reflected in fixture expectations and `test_motif_bridge.sh`.
- Preserve output formatting conventions used by all implementations:
  - Matrix values printed to 6 decimal places.
  - HOMER headers use 6 tab-separated fields (`>id<TAB>description<TAB>threshold<TAB>0<TAB>0<TAB>0`).
  - MEME output uses fixed header scaffolding, with generated defaults (for example `nsites=20`, `E=0` in homer2meme output).
- `-e` extraction is expected to match by motif ID and by motif name/description.
- In `homer2meme` auto mode, log-odds vs probability detection follows the row-sum heuristic (`sum in [0.98, 1.02]` treated as probability).
