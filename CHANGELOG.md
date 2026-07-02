# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Zenodo archive metadata in `.zenodo.json`.

### Changed
- Expanded `CITATION.cff` with repository URL, abstract, and citation keywords.

## [0.3.0] - 2026-07-03

### Added
- `--strict` mode for Python, Perl, and Rust CLIs to fail on malformed matrix rows, invalid probability values, probability row-sum drift, MEME `alength=`/alphabet conflicts, and ambiguous HOMER auto-detection rows.
- `CITATION.cff` software citation metadata.
- Python package URLs and `dev` optional dependencies in `pyproject.toml`.
- `--version` flag for Python, Perl, and Rust CLIs.
- Regression coverage for CLI version reporting, Perl `VERSION_FROM` parsing, and MEME `alength=` conflicts.
- Python and Rust unit tests for core IC, log-odds, reverse-complement, trimming, and format-writer behavior.
- Explicit pytest configuration for the `tests/` suite.
- README caveats for downstream interpretation of recalculated HOMER thresholds, `--input-format auto`, fixed MEME metadata, background handling, and formatted row-sum drift.
- Warnings when recalculated HOMER thresholds clip to `0`, and when nonnegative rows fall near the `--input-format auto` row-sum boundary.
- Default-off safety flags: `--nsites`, `--evalue`, `--renormalize`, `--keep-threshold`, and comma-separated background vectors for `-b`.
- Regression coverage for zero-threshold warnings, auto-detection gray-zone warnings, MEME metadata overrides, row renormalization, background vectors, invalid background vectors, tied maximum scoring with asymmetric backgrounds, and JSON metadata preservation.

### Fixed
- `homer2meme` now writes valid comma-separated background vectors from `-b` into the MEME `Background letter frequencies` header instead of always writing a uniform background.
- JSON output now uses top-level `"format": "motif-bridge-json"` instead of hard-coding `"source": "meme"`.
- Perl PROTEIN information-content calculations now use the 20-letter maximum entropy.
- MEME `ALPHABET=` is auto-detected across Python, Perl, and Rust readers, with `--alphabet` retained as an override.
- Rust MEME parsing now warns on malformed matrix rows, matching Python and Perl behavior.
- HOMER score calculation is centralized on the motif model instead of duplicated in Python I/O code.
- Rust JSON string output now uses `serde_json::to_string` escaping.
- PROTEIN MEME output no longer emits DNA/RNA `strands: + -` metadata.
- Reverse-complement requests for unsupported alphabets now skip motifs with warnings instead of failing hard.
- Perl package version discovery now reads `our $VERSION` from `perl_scripts/meme2homer.pl`.
- Perl stderr is consistently emitted as UTF-8 across both converters.
- MEME readers now warn on `alength=`/`ALPHABET=` width conflicts and keep the alphabet-derived width instead of dropping valid rows.
- Python CLI version reporting now prefers the source-tree `pyproject.toml` version when running from a checkout, avoiding stale editable-install metadata.
- Shell regression checks now run inline Python validators inside `if` conditions so failures are reported by the suite under `set -euo pipefail`.
- Version regression tests compare CLI output to package metadata instead of hard-coding one release number.
- Background vectors now reject length mismatches, values outside `(0,1]`, and sums farther than `1e-3` from `1.0` instead of silently affecting thresholds or log-odds restoration.
- JSON output/input now carries optional `threshold`, `nsites`, and `evalue` motif metadata when present.
- Rust HOMER threshold scoring now matches Python and Perl by choosing the first maximum-probability column when a row has tied maxima.

## [0.2.0] - 2026-05-21

### Added
- CLI regression test suite (`test_motif_bridge.sh`) with 84 checks across Python, Perl, and Rust implementations
- Fixture files for testing (`fixtures/`)
- GitHub Actions CI workflow (`.github/workflows/ci.yml`)
- `.gitignore` for Rust, Python, and OS artifacts
- JSON output format (`-f json`) for meme2homer across all three languages
- JSON input format (`-f json`) for homer2meme (Python/Rust/Perl)
- `--input-format` flag to explicitly specify matrix type (auto/logodds/probability)
- `--alphabet` flag for RNA (ACGU) and PROTEIN motif support
- `--alphabet` and `--background` flags for homer2meme to align conversion behavior with meme2homer
- Test stages 10 (--input-format), 11 (--alphabet), 12 (motif operations), 13 (MOTIF word boundary), and 14 (negative matrix warnings) added to the test suite

### Changed
- Refactored Rust error handling to use `Result` types and `?` operator for better testability
- Applied `cargo fmt` and `cargo clippy` standards to Rust code
- Added `serde_json` dependency to Rust crate for JSON parsing
- Replaced fragile regex-based JSON parsing in Perl with `JSON::PP` core module
- Unified Python version requirement to 3.8+ across docstrings and pyproject.toml

## [0.1.0] - 2026-05-16

### Added
- Initial release of motif-bridge
- Bidirectional MEME ↔ HOMER motif converters in Perl, Python, and Rust
- Support for gzip compressed input and stdin piping
- Single motif extraction via `-e` flag
- Published to Bioconda and Crates.io
