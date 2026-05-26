# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `--version` flag for Python, Perl, and Rust CLIs.
- Regression coverage for CLI version reporting, Perl `VERSION_FROM` parsing, and MEME `alength=` conflicts.
- Python and Rust unit tests for core IC, log-odds, reverse-complement, trimming, and format-writer behavior.

### Fixed
- Perl package version discovery now reads `our $VERSION` from `perl_scripts/meme2homer.pl`.
- Perl stderr is consistently emitted as UTF-8 across both converters.
- MEME readers now warn on `alength=`/`ALPHABET=` width conflicts and keep the alphabet-derived width instead of dropping valid rows.
- Shell regression checks now run inline Python validators inside `if` conditions so failures are reported by the suite under `set -euo pipefail`.

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
