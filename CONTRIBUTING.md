# Contributing to motif-bridge

Thank you for your interest in contributing! This project provides identical functionality in three languages (Perl, Python, Rust) to serve different user environments.

## Development Setup

### Prerequisites
- **Python 3.8+** (for Python scripts)
- **Perl 5** with `IO::Uncompress::Gunzip` (for Perl scripts)
- **Rust 1.70+** (for Rust binaries)

### Running Tests
```bash
bash test_motif_bridge.sh
```
This runs 84 regression checks covering cross-language consistency, stdin/gzip support, round-trip conversion, JSON I/O, and format compliance.

To run specific test stages:
```bash
TEST_STAGE=1 bash test_motif_bridge.sh          # Stage 1 only
TEST_STAGE=1,3,5 bash test_motif_bridge.sh      # Stages 1, 3, 5
```

### Local Development with micromamba
```bash
micromamba create -n motif-bridge python=3.10 ruff=0.9.6 perl -c conda-forge -y
micromamba run -n motif-bridge ruff check python_scripts/
micromamba run -n motif-bridge ruff format --check python_scripts/
micromamba run -n motif-bridge bash test_motif_bridge.sh
```

### Code Quality
- **Rust**: `cd rust_scripts && cargo fmt && cargo clippy -- -D warnings`
- **Python**: `ruff check python_scripts/ && ruff format --check python_scripts/`
- **Perl**: No automated linter configured; follow existing style.

## Pull Request Guidelines
1. Ensure `test_motif_bridge.sh` passes locally
2. Update all three language implementations if changing conversion logic
3. Add fixture files and test cases for new features
4. Update `CHANGELOG.md` under `[Unreleased]`

## Adding a New Feature
If you add a feature (e.g., new CLI flag, format support):
1. Implement in **Python** first (reference implementation)
2. Port to **Perl** and **Rust**
3. Add cross-language diff tests in `test_motif_bridge.sh`
4. Update README.md documentation

## Reporting Issues
- Include the command run, input file (or minimal reproducible example), and expected vs actual output
- Specify which language implementation you used (Perl/Python/Rust)
