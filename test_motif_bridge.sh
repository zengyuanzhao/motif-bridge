#!/usr/bin/env bash
# test_motif_bridge.sh - CLI regression tests for motif-bridge
#
# Tests meme2homer and homer2meme across Python, Perl, and Rust implementations.
# Requires: python3 (or python), perl, cargo (for Rust build), diff, gzip
#
# Usage:
#   bash test_motif_bridge.sh              # Run all tests
#   TEST_STAGE=1 bash test_motif_bridge.sh # Run only stage 1
#   TEST_STAGE=1,3,5 bash test_motif_bridge.sh  # Run stages 1, 3, 5
#   TEST_STAGES=python bash test_motif_bridge.sh # Run only Python-related tests
#
# Available stages:
#   1  - meme2homer cross-language consistency
#   2  - homer2meme cross-language consistency
#   3  - Single motif extraction (-e flag)
#   4  - stdin support
#   5  - gzip input support
#   6  - Round-trip consistency
#   7  - Log-odds to probability conversion
#   8  - Format compliance
#   9  - JSON output format
#   10 - --input-format explicit specification
#   11 - --alphabet support
#   12 - Motif Operations (--rc, --trim-edges, --min-ic)
#   13 - MEME MOTIF word-boundary parsing
#   14 - Negative matrix value warnings
#   15 - Version and parser metadata regressions
#
# Exit codes:
#   0 - all tests passed
#   1 - one or more tests failed

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIXTURES="$SCRIPT_DIR/fixtures"
WORK_DIR=$(mktemp -d)
trap 'rm -rf "$WORK_DIR"' EXIT
export PYTHONPATH="$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}"

PASS=0
FAIL=0
TOTAL=0

# ---------------------------------------------------------------------------
# Stage filtering
# ---------------------------------------------------------------------------

STAGE_FILTER="${TEST_STAGE:-${TEST_STAGES:-all}}"

should_run_stage() {
    local stage="$1"
    if [[ "$STAGE_FILTER" == "all" ]]; then
        return 0
    fi
    # Check if stage number is in the comma-separated list
    IFS=',' read -ra stages <<< "$STAGE_FILTER"
    for s in "${stages[@]}"; do
        if [[ "$s" == "$stage" ]]; then
            return 0
        fi
    done
    return 1
}

run_stage() {
    local stage="$1"
    local label="$2"
    if should_run_stage "$stage"; then
        echo "=== Test $stage: $label ==="
        return 0
    else
        echo "=== Test $stage: $label === (skipped)"
        return 1
    fi
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

pass() {
    PASS=$((PASS + 1))
    TOTAL=$((TOTAL + 1))
    echo "  PASS: $1"
}

fail() {
    FAIL=$((FAIL + 1))
    TOTAL=$((TOTAL + 1))
    echo "  FAIL: $1"
    if [ -n "${2:-}" ]; then
        echo "        $2"
    fi
}

check_diff() {
    local file1="$1" file2="$2" label="$3"
    if diff -q "$file1" "$file2" > /dev/null 2>&1; then
        pass "$label"
    else
        fail "$label" "$(diff "$file1" "$file2" | head -5)"
    fi
}

check_exit() {
    local exit_code="$1" label="$2" expected="${3:-0}"
    if [ "$exit_code" -eq "$expected" ]; then
        pass "$label"
    else
        fail "$label" "expected exit $expected, got $exit_code"
    fi
}

ensure_m2h_baseline() {
    if [ ! -f "$WORK_DIR/py_m2h.homer" ]; then
        python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test.meme" -j JASPAR2026 > "$WORK_DIR/py_m2h.homer" 2>&1
        perl "$PERL/meme2homer.pl" -i "$FIXTURES/test.meme" -j JASPAR2026 > "$WORK_DIR/pl_m2h.homer" 2>&1
        if [ -n "$RUST_BIN" ]; then
            "$RUST_BIN/meme2homer" -i "$FIXTURES/test.meme" -j JASPAR2026 > "$WORK_DIR/rs_m2h.homer" 2>&1
        fi
    fi
}

ensure_h2m_baseline() {
    if [ ! -f "$WORK_DIR/py_h2m.meme" ]; then
        python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test.homer" > "$WORK_DIR/py_h2m.meme" 2>&1
        perl "$PERL/homer2meme.pl" -i "$FIXTURES/test.homer" > "$WORK_DIR/pl_h2m.meme" 2>&1
        if [ -n "$RUST_BIN" ]; then
            "$RUST_BIN/homer2meme" -i "$FIXTURES/test.homer" > "$WORK_DIR/rs_h2m.meme" 2>&1
        fi
    fi
}

# ---------------------------------------------------------------------------
# Prerequisites
# ---------------------------------------------------------------------------

echo "=== Prerequisites ==="

if ! command -v python3 > /dev/null 2>&1; then
    if command -v python > /dev/null 2>&1; then
        python3() { python "$@"; }
    else
        echo "SKIP: python3 not found"
        exit 1
    fi
fi
command -v perl > /dev/null 2>&1 || { echo "SKIP: perl not found"; exit 1; }
command -v diff > /dev/null 2>&1 || { echo "SKIP: diff not found"; exit 1; }

PYTHON="$SCRIPT_DIR/python_scripts"
PERL="$SCRIPT_DIR/perl_scripts"
RUST_BIN=""

if command -v cargo > /dev/null 2>&1; then
    echo "  Building Rust binaries..."
    (cd "$SCRIPT_DIR/rust_scripts" && cargo build --release > /dev/null 2>&1)
    RUST_BIN="$SCRIPT_DIR/rust_scripts/target/release"
fi

echo "  python3: $(python3 --version 2>&1)"
echo "  perl:    $(perl --version 2>&1 | head -1)"
if [ -n "$RUST_BIN" ]; then
    echo "  rust:    $(rustc --version 2>&1)"
else
    echo "  rust:    not available (cargo not found)"
fi
echo ""

# ---------------------------------------------------------------------------
# Test 1: meme2homer cross-language consistency
# ---------------------------------------------------------------------------

if run_stage 1 "meme2homer cross-language consistency"; then

python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test.meme" -j JASPAR2026 > "$WORK_DIR/py_m2h.homer" 2>&1
perl "$PERL/meme2homer.pl" -i "$FIXTURES/test.meme" -j JASPAR2026 > "$WORK_DIR/pl_m2h.homer" 2>&1
check_diff "$WORK_DIR/py_m2h.homer" "$WORK_DIR/pl_m2h.homer" "Python vs Perl meme2homer"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test.meme" -j JASPAR2026 > "$WORK_DIR/rs_m2h.homer" 2>&1
    check_diff "$WORK_DIR/py_m2h.homer" "$WORK_DIR/rs_m2h.homer" "Python vs Rust meme2homer"
fi

# Compare with expected output
check_diff "$WORK_DIR/py_m2h.homer" "$FIXTURES/expected_meme2homer.homer" "meme2homer vs expected"

echo ""
fi

# ---------------------------------------------------------------------------
# Test 2: homer2meme cross-language consistency
# ---------------------------------------------------------------------------

if run_stage 2 "homer2meme cross-language consistency"; then

python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test.homer" > "$WORK_DIR/py_h2m.meme" 2>&1
perl "$PERL/homer2meme.pl" -i "$FIXTURES/test.homer" > "$WORK_DIR/pl_h2m.meme" 2>&1
check_diff "$WORK_DIR/py_h2m.meme" "$WORK_DIR/pl_h2m.meme" "Python vs Perl homer2meme"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/homer2meme" -i "$FIXTURES/test.homer" > "$WORK_DIR/rs_h2m.meme" 2>&1
    check_diff "$WORK_DIR/py_h2m.meme" "$WORK_DIR/rs_h2m.meme" "Python vs Rust homer2meme"
fi

check_diff "$WORK_DIR/py_h2m.meme" "$FIXTURES/expected_homer2meme.meme" "homer2meme vs expected"

echo ""
fi

# ---------------------------------------------------------------------------
# Test 3: Single motif extraction (-e flag)
# ---------------------------------------------------------------------------

if run_stage 3 "Single motif extraction (-e flag)"; then

python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test.meme" -j JASPAR2026 -e MA0021.1 > "$WORK_DIR/py_extract.homer" 2>&1
perl "$PERL/meme2homer.pl" -i "$FIXTURES/test.meme" -j JASPAR2026 -e MA0021.1 > "$WORK_DIR/pl_extract.homer" 2>&1
check_diff "$WORK_DIR/py_extract.homer" "$WORK_DIR/pl_extract.homer" "Python vs Perl extract"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test.meme" -j JASPAR2026 -e MA0021.1 > "$WORK_DIR/rs_extract.homer" 2>&1
    check_diff "$WORK_DIR/py_extract.homer" "$WORK_DIR/rs_extract.homer" "Python vs Rust extract"
fi

check_diff "$WORK_DIR/py_extract.homer" "$FIXTURES/expected_extract.homer" "extract vs expected"

# homer2meme extract
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test.homer" -e MA0021.1 > "$WORK_DIR/py_h2m_extract.meme" 2>&1
perl "$PERL/homer2meme.pl" -i "$FIXTURES/test.homer" -e MA0021.1 > "$WORK_DIR/pl_h2m_extract.meme" 2>&1
check_diff "$WORK_DIR/py_h2m_extract.meme" "$WORK_DIR/pl_h2m_extract.meme" "Python vs Perl homer2meme extract"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/homer2meme" -i "$FIXTURES/test.homer" -e MA0021.1 > "$WORK_DIR/rs_h2m_extract.meme" 2>&1
    check_diff "$WORK_DIR/py_h2m_extract.meme" "$WORK_DIR/rs_h2m_extract.meme" "Python vs Rust homer2meme extract"
fi

echo ""
fi

# ---------------------------------------------------------------------------
# Test 4: stdin support
# ---------------------------------------------------------------------------

if run_stage 4 "stdin support"; then

ensure_m2h_baseline
ensure_h2m_baseline

cat "$FIXTURES/test.meme" | python3 "$PYTHON/meme2homer.py" -i - -j JASPAR2026 > "$WORK_DIR/stdin_py.homer" 2>&1
check_diff "$WORK_DIR/stdin_py.homer" "$WORK_DIR/py_m2h.homer" "Python meme2homer stdin"

cat "$FIXTURES/test.meme" | perl "$PERL/meme2homer.pl" -i - -j JASPAR2026 > "$WORK_DIR/stdin_pl.homer" 2>&1
check_diff "$WORK_DIR/stdin_pl.homer" "$WORK_DIR/pl_m2h.homer" "Perl meme2homer stdin"

if [ -n "$RUST_BIN" ]; then
    cat "$FIXTURES/test.meme" | "$RUST_BIN/meme2homer" -i - -j JASPAR2026 > "$WORK_DIR/stdin_rs.homer" 2>&1
    check_diff "$WORK_DIR/stdin_rs.homer" "$WORK_DIR/rs_m2h.homer" "Rust meme2homer stdin"
fi

# homer2meme stdin
cat "$FIXTURES/test.homer" | python3 "$PYTHON/homer2meme.py" -i - > "$WORK_DIR/stdin_py_h2m.meme" 2>&1
check_diff "$WORK_DIR/stdin_py_h2m.meme" "$WORK_DIR/py_h2m.meme" "Python homer2meme stdin"

cat "$FIXTURES/test.homer" | perl "$PERL/homer2meme.pl" -i - > "$WORK_DIR/stdin_pl_h2m.meme" 2>&1
check_diff "$WORK_DIR/stdin_pl_h2m.meme" "$WORK_DIR/pl_h2m.meme" "Perl homer2meme stdin"

echo ""
fi

# ---------------------------------------------------------------------------
# Test 5: gzip input support
# ---------------------------------------------------------------------------

if run_stage 5 "gzip input support"; then

ensure_m2h_baseline
ensure_h2m_baseline

python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test.meme.gz" -j JASPAR2026 > "$WORK_DIR/gz_py.homer" 2>&1
check_diff "$WORK_DIR/gz_py.homer" "$WORK_DIR/py_m2h.homer" "Python meme2homer gzip"

perl "$PERL/meme2homer.pl" -i "$FIXTURES/test.meme.gz" -j JASPAR2026 > "$WORK_DIR/gz_pl.homer" 2>&1
check_diff "$WORK_DIR/gz_pl.homer" "$WORK_DIR/pl_m2h.homer" "Perl meme2homer gzip"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test.meme.gz" -j JASPAR2026 > "$WORK_DIR/gz_rs.homer" 2>&1
    check_diff "$WORK_DIR/gz_rs.homer" "$WORK_DIR/rs_m2h.homer" "Rust meme2homer gzip"
fi

# homer2meme gzip
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test.homer.gz" > "$WORK_DIR/gz_py_h2m.meme" 2>&1
check_diff "$WORK_DIR/gz_py_h2m.meme" "$WORK_DIR/py_h2m.meme" "Python homer2meme gzip"

perl "$PERL/homer2meme.pl" -i "$FIXTURES/test.homer.gz" > "$WORK_DIR/gz_pl_h2m.meme" 2>&1
check_diff "$WORK_DIR/gz_pl_h2m.meme" "$WORK_DIR/pl_h2m.meme" "Perl homer2meme gzip"

echo ""
fi

# ---------------------------------------------------------------------------
# Test 6: Round-trip consistency (meme->homer->meme)
# ---------------------------------------------------------------------------

if run_stage 6 "Round-trip consistency"; then

# meme -> homer -> meme
python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test.meme" -j JASPAR2026 > "$WORK_DIR/rt_m2h.homer" 2>&1
python3 "$PYTHON/homer2meme.py" -i "$WORK_DIR/rt_m2h.homer" > "$WORK_DIR/rt_m2h2m.meme" 2>&1

# Compare matrix values (ignore header metadata differences)
if python3 -c "
import sys
def parse_meme(path):
    motifs = []
    with open(path) as f:
        in_matrix = False
        current = None
        for line in f:
            s = line.strip()
            if s.startswith('MOTIF'):
                current = {'id': s, 'matrix': []}
                motifs.append(current)
                in_matrix = False
            elif s.startswith('letter-probability'):
                in_matrix = True
            elif s.startswith('//') or s.startswith('MOTIF') and current:
                in_matrix = False
            elif in_matrix and s and (s[0].isdigit() or s.startswith('.')):
                current['matrix'].append([float(x) for x in s.split()])
    return motifs

orig = parse_meme('$FIXTURES/test.meme')
rt = parse_meme('$WORK_DIR/rt_m2h2m.meme')
ok = True
for o, r in zip(orig, rt):
    if len(o['matrix']) != len(r['matrix']):
        print(f'Matrix length mismatch for {o[\"id\"]}: {len(o[\"matrix\"])} vs {len(r[\"matrix\"])}')
        ok = False
        continue
    for i, (orow, rrow) in enumerate(zip(o['matrix'], r['matrix'])):
        for j, (ov, rv) in enumerate(zip(orow, rrow)):
            if abs(ov - rv) > 1e-5:
                print(f'Value mismatch at {o[\"id\"]} row {i} col {j}: {ov} vs {rv}')
                ok = False
if ok:
    print('OK')
    sys.exit(0)
else:
    sys.exit(1)
" > "$WORK_DIR/rt_result.txt" 2>&1; then
    pass "Round-trip meme->homer->meme matrix consistency"
else
    fail "Round-trip meme->homer->meme matrix consistency" "$(cat "$WORK_DIR/rt_result.txt")"
fi

# homer -> meme -> homer
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test.homer" > "$WORK_DIR/rt_h2m.meme" 2>&1
python3 "$PYTHON/meme2homer.py" -i "$WORK_DIR/rt_h2m.meme" -j JASPAR2026 > "$WORK_DIR/rt_h2m2h.homer" 2>&1

# Compare matrix values
if python3 -c "
import sys
def parse_homer(path):
    motifs = []
    with open(path) as f:
        current = None
        for line in f:
            s = line.strip()
            if s.startswith('>'):
                current = {'header': s, 'matrix': []}
                motifs.append(current)
            elif current and s and (s[0].isdigit() or s.startswith('.')):
                current['matrix'].append([float(x) for x in s.split()])
    return motifs

orig = parse_homer('$FIXTURES/test.homer')
rt = parse_homer('$WORK_DIR/rt_h2m2h.homer')
ok = True
for o, r in zip(orig, rt):
    if len(o['matrix']) != len(r['matrix']):
        print(f'Matrix length mismatch: {len(o[\"matrix\"])} vs {len(r[\"matrix\"])}')
        ok = False
        continue
    for i, (orow, rrow) in enumerate(zip(o['matrix'], r['matrix'])):
        for j, (ov, rv) in enumerate(zip(orow, rrow)):
            if abs(ov - rv) > 1e-5:
                print(f'Value mismatch at row {i} col {j}: {ov} vs {rv}')
                ok = False
if ok:
    print('OK')
    sys.exit(0)
else:
    sys.exit(1)
" > "$WORK_DIR/rt_result2.txt" 2>&1; then
    pass "Round-trip homer->meme->homer matrix consistency"
else
    fail "Round-trip homer->meme->homer matrix consistency" "$(cat "$WORK_DIR/rt_result2.txt")"
fi

echo ""
fi

# ---------------------------------------------------------------------------
# Test 7: Log-odds to probability conversion
# ---------------------------------------------------------------------------

if run_stage 7 "Log-odds to probability conversion"; then

python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test_logodds.homer" > "$WORK_DIR/py_logodds.meme" 2>&1
perl "$PERL/homer2meme.pl" -i "$FIXTURES/test_logodds.homer" > "$WORK_DIR/pl_logodds.meme" 2>&1
check_diff "$WORK_DIR/py_logodds.meme" "$WORK_DIR/pl_logodds.meme" "Python vs Perl log-odds conversion"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/homer2meme" -i "$FIXTURES/test_logodds.homer" > "$WORK_DIR/rs_logodds.meme" 2>&1
    check_diff "$WORK_DIR/py_logodds.meme" "$WORK_DIR/rs_logodds.meme" "Python vs Rust log-odds conversion"
fi

# Verify output rows sum to ~1.0 (probability)
if python3 -c "
import sys
ok = True
with open('$WORK_DIR/py_logodds.meme') as f:
    in_matrix = False
    for line in f:
        s = line.strip()
        if s.startswith('letter-probability'):
            in_matrix = True
            continue
        if s.startswith('MOTIF') or s.startswith('//') or s.startswith('MEME') or s.startswith('ALPHABET') or s.startswith('strands') or s.startswith('Background') or s.startswith('A 0.25') or s == '':
            in_matrix = False
            continue
        if in_matrix and s:
            vals = [float(x) for x in s.split()]
            if abs(sum(vals) - 1.0) > 0.02:
                print(f'Row sum not near 1.0: {sum(vals):.6f} for row: {s}')
                ok = False
if ok:
    print('OK')
    sys.exit(0)
else:
    sys.exit(1)
" > "$WORK_DIR/logodds_check.txt" 2>&1; then
    pass "Log-odds output rows sum to ~1.0"
else
    fail "Log-odds output rows sum to ~1.0" "$(cat "$WORK_DIR/logodds_check.txt")"
fi

# Verify --background affects log-odds conversion consistently
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test_logodds.homer" --input-format logodds -b 0.2 > "$WORK_DIR/py_logodds_bg.meme" 2>&1
perl "$PERL/homer2meme.pl" -i "$FIXTURES/test_logodds.homer" --input-format logodds -b 0.2 > "$WORK_DIR/pl_logodds_bg.meme" 2>&1
check_diff "$WORK_DIR/py_logodds_bg.meme" "$WORK_DIR/pl_logodds_bg.meme" "Python vs Perl log-odds background"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/homer2meme" -i "$FIXTURES/test_logodds.homer" --input-format logodds -b 0.2 > "$WORK_DIR/rs_logodds_bg.meme" 2>&1
    check_diff "$WORK_DIR/py_logodds_bg.meme" "$WORK_DIR/rs_logodds_bg.meme" "Python vs Rust log-odds background"
fi

check_diff "$WORK_DIR/py_logodds_bg.meme" "$FIXTURES/expected_logodds_bg_0p2.meme" "--background log-odds conversion"

echo ""
fi

# ---------------------------------------------------------------------------
# Test 8: Format compliance
# ---------------------------------------------------------------------------

if run_stage 8 "Format compliance"; then

# MEME output: check header presence
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test.homer" > "$WORK_DIR/format_meme.meme" 2>&1
if python3 -c "
import sys
with open('$WORK_DIR/format_meme.meme') as f:
    content = f.read()
checks = [
    ('MEME version 4' in content, 'MEME version header'),
    ('ALPHABET= ACGT' in content, 'ALPHABET line'),
    ('strands: + -' in content, 'strands line'),
    ('Background letter frequencies' in content, 'Background frequencies header'),
    ('A 0.25 C 0.25 G 0.25 T 0.25' in content, 'Background frequencies values'),
    ('MOTIF' in content, 'MOTIF keyword'),
    ('letter-probability matrix:' in content, 'letter-probability matrix header'),
]
ok = True
for result, label in checks:
    if not result:
        print(f'Missing: {label}')
        ok = False
if ok:
    print('OK')
    sys.exit(0)
else:
    sys.exit(1)
" > "$WORK_DIR/meme_format_check.txt" 2>&1; then
    pass "MEME format header compliance"
else
    fail "MEME format header compliance" "$(cat "$WORK_DIR/meme_format_check.txt")"
fi

# HOMER output: check 6 tab-separated header fields and 4-column matrix rows
python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test.meme" -j JASPAR2026 > "$WORK_DIR/format_homer.homer" 2>&1
if python3 -c "
import sys
ok = True
with open('$WORK_DIR/format_homer.homer') as f:
    for line in f:
        s = line.strip()
        if s.startswith('>'):
            fields = s[1:].split('\t')
            if len(fields) != 6:
                print(f'HOMER header has {len(fields)} fields (expected 6): {s}')
                ok = False
        elif s and (s[0].isdigit() or s.startswith('.')):
            vals = s.split('\t')
            if len(vals) != 4:
                print(f'Matrix row has {len(vals)} cols (expected 4): {s}')
                ok = False
if ok:
    print('OK')
    sys.exit(0)
else:
    sys.exit(1)
" > "$WORK_DIR/homer_format_check.txt" 2>&1; then
    pass "HOMER format compliance (6 header fields, 4-col matrix)"
else
    fail "HOMER format compliance (6 header fields, 4-col matrix)" "$(cat "$WORK_DIR/homer_format_check.txt")"
fi

echo ""
fi

# ---------------------------------------------------------------------------
# Test 9: JSON output format
# ---------------------------------------------------------------------------

if run_stage 9 "JSON output format"; then

python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test.meme" -j JASPAR2026 -f json > "$WORK_DIR/py_json.json" 2>&1
perl "$PERL/meme2homer.pl" -i "$FIXTURES/test.meme" -j JASPAR2026 -f json > "$WORK_DIR/pl_json.json" 2>&1
check_diff "$WORK_DIR/py_json.json" "$WORK_DIR/pl_json.json" "Python vs Perl JSON output"

# Unicode descriptions should be preserved consistently
python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test_unicode_desc.meme" -j JASPAR2026 -f json > "$WORK_DIR/py_unicode.json" 2>&1
perl "$PERL/meme2homer.pl" -i "$FIXTURES/test_unicode_desc.meme" -j JASPAR2026 -f json > "$WORK_DIR/pl_unicode.json" 2>&1
check_diff "$WORK_DIR/py_unicode.json" "$WORK_DIR/pl_unicode.json" "Python vs Perl JSON unicode"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test_unicode_desc.meme" -j JASPAR2026 -f json > "$WORK_DIR/rs_unicode.json" 2>&1
    check_diff "$WORK_DIR/py_unicode.json" "$WORK_DIR/rs_unicode.json" "Python vs Rust JSON unicode"
fi

if python3 -c "
import json, sys
with open('$WORK_DIR/py_unicode.json', encoding='utf-8') as f:
    data = json.load(f)
desc = data.get('motifs', [{}])[0].get('description', '')
if desc != 'alpha_α/JASPAR2026':
    print(f'Unexpected description: {desc}')
    sys.exit(1)
print('OK')
sys.exit(0)
" > "$WORK_DIR/unicode_check.txt" 2>&1; then
    pass "Unicode descriptions preserved"
else
    fail "Unicode descriptions preserved" "$(cat "$WORK_DIR/unicode_check.txt")"
fi

# Validate JSON structure
if python3 -c "
import json, sys
try:
    with open('$WORK_DIR/py_json.json') as f:
        data = json.load(f)
    if data.get('version') != '1.0':
        print('Missing version')
        sys.exit(1)
    if data.get('source') != 'meme':
        print('Missing source')
        sys.exit(1)
    motifs = data.get('motifs', [])
    if len(motifs) != 2:
        print(f'Expected 2 motifs, got {len(motifs)}')
        sys.exit(1)
    for m in motifs:
        if 'id' not in m or 'matrix' not in m:
            print('Missing id or matrix in motif')
            sys.exit(1)
        if len(m['matrix'][0]) != 4:
            print('Matrix row should have 4 columns')
            sys.exit(1)
    print('OK')
    sys.exit(0)
except Exception as e:
    print(f'Error: {e}')
    sys.exit(1)
" > "$WORK_DIR/json_check.txt" 2>&1; then
    pass "JSON output structure validation"
else
    fail "JSON output structure validation" "$(cat "$WORK_DIR/json_check.txt")"
fi

# JSON round-trip: meme -> json -> meme
python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test.meme" -j JASPAR2026 -f json > "$WORK_DIR/rt_json.json" 2>&1
python3 "$PYTHON/homer2meme.py" -i "$WORK_DIR/rt_json.json" -f json > "$WORK_DIR/rt_json_meme.meme" 2>&1

# Compare matrix values between original MEME and round-trip MEME
if python3 -c "
import sys
def parse_meme(path):
    motifs = []
    with open(path) as f:
        in_matrix = False
        current = None
        for line in f:
            s = line.strip()
            if s.startswith('MOTIF'):
                current = {'id': s, 'matrix': []}
                motifs.append(current)
                in_matrix = False
            elif s.startswith('letter-probability'):
                in_matrix = True
            elif s.startswith('//') or (s.startswith('MOTIF') and current):
                in_matrix = False
            elif in_matrix and s and (s[0].isdigit() or s.startswith('.')):
                current['matrix'].append([float(x) for x in s.split()])
    return motifs

orig = parse_meme('$FIXTURES/test.meme')
rt = parse_meme('$WORK_DIR/rt_json_meme.meme')
ok = True
if len(orig) != len(rt):
    print(f'Motif count mismatch: {len(orig)} vs {len(rt)}')
    ok = False
else:
    for o, r in zip(orig, rt):
        if len(o['matrix']) != len(r['matrix']):
            print(f'Matrix length mismatch for {o[\"id\"]}: {len(o[\"matrix\"])} vs {len(r[\"matrix\"])}')
            ok = False
            continue
        for i, (orow, rrow) in enumerate(zip(o['matrix'], r['matrix'])):
            for j, (ov, rv) in enumerate(zip(orow, rrow)):
                if abs(ov - rv) > 1e-5:
                    print(f'Value mismatch at {o[\"id\"]} row {i} col {j}: {ov} vs {rv}')
                    ok = False
if ok:
    print('OK')
    sys.exit(0)
else:
    sys.exit(1)
" > "$WORK_DIR/json_rt_check.txt" 2>&1; then
    pass "JSON round-trip meme->json->meme matrix consistency"
else
    fail "JSON round-trip meme->json->meme matrix consistency" "$(cat "$WORK_DIR/json_rt_check.txt")"
fi

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test.meme" -j JASPAR2026 -f json > "$WORK_DIR/rs_json.json" 2>&1
    check_diff "$WORK_DIR/py_json.json" "$WORK_DIR/rs_json.json" "Python vs Rust JSON output"
fi

# Test homer2meme JSON input
python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test.meme" -j JASPAR2026 -f json > "$WORK_DIR/json_input.json" 2>&1
python3 "$PYTHON/homer2meme.py" -i "$WORK_DIR/json_input.json" -f json > "$WORK_DIR/json_to_meme.meme" 2>&1
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test.homer" > "$WORK_DIR/direct_to_meme.meme" 2>&1

# Compare matrix values
if python3 -c "
import sys
def parse_meme(path):
    motifs = []
    with open(path) as f:
        in_matrix = False
        current = None
        for line in f:
            s = line.strip()
            if s.startswith('MOTIF'):
                current = {'id': s, 'matrix': []}
                motifs.append(current)
                in_matrix = False
            elif s.startswith('letter-probability'):
                in_matrix = True
            elif s.startswith('//') or (s.startswith('MOTIF') and current):
                in_matrix = False
            elif in_matrix and s and (s[0].isdigit() or s.startswith('.')):
                current['matrix'].append([float(x) for x in s.split()])
    return motifs

j = parse_meme('$WORK_DIR/json_to_meme.meme')
d = parse_meme('$WORK_DIR/direct_to_meme.meme')
ok = True
for a, b in zip(j, d):
    if len(a['matrix']) != len(b['matrix']):
        print(f'Matrix length mismatch')
        ok = False
        continue
    for i, (ar, br) in enumerate(zip(a['matrix'], b['matrix'])):
        for av, bv in zip(ar, br):
            if abs(av - bv) > 1e-5:
                print(f'Value mismatch at row {i}: {av} vs {bv}')
                ok = False
if ok:
    print('OK')
    sys.exit(0)
else:
    sys.exit(1)
" > "$WORK_DIR/json_input_check.txt" 2>&1; then
    pass "homer2meme JSON input produces correct MEME"
else
    fail "homer2meme JSON input produces correct MEME" "$(cat "$WORK_DIR/json_input_check.txt")"
fi

# Mixed alphabet JSON input should skip mismatched motifs
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test_mixed_alphabet.json" -f json > "$WORK_DIR/py_mixed_alphabet.meme" 2> "$WORK_DIR/py_mixed_alphabet.err"
perl "$PERL/homer2meme.pl" -i "$FIXTURES/test_mixed_alphabet.json" -f json > "$WORK_DIR/pl_mixed_alphabet.meme" 2> "$WORK_DIR/pl_mixed_alphabet.err"
check_diff "$WORK_DIR/py_mixed_alphabet.meme" "$WORK_DIR/pl_mixed_alphabet.meme" "Python vs Perl mixed-alphabet JSON"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/homer2meme" -i "$FIXTURES/test_mixed_alphabet.json" -f json > "$WORK_DIR/rs_mixed_alphabet.meme" 2> "$WORK_DIR/rs_mixed_alphabet.err"
    check_diff "$WORK_DIR/py_mixed_alphabet.meme" "$WORK_DIR/rs_mixed_alphabet.meme" "Python vs Rust mixed-alphabet JSON"
fi

check_diff "$WORK_DIR/py_mixed_alphabet.meme" "$FIXTURES/expected_mixed_alphabet.meme" "Mixed-alphabet JSON skips mismatched motifs"

if grep -q "Warning: skipping motif" "$WORK_DIR/py_mixed_alphabet.err"; then
    pass "Python warns on mixed-alphabet JSON"
else
    fail "Python warns on mixed-alphabet JSON" "$(cat "$WORK_DIR/py_mixed_alphabet.err")"
fi

if grep -q "Warning: skipping motif" "$WORK_DIR/pl_mixed_alphabet.err"; then
    pass "Perl warns on mixed-alphabet JSON"
else
    fail "Perl warns on mixed-alphabet JSON" "$(cat "$WORK_DIR/pl_mixed_alphabet.err")"
fi

if [ -n "$RUST_BIN" ]; then
    if grep -q "Warning: skipping motif" "$WORK_DIR/rs_mixed_alphabet.err"; then
        pass "Rust warns on mixed-alphabet JSON"
    else
        fail "Rust warns on mixed-alphabet JSON" "$(cat "$WORK_DIR/rs_mixed_alphabet.err")"
    fi
fi

echo ""
fi

# ---------------------------------------------------------------------------
# Test 10: --input-format explicit specification
# ---------------------------------------------------------------------------

if run_stage 10 "--input-format explicit specification"; then

# Force probability format (should skip auto-detection)
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test.homer" --input-format probability > "$WORK_DIR/force_prob.meme" 2>&1
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test.homer" > "$WORK_DIR/auto_prob.meme" 2>&1
check_diff "$WORK_DIR/force_prob.meme" "$WORK_DIR/auto_prob.meme" "Python --input-format probability vs auto"

# Force log-odds format
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test_logodds.homer" --input-format logodds > "$WORK_DIR/force_logodds.meme" 2>&1
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test_logodds.homer" > "$WORK_DIR/auto_logodds.meme" 2>&1
check_diff "$WORK_DIR/force_logodds.meme" "$WORK_DIR/auto_logodds.meme" "Python --input-format logodds vs auto"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/homer2meme" -i "$FIXTURES/test.homer" --input-format probability > "$WORK_DIR/rs_force_prob.meme" 2>&1
    check_diff "$WORK_DIR/force_prob.meme" "$WORK_DIR/rs_force_prob.meme" "Python vs Rust --input-format probability"

    "$RUST_BIN/homer2meme" -i "$FIXTURES/test_logodds.homer" --input-format logodds > "$WORK_DIR/rs_force_logodds.meme" 2>&1
    check_diff "$WORK_DIR/force_logodds.meme" "$WORK_DIR/rs_force_logodds.meme" "Python vs Rust --input-format logodds"
fi

echo ""
fi

# Test 11: --alphabet support
# ---------------------------------------------------------------------------

if run_stage 11 "--alphabet support"; then

python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test_rna.meme" --alphabet ACGU > "$WORK_DIR/py_rna.homer" 2>&1
perl "$PERL/meme2homer.pl" -i "$FIXTURES/test_rna.meme" --alphabet ACGU > "$WORK_DIR/pl_rna.homer" 2>&1
check_diff "$WORK_DIR/py_rna.homer" "$WORK_DIR/pl_rna.homer" "Python vs Perl --alphabet ACGU"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test_rna.meme" --alphabet ACGU > "$WORK_DIR/rs_rna.homer" 2>&1
    check_diff "$WORK_DIR/py_rna.homer" "$WORK_DIR/rs_rna.homer" "Python vs Rust --alphabet ACGU"
fi

python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test_rna.meme" --alphabet ACGU -f json > "$WORK_DIR/py_rna.json" 2>&1
perl "$PERL/meme2homer.pl" -i "$FIXTURES/test_rna.meme" --alphabet ACGU -f json > "$WORK_DIR/pl_rna.json" 2>&1

if python3 -c "
import json, sys
with open('$WORK_DIR/py_rna.json') as f:
    py = json.load(f)
with open('$WORK_DIR/pl_rna.json') as f:
    pl = json.load(f)
if py['motifs'][0].get('alphabet') != 'ACGU' or pl['motifs'][0].get('alphabet') != 'ACGU':
    print('Alphabet key missing or incorrect')
    sys.exit(1)
print('OK')
sys.exit(0)
" > "$WORK_DIR/json_alphabet_diff_pl.txt" 2>&1; then
    pass "Python vs Perl JSON --alphabet ACGU"
else
    fail "Python vs Perl JSON --alphabet ACGU" "$(cat "$WORK_DIR/json_alphabet_diff_pl.txt")"
fi

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test_rna.meme" --alphabet ACGU -f json > "$WORK_DIR/rs_rna.json" 2>&1
    if python3 -c "
import json, sys
with open('$WORK_DIR/py_rna.json') as f:
    py = json.load(f)
with open('$WORK_DIR/rs_rna.json') as f:
    rs = json.load(f)
if py['motifs'][0].get('alphabet') != 'ACGU' or rs['motifs'][0].get('alphabet') != 'ACGU':
    print('Alphabet key missing or incorrect')
    sys.exit(1)
print('OK')
sys.exit(0)
" > "$WORK_DIR/json_alphabet_diff_rs.txt" 2>&1; then
        pass "Python vs Rust JSON --alphabet ACGU"
    else
        fail "Python vs Rust JSON --alphabet ACGU" "$(cat "$WORK_DIR/json_alphabet_diff_rs.txt")"
    fi
fi

# meme2homer should auto-detect ALPHABET from MEME header when --alphabet is omitted
python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test_rna.meme" -f json > "$WORK_DIR/py_rna_auto.json" 2>&1
perl "$PERL/meme2homer.pl" -i "$FIXTURES/test_rna.meme" -f json > "$WORK_DIR/pl_rna_auto.json" 2>&1
check_diff "$WORK_DIR/py_rna_auto.json" "$WORK_DIR/pl_rna_auto.json" "Python vs Perl auto-detect ALPHABET"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test_rna.meme" -f json > "$WORK_DIR/rs_rna_auto.json" 2>&1
    check_diff "$WORK_DIR/py_rna_auto.json" "$WORK_DIR/rs_rna_auto.json" "Python vs Rust auto-detect ALPHABET"
fi

if python3 -c "
import json, sys
with open('$WORK_DIR/py_rna_auto.json') as f:
    py = json.load(f)
if py['motifs'][0].get('alphabet') != 'ACGU':
    print('Alphabet auto-detect missing or incorrect')
    sys.exit(1)
print('OK')
sys.exit(0)
" > "$WORK_DIR/json_alphabet_auto.txt" 2>&1; then
    pass "meme2homer auto-detects ALPHABET"
else
    fail "meme2homer auto-detects ALPHABET" "$(cat "$WORK_DIR/json_alphabet_auto.txt")"
fi

# homer2meme should respect --alphabet for MEME header
python3 "$PYTHON/homer2meme.py" -i "$WORK_DIR/py_rna.homer" --alphabet ACGU > "$WORK_DIR/py_rna_back.meme" 2>&1
perl "$PERL/homer2meme.pl" -i "$WORK_DIR/py_rna.homer" --alphabet ACGU > "$WORK_DIR/pl_rna_back.meme" 2>&1
check_diff "$WORK_DIR/py_rna_back.meme" "$WORK_DIR/pl_rna_back.meme" "Python vs Perl homer2meme --alphabet ACGU"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/homer2meme" -i "$WORK_DIR/py_rna.homer" --alphabet ACGU > "$WORK_DIR/rs_rna_back.meme" 2>&1
    check_diff "$WORK_DIR/py_rna_back.meme" "$WORK_DIR/rs_rna_back.meme" "Python vs Rust homer2meme --alphabet ACGU"
fi

if python3 -c "
import sys
lines = [l.strip() for l in open('$WORK_DIR/py_rna_back.meme')]
alphabet = next((l for l in lines if l.startswith('ALPHABET=')), '')
if alphabet != 'ALPHABET= ACGU':
    print(f'Unexpected alphabet line: {alphabet}')
    sys.exit(1)
try:
    idx = lines.index('Background letter frequencies')
    bg = lines[idx + 1]
except ValueError:
    print('Missing background letter frequencies')
    sys.exit(1)
if 'U 0.25' not in bg or 'T ' in bg:
    print(f'Unexpected background line: {bg}')
    sys.exit(1)
print('OK')
sys.exit(0)
" > "$WORK_DIR/rna_header_check.txt" 2>&1; then
    pass "homer2meme ACGU header/background"
else
    fail "homer2meme ACGU header/background" "$(cat "$WORK_DIR/rna_header_check.txt")"
fi

echo ""
fi

# ---------------------------------------------------------------------------
# Test 12: Motif Operations (--rc, --trim-edges, --min-ic)
# ---------------------------------------------------------------------------

if run_stage 12 "Motif Operations (--rc, --trim-edges, --min-ic)"; then

# Test --rc (reverse complement)
python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test.meme" -j JASPAR2026 --rc > "$WORK_DIR/py_rc.homer" 2>/dev/null
perl "$PERL/meme2homer.pl" -i "$FIXTURES/test.meme" -j JASPAR2026 --rc > "$WORK_DIR/pl_rc.homer" 2>/dev/null
check_diff "$WORK_DIR/py_rc.homer" "$WORK_DIR/pl_rc.homer" "Python vs Perl --rc"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test.meme" -j JASPAR2026 --rc > "$WORK_DIR/rs_rc.homer" 2>/dev/null
    check_diff "$WORK_DIR/py_rc.homer" "$WORK_DIR/rs_rc.homer" "Python vs Rust --rc"
fi

# Verify RC matrix is actually reversed and columns swapped
if python3 -c "
import sys
with open('$WORK_DIR/py_rc.homer') as f:
    content = f.read()
# Check RC suffix in headers
headers = [l.strip()[1:].split('\t')[0] for l in content.split('\n') if l.strip().startswith('>')]
for h in headers:
    if not h.endswith('_RC'):
        print(f'RC suffix missing in header: {h}')
        sys.exit(1)
# Check row count matches original
orig_rows = 0
with open('$FIXTURES/test.meme') as f:
    in_matrix = False
    for line in f:
        s = line.strip()
        if s.startswith('letter-probability'):
            in_matrix = True
        elif s.startswith('MOTIF') or s.startswith('//'):
            in_matrix = False
        elif in_matrix and s and (s[0].isdigit() or s.startswith('.')):
            orig_rows += 1
rc_rows = 0
in_matrix = False
for line in content.split('\n'):
    s = line.strip()
    if s.startswith('>'):
        in_matrix = True
        continue
    if s.startswith('>') == False and in_matrix and s and (s[0].isdigit() or s.startswith('.')):
        rc_rows += 1
if orig_rows != rc_rows:
    print(f'Row count mismatch: original {orig_rows} vs RC {rc_rows}')
    sys.exit(1)
print('OK')
sys.exit(0)
" > "$WORK_DIR/rc_check.txt" 2>&1; then
    pass "Reverse complement matrix correctness"
else
    fail "Reverse complement matrix correctness" "$(cat "$WORK_DIR/rc_check.txt")"
fi

# Test --trim-edges
python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test.meme" -j JASPAR2026 --trim-edges 0.5 > "$WORK_DIR/py_trim.homer" 2>/dev/null
perl "$PERL/meme2homer.pl" -i "$FIXTURES/test.meme" -j JASPAR2026 --trim-edges 0.5 > "$WORK_DIR/pl_trim.homer" 2>/dev/null
check_diff "$WORK_DIR/py_trim.homer" "$WORK_DIR/pl_trim.homer" "Python vs Perl --trim-edges"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test.meme" -j JASPAR2026 --trim-edges 0.5 > "$WORK_DIR/rs_trim.homer" 2>/dev/null
    check_diff "$WORK_DIR/py_trim.homer" "$WORK_DIR/rs_trim.homer" "Python vs Rust --trim-edges"
fi

# Test --min-ic
python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test.meme" -j JASPAR2026 --min-ic 100.0 > "$WORK_DIR/py_minic.homer" 2>/dev/null
perl "$PERL/meme2homer.pl" -i "$FIXTURES/test.meme" -j JASPAR2026 --min-ic 100.0 > "$WORK_DIR/pl_minic.homer" 2>/dev/null
check_diff "$WORK_DIR/py_minic.homer" "$WORK_DIR/pl_minic.homer" "Python vs Perl --min-ic"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test.meme" -j JASPAR2026 --min-ic 100.0 > "$WORK_DIR/rs_minic.homer" 2>/dev/null
    check_diff "$WORK_DIR/py_minic.homer" "$WORK_DIR/rs_minic.homer" "Python vs Rust --min-ic"
fi

# Verify --min-ic filters out all motifs (threshold too high)
if python3 -c "
import sys
with open('$WORK_DIR/py_minic.homer') as f:
    content = f.read()
if content.strip() != '':
    print('Expected empty output for high min-ic threshold')
    sys.exit(1)
print('OK')
sys.exit(0)
" > "$WORK_DIR/minic_check.txt" 2>&1; then
    pass "--min-ic filters motifs correctly"
else
    fail "--min-ic filters motifs correctly" "$(cat "$WORK_DIR/minic_check.txt")"
fi

# Test --trim-edges on homer2meme
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test.homer" --trim-edges 0.5 > "$WORK_DIR/py_h2m_trim.meme" 2>/dev/null
perl "$PERL/homer2meme.pl" -i "$FIXTURES/test.homer" --trim-edges 0.5 > "$WORK_DIR/pl_h2m_trim.meme" 2>/dev/null
check_diff "$WORK_DIR/py_h2m_trim.meme" "$WORK_DIR/pl_h2m_trim.meme" "Python vs Perl homer2meme --trim-edges"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/homer2meme" -i "$FIXTURES/test.homer" --trim-edges 0.5 > "$WORK_DIR/rs_h2m_trim.meme" 2>/dev/null
    check_diff "$WORK_DIR/py_h2m_trim.meme" "$WORK_DIR/rs_h2m_trim.meme" "Python vs Rust homer2meme --trim-edges"
fi

# Test --rc on homer2meme
python3 "$PYTHON/homer2meme.py" -i "$FIXTURES/test.homer" --rc > "$WORK_DIR/py_h2m_rc.meme" 2>/dev/null
perl "$PERL/homer2meme.pl" -i "$FIXTURES/test.homer" --rc > "$WORK_DIR/pl_h2m_rc.meme" 2>/dev/null
check_diff "$WORK_DIR/py_h2m_rc.meme" "$WORK_DIR/pl_h2m_rc.meme" "Python vs Perl homer2meme --rc"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/homer2meme" -i "$FIXTURES/test.homer" --rc > "$WORK_DIR/rs_h2m_rc.meme" 2>/dev/null
    check_diff "$WORK_DIR/py_h2m_rc.meme" "$WORK_DIR/rs_h2m_rc.meme" "Python vs Rust homer2meme --rc"
fi

echo ""
fi

# ---------------------------------------------------------------------------
# Test 13: MEME MOTIF word-boundary parsing
# ---------------------------------------------------------------------------

if run_stage 13 "MEME MOTIF word-boundary parsing"; then

python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test_motif_boundary.meme" -j JASPAR2026 > "$WORK_DIR/py_boundary.homer" 2> "$WORK_DIR/py_boundary.err"
perl "$PERL/meme2homer.pl" -i "$FIXTURES/test_motif_boundary.meme" -j JASPAR2026 > "$WORK_DIR/pl_boundary.homer" 2> "$WORK_DIR/pl_boundary.err"
check_diff "$WORK_DIR/py_boundary.homer" "$WORK_DIR/pl_boundary.homer" "Python vs Perl MOTIF boundary"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test_motif_boundary.meme" -j JASPAR2026 > "$WORK_DIR/rs_boundary.homer" 2> "$WORK_DIR/rs_boundary.err"
    check_diff "$WORK_DIR/py_boundary.homer" "$WORK_DIR/rs_boundary.homer" "Python vs Rust MOTIF boundary"
fi

py_count=$(grep -c '^>' "$WORK_DIR/py_boundary.homer" || true)
if [ "$py_count" -eq 1 ]; then
    pass "Python ignores MOTIFY line"
else
    fail "Python ignores MOTIFY line" "expected 1 motif, got $py_count"
fi

py_rows=$(grep -E '^[0-9.-]' "$WORK_DIR/py_boundary.homer" | wc -l | tr -d ' ')
if [ "$py_rows" -eq 2 ]; then
    pass "Python MOTIF row count"
else
    fail "Python MOTIF row count" "expected 2 rows, got $py_rows"
fi

pl_count=$(grep -c '^>' "$WORK_DIR/pl_boundary.homer" || true)
if [ "$pl_count" -eq 1 ]; then
    pass "Perl ignores MOTIFY line"
else
    fail "Perl ignores MOTIFY line" "expected 1 motif, got $pl_count"
fi

pl_rows=$(grep -E '^[0-9.-]' "$WORK_DIR/pl_boundary.homer" | wc -l | tr -d ' ')
if [ "$pl_rows" -eq 2 ]; then
    pass "Perl MOTIF row count"
else
    fail "Perl MOTIF row count" "expected 2 rows, got $pl_rows"
fi

if [ -n "$RUST_BIN" ]; then
    rs_count=$(grep -c '^>' "$WORK_DIR/rs_boundary.homer" || true)
    if [ "$rs_count" -eq 1 ]; then
        pass "Rust ignores MOTIFY line"
    else
        fail "Rust ignores MOTIFY line" "expected 1 motif, got $rs_count"
    fi

    rs_rows=$(grep -E '^[0-9.-]' "$WORK_DIR/rs_boundary.homer" | wc -l | tr -d ' ')
    if [ "$rs_rows" -eq 2 ]; then
        pass "Rust MOTIF row count"
    else
        fail "Rust MOTIF row count" "expected 2 rows, got $rs_rows"
    fi
fi

echo ""
fi

# ---------------------------------------------------------------------------
# Test 14: Negative matrix value warnings
# ---------------------------------------------------------------------------

if run_stage 14 "Negative matrix value warnings"; then

python3 "$PYTHON/meme2homer.py" -i "$FIXTURES/test_negative.meme" -j JASPAR2026 > "$WORK_DIR/py_neg.homer" 2> "$WORK_DIR/py_neg.err"
perl "$PERL/meme2homer.pl" -i "$FIXTURES/test_negative.meme" -j JASPAR2026 > "$WORK_DIR/pl_neg.homer" 2> "$WORK_DIR/pl_neg.err"
check_diff "$WORK_DIR/py_neg.homer" "$WORK_DIR/pl_neg.homer" "Python vs Perl negative values"

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$FIXTURES/test_negative.meme" -j JASPAR2026 > "$WORK_DIR/rs_neg.homer" 2> "$WORK_DIR/rs_neg.err"
    check_diff "$WORK_DIR/py_neg.homer" "$WORK_DIR/rs_neg.homer" "Python vs Rust negative values"
fi

if grep -q "Warning: negative value in matrix row" "$WORK_DIR/py_neg.err"; then
    pass "Python warns on negative values"
else
    fail "Python warns on negative values" "$(cat "$WORK_DIR/py_neg.err")"
fi

if grep -q "Warning: negative value in matrix row" "$WORK_DIR/pl_neg.err"; then
    pass "Perl warns on negative values"
else
    fail "Perl warns on negative values" "$(cat "$WORK_DIR/pl_neg.err")"
fi

if [ -n "$RUST_BIN" ]; then
    if grep -q "Warning: negative value in matrix row" "$WORK_DIR/rs_neg.err"; then
        pass "Rust warns on negative values"
    else
        fail "Rust warns on negative values" "$(cat "$WORK_DIR/rs_neg.err")"
    fi
fi

echo ""
fi

# ---------------------------------------------------------------------------
# Test 15: Version and parser metadata regressions
# ---------------------------------------------------------------------------

if run_stage 15 "Version and parser metadata regressions"; then

PYPROJECT_VERSION=$(python3 - "$SCRIPT_DIR/pyproject.toml" <<'PY'
import re
import sys
from pathlib import Path

match = re.search(r'^version\s*=\s*"([^"]+)"', Path(sys.argv[1]).read_text(), re.M)
if not match:
    raise SystemExit("missing pyproject version")
print(match.group(1))
PY
)
CARGO_VERSION=$(python3 - "$SCRIPT_DIR/rust_scripts/Cargo.toml" <<'PY'
import re
import sys
from pathlib import Path

match = re.search(r'^version\s*=\s*"([^"]+)"', Path(sys.argv[1]).read_text(), re.M)
if not match:
    raise SystemExit("missing Cargo version")
print(match.group(1))
PY
)
PERL_M2H_VERSION=$(perl -MExtUtils::MM_Unix -e "print ExtUtils::MM_Unix->parse_version('$SCRIPT_DIR/perl_scripts/meme2homer.pl')")
PERL_H2M_VERSION=$(perl -MExtUtils::MM_Unix -e "print ExtUtils::MM_Unix->parse_version('$SCRIPT_DIR/perl_scripts/homer2meme.pl')")

if python3 "$PYTHON/meme2homer.py" --version > "$WORK_DIR/py_m2h_version.txt" 2>&1 \
    && grep -Fq "$PYPROJECT_VERSION" "$WORK_DIR/py_m2h_version.txt"; then
    pass "Python meme2homer --version"
else
    fail "Python meme2homer --version" "$(cat "$WORK_DIR/py_m2h_version.txt" 2>/dev/null || true)"
fi

if python3 "$PYTHON/homer2meme.py" --version > "$WORK_DIR/py_h2m_version.txt" 2>&1 \
    && grep -Fq "$PYPROJECT_VERSION" "$WORK_DIR/py_h2m_version.txt"; then
    pass "Python homer2meme --version"
else
    fail "Python homer2meme --version" "$(cat "$WORK_DIR/py_h2m_version.txt" 2>/dev/null || true)"
fi

if perl "$PERL/meme2homer.pl" --version > "$WORK_DIR/pl_m2h_version.txt" 2>&1 \
    && grep -Fq "$PERL_M2H_VERSION" "$WORK_DIR/pl_m2h_version.txt"; then
    pass "Perl meme2homer --version"
else
    fail "Perl meme2homer --version" "$(cat "$WORK_DIR/pl_m2h_version.txt" 2>/dev/null || true)"
fi

if perl "$PERL/homer2meme.pl" --version > "$WORK_DIR/pl_h2m_version.txt" 2>&1 \
    && grep -Fq "$PERL_H2M_VERSION" "$WORK_DIR/pl_h2m_version.txt"; then
    pass "Perl homer2meme --version"
else
    fail "Perl homer2meme --version" "$(cat "$WORK_DIR/pl_h2m_version.txt" 2>/dev/null || true)"
fi

if [ -n "$PERL_M2H_VERSION" ]; then
    pass "Makefile.PL VERSION_FROM target parses"
else
    fail "Makefile.PL VERSION_FROM target parses" "$(cat "$WORK_DIR/pl_version_from.txt")"
fi

if [ -n "$RUST_BIN" ]; then
    if "$RUST_BIN/meme2homer" --version > "$WORK_DIR/rs_m2h_version.txt" 2>&1 \
        && grep -Fq "$CARGO_VERSION" "$WORK_DIR/rs_m2h_version.txt"; then
        pass "Rust meme2homer --version"
    else
        fail "Rust meme2homer --version" "$(cat "$WORK_DIR/rs_m2h_version.txt" 2>/dev/null || true)"
    fi

    if "$RUST_BIN/homer2meme" --version > "$WORK_DIR/rs_h2m_version.txt" 2>&1 \
        && grep -Fq "$CARGO_VERSION" "$WORK_DIR/rs_h2m_version.txt"; then
        pass "Rust homer2meme --version"
    else
        fail "Rust homer2meme --version" "$(cat "$WORK_DIR/rs_h2m_version.txt" 2>/dev/null || true)"
    fi
fi

cat > "$WORK_DIR/protein_bad_alength.meme" <<'EOF'
MEME version 4

ALPHABET= PROTEIN

Background letter frequencies
A 0.05 C 0.05 D 0.05 E 0.05 F 0.05 G 0.05 H 0.05 I 0.05 K 0.05 L 0.05 M 0.05 N 0.05 P 0.05 Q 0.05 R 0.05 S 0.05 T 0.05 V 0.05 W 0.05 Y 0.05

MOTIF P0001 protein_motif

letter-probability matrix: alength= 4 w= 2 nsites= 20 E= 0
  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05
  0.10  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.00  0.05

//
EOF

python3 "$PYTHON/meme2homer.py" -i "$WORK_DIR/protein_bad_alength.meme" -f json > "$WORK_DIR/py_bad_alength.json" 2> "$WORK_DIR/py_bad_alength.err"
perl "$PERL/meme2homer.pl" -i "$WORK_DIR/protein_bad_alength.meme" -f json > "$WORK_DIR/pl_bad_alength.json" 2> "$WORK_DIR/pl_bad_alength.err"
check_diff "$WORK_DIR/py_bad_alength.json" "$WORK_DIR/pl_bad_alength.json" "Python vs Perl alength conflict"

if grep -q "alength=4 conflicts with alphabet PROTEIN" "$WORK_DIR/py_bad_alength.err"; then
    pass "Python warns on alength conflict"
else
    fail "Python warns on alength conflict" "$(cat "$WORK_DIR/py_bad_alength.err")"
fi

if grep -q "alength=4 conflicts with alphabet PROTEIN" "$WORK_DIR/pl_bad_alength.err"; then
    pass "Perl warns on alength conflict"
else
    fail "Perl warns on alength conflict" "$(cat "$WORK_DIR/pl_bad_alength.err")"
fi

if grep -q '"alphabet": "PROTEIN"' "$WORK_DIR/py_bad_alength.json"; then
    pass "Python keeps alphabet-derived PROTEIN rows"
else
    fail "Python keeps alphabet-derived PROTEIN rows" "$(cat "$WORK_DIR/py_bad_alength.json")"
fi

if [ -n "$RUST_BIN" ]; then
    "$RUST_BIN/meme2homer" -i "$WORK_DIR/protein_bad_alength.meme" -f json > "$WORK_DIR/rs_bad_alength.json" 2> "$WORK_DIR/rs_bad_alength.err"
    check_diff "$WORK_DIR/py_bad_alength.json" "$WORK_DIR/rs_bad_alength.json" "Python vs Rust alength conflict"

    if grep -q "alength=4 conflicts with alphabet PROTEIN" "$WORK_DIR/rs_bad_alength.err"; then
        pass "Rust warns on alength conflict"
    else
        fail "Rust warns on alength conflict" "$(cat "$WORK_DIR/rs_bad_alength.err")"
    fi
fi

echo ""
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

echo "========================================"
echo "  Results: $PASS passed, $FAIL failed (out of $TOTAL)"
echo "========================================"

if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
exit 0
