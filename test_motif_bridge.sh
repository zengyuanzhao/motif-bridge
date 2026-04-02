#!/usr/bin/env bash
# =============================================================================
# motif-bridge Comprehensive Test Script
# Testing: Perl / Python / Rust implementations, meme2homer / homer2meme conversion
# Usage: bash test_motif_bridge.sh [repo_root]
# Example: bash test_motif_bridge.sh ~/lab/motif-bridge
# =============================================================================

set -euo pipefail

# --------------------------------------------------------------------------- #
# Path Configuration
# --------------------------------------------------------------------------- #
REPO_ROOT="${1:-$(pwd)}"
DATA_DIR="$REPO_ROOT/data"
OUT_DIR="$REPO_ROOT/test_output"
RUST_BIN="$REPO_ROOT/rust_scripts/target/release"

HOMER_INPUT="$DATA_DIR/homer.known.motifs"
MEME_SMALL="$DATA_DIR/JASPAR2024_small.meme"
MEME_LARGE="$DATA_DIR/JASPAR2024_vertebrates.meme"

PERL_M2H="$REPO_ROOT/perl_scripts/meme2homer.pl"
PERL_H2M="$REPO_ROOT/perl_scripts/homer2meme.pl"
PY_M2H="$REPO_ROOT/python_scripts/meme2homer.py"
PY_H2M="$REPO_ROOT/python_scripts/homer2meme.py"
RUST_M2H="$RUST_BIN/meme2homer"
RUST_H2M="$RUST_BIN/homer2meme"

# --------------------------------------------------------------------------- #
# Color Output Functions
# --------------------------------------------------------------------------- #
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

PASS=0; FAIL=0; SKIP=0
declare -a FAILURES=()

log_section() { echo -e "\n${BOLD}${CYAN}══════════════════════════════════════════${RESET}"; echo -e "${BOLD}${CYAN}  $1${RESET}"; echo -e "${BOLD}${CYAN}══════════════════════════════════════════${RESET}"; }
log_ok()   { echo -e "  ${GREEN}✔${RESET}  $1"; ((++PASS)); }
log_fail() { echo -e "  ${RED}✘${RESET}  $1"; ((++FAIL)); FAILURES+=("$1"); }
log_skip() { echo -e "  ${YELLOW}⊘${RESET}  $1 ${YELLOW}[SKIP]${RESET}"; ((++SKIP)); }
log_info() { echo -e "  ${YELLOW}ℹ${RESET}  $1"; }

# --------------------------------------------------------------------------- #
# Helper Functions
# --------------------------------------------------------------------------- #
check_output() {
    local label="$1" file="$2" min_lines="${3:-5}"
    if [[ -s "$file" ]]; then
        local lines
        lines=$(wc -l < "$file")
        if (( lines >= min_lines )); then
            log_ok "$label → $(wc -l < "$file") lines, $(wc -c < "$file") bytes"
        else
            log_fail "$label → output has too few lines ($lines lines, expected >=$min_lines)"
        fi
    else
        log_fail "$label → output file is empty or does not exist"
    fi
}

check_motif_count() {
    local label="$1" file="$2" pattern="$3" min_count="${4:-1}"
    local count
    count=$(grep -c "$pattern" "$file" 2>/dev/null || echo 0)
    if (( count >= min_count )); then
        log_ok "$label → detected $count motifs"
    else
        log_fail "$label → motif count insufficient ($count, expected >=$min_count)"
    fi
}

check_roundtrip() {
    local label="$1" orig="$2" roundtrip="$3"
    local orig_n rt_n
    orig_n=$(grep -c "^MOTIF" "$orig" 2>/dev/null || echo 0)
    rt_n=$(grep -c "^MOTIF" "$roundtrip" 2>/dev/null || echo 0)
    if (( orig_n > 0 && rt_n == orig_n )); then
        log_ok "$label round-trip consistent (MOTIF count: $orig_n)"
    elif (( orig_n > 0 && rt_n > 0 )); then
        log_fail "$label round-trip lossy (original $orig_n, restored $rt_n)"
    else
        log_fail "$label round-trip check failed (orig=$orig_n rt=$rt_n)"
    fi
}

time_cmd() {
    local label="$1"; shift
    local start end elapsed
    start=$(date +%s%3N)
    "$@" > /dev/null 2>&1
    end=$(date +%s%3N)
    elapsed=$(( end - start ))
    echo -e "    ${YELLOW}Time: ${elapsed} ms${RESET}"
}

# --------------------------------------------------------------------------- #
# Environment Check
# --------------------------------------------------------------------------- #
log_section "0. Environment Check"

[[ -d "$DATA_DIR" ]]      && log_ok "data/ directory exists" || { log_fail "data/ directory not found: $DATA_DIR"; exit 1; }
[[ -f "$HOMER_INPUT" ]]   && log_ok "homer.known.motifs exists ($(wc -l < "$HOMER_INPUT") lines)" || log_fail "homer.known.motifs not found"
[[ -f "$MEME_SMALL" ]]    && log_ok "JASPAR2024_small.meme exists ($(wc -l < "$MEME_SMALL") lines)" || log_fail "JASPAR2024_small.meme not found"
[[ -f "$MEME_LARGE" ]]    && log_ok "JASPAR2024_vertebrates.meme exists ($(wc -l < "$MEME_LARGE") lines)" || log_fail "JASPAR2024_vertebrates.meme not found"

command -v perl   &>/dev/null && log_ok "Perl:   $(perl -e 'print $]')" || log_skip "Perl not installed"
command -v python3 &>/dev/null && log_ok "Python: $(python3 --version 2>&1)" || log_skip "Python3 not installed"

if [[ -x "$RUST_M2H" && -x "$RUST_H2M" ]]; then
    log_ok "Rust binaries compiled: $RUST_BIN"
    RUST_AVAIL=1
else
    log_info "Rust binaries not found, attempting to compile..."
    if command -v cargo &>/dev/null; then
        echo "    cargo build --release ..."
        (cd "$REPO_ROOT/rust_scripts" && cargo build --release -q) \
            && log_ok "Rust compiled successfully" && RUST_AVAIL=1 \
            || { log_fail "Rust compilation failed"; RUST_AVAIL=0; }
    else
        log_skip "cargo not installed, skipping Rust tests"
        RUST_AVAIL=0
    fi
fi

mkdir -p "$OUT_DIR"

# --------------------------------------------------------------------------- #
# 1. meme2homer — Small File
# --------------------------------------------------------------------------- #
log_section "1. meme2homer (JASPAR2024_small.meme)"

MEME_MOTIF_N=$(grep -c "^MOTIF" "$MEME_SMALL")
log_info "Input MOTIF count: $MEME_MOTIF_N"

if command -v perl &>/dev/null; then
    perl "$PERL_M2H" -i "$MEME_SMALL" -j JASPAR2024 > "$OUT_DIR/perl_m2h_small.homer" 2>/dev/null
    check_output  "Perl   meme2homer small" "$OUT_DIR/perl_m2h_small.homer"
    check_motif_count "Perl   meme2homer motif count" "$OUT_DIR/perl_m2h_small.homer" "^>" "$MEME_MOTIF_N"
fi

if command -v python3 &>/dev/null; then
    python3 "$PY_M2H" -i "$MEME_SMALL" -j JASPAR2024 > "$OUT_DIR/py_m2h_small.homer" 2>/dev/null
    check_output  "Python meme2homer small" "$OUT_DIR/py_m2h_small.homer"
    check_motif_count "Python meme2homer motif count" "$OUT_DIR/py_m2h_small.homer" "^>" "$MEME_MOTIF_N"
fi

if (( RUST_AVAIL )); then
    "$RUST_M2H" -i "$MEME_SMALL" -j JASPAR2024 > "$OUT_DIR/rust_m2h_small.homer" 2>/dev/null
    check_output  "Rust   meme2homer small" "$OUT_DIR/rust_m2h_small.homer"
    check_motif_count "Rust   meme2homer motif count" "$OUT_DIR/rust_m2h_small.homer" "^>" "$MEME_MOTIF_N"
fi

# --------------------------------------------------------------------------- #
# 2. homer2meme — HOMER Known Motif Library
# --------------------------------------------------------------------------- #
log_section "2. homer2meme (homer.known.motifs)"

HOMER_MOTIF_N=$(grep -c "^>" "$HOMER_INPUT")
log_info "Input MOTIF count: $HOMER_MOTIF_N"

if command -v perl &>/dev/null; then
    perl "$PERL_H2M" -i "$HOMER_INPUT" > "$OUT_DIR/perl_h2m.meme" 2>/dev/null
    check_output      "Perl   homer2meme" "$OUT_DIR/perl_h2m.meme"
    check_motif_count "Perl   homer2meme motif count" "$OUT_DIR/perl_h2m.meme" "^MOTIF" "$HOMER_MOTIF_N"
fi

if command -v python3 &>/dev/null; then
    python3 "$PY_H2M" -i "$HOMER_INPUT" > "$OUT_DIR/py_h2m.meme" 2>/dev/null
    check_output      "Python homer2meme" "$OUT_DIR/py_h2m.meme"
    check_motif_count "Python homer2meme motif count" "$OUT_DIR/py_h2m.meme" "^MOTIF" "$HOMER_MOTIF_N"
fi

if (( RUST_AVAIL )); then
    "$RUST_H2M" -i "$HOMER_INPUT" > "$OUT_DIR/rust_h2m.meme" 2>/dev/null
    check_output      "Rust   homer2meme" "$OUT_DIR/rust_h2m.meme"
    check_motif_count "Rust   homer2meme motif count" "$OUT_DIR/rust_h2m.meme" "^MOTIF" "$HOMER_MOTIF_N"
fi

# --------------------------------------------------------------------------- #
# 3. Cross-Implementation Output Consistency
# --------------------------------------------------------------------------- #
log_section "3. Cross-Implementation Consistency"

if [[ -s "$OUT_DIR/perl_m2h_small.homer" && -s "$OUT_DIR/py_m2h_small.homer" ]]; then
    if diff -q "$OUT_DIR/perl_m2h_small.homer" "$OUT_DIR/py_m2h_small.homer" &>/dev/null; then
        log_ok "meme2homer small: Perl == Python (byte-identical)"
    else
        DIFF_LINES=$(diff "$OUT_DIR/perl_m2h_small.homer" "$OUT_DIR/py_m2h_small.homer" | grep "^[<>]" | wc -l)
        log_fail "meme2homer small: Perl ≠ Python ($DIFF_LINES lines differ)"
    fi
fi

if (( RUST_AVAIL )) && [[ -s "$OUT_DIR/perl_m2h_small.homer" && -s "$OUT_DIR/rust_m2h_small.homer" ]]; then
    if diff -q "$OUT_DIR/perl_m2h_small.homer" "$OUT_DIR/rust_m2h_small.homer" &>/dev/null; then
        log_ok "meme2homer small: Perl == Rust (byte-identical)"
    else
        DIFF_LINES=$(diff "$OUT_DIR/perl_m2h_small.homer" "$OUT_DIR/rust_m2h_small.homer" | grep "^[<>]" | wc -l)
        log_fail "meme2homer small: Perl ≠ Rust ($DIFF_LINES lines differ)"
    fi
fi

if [[ -s "$OUT_DIR/perl_h2m.meme" && -s "$OUT_DIR/py_h2m.meme" ]]; then
    P_N=$(grep -c "^MOTIF" "$OUT_DIR/perl_h2m.meme")
    Y_N=$(grep -c "^MOTIF" "$OUT_DIR/py_h2m.meme")
    if (( P_N == Y_N )); then
        log_ok "homer2meme: Perl == Python (MOTIF count: $P_N)"
    else
        log_fail "homer2meme: Perl($P_N) ≠ Python($Y_N) motif count mismatch"
    fi
fi

if (( RUST_AVAIL )) && [[ -s "$OUT_DIR/perl_h2m.meme" && -s "$OUT_DIR/rust_h2m.meme" ]]; then
    P_N=$(grep -c "^MOTIF" "$OUT_DIR/perl_h2m.meme")
    R_N=$(grep -c "^MOTIF" "$OUT_DIR/rust_h2m.meme")
    if (( P_N == R_N )); then
        log_ok "homer2meme: Perl == Rust (MOTIF count: $P_N)"
    else
        log_fail "homer2meme: Perl($P_N) ≠ Rust($R_N) motif count mismatch"
    fi
fi

# --------------------------------------------------------------------------- #
# 4. Round-trip Conversion
# --------------------------------------------------------------------------- #
log_section "4. Round-trip Conversion"

if command -v python3 &>/dev/null && [[ -s "$OUT_DIR/py_m2h_small.homer" ]]; then
    python3 "$PY_H2M" -i "$OUT_DIR/py_m2h_small.homer" > "$OUT_DIR/roundtrip_py_meme.meme" 2>/dev/null
    check_roundtrip "Python meme→homer→meme" "$MEME_SMALL" "$OUT_DIR/roundtrip_py_meme.meme"
fi

if command -v python3 &>/dev/null && [[ -s "$OUT_DIR/py_h2m.meme" ]]; then
    python3 "$PY_M2H" -i "$OUT_DIR/py_h2m.meme" > "$OUT_DIR/roundtrip_py_homer.homer" 2>/dev/null
    ORIG_N=$(grep -c "^>" "$HOMER_INPUT")
    RT_N=$(grep -c "^>" "$OUT_DIR/roundtrip_py_homer.homer" 2>/dev/null || echo 0)
    if (( ORIG_N == RT_N )); then
        log_ok "Python homer→meme→homer round-trip consistent ($ORIG_N motifs)"
    else
        log_fail "Python homer→meme→homer round-trip lossy (orig=$ORIG_N rt=$RT_N)"
    fi
fi

if (( RUST_AVAIL )) && [[ -s "$OUT_DIR/rust_m2h_small.homer" ]]; then
    "$RUST_H2M" -i "$OUT_DIR/rust_m2h_small.homer" > "$OUT_DIR/roundtrip_rust_meme.meme" 2>/dev/null
    check_roundtrip "Rust   meme→homer→meme" "$MEME_SMALL" "$OUT_DIR/roundtrip_rust_meme.meme"
fi

# --------------------------------------------------------------------------- #
# 5. Single Motif Extraction (-e option)
# --------------------------------------------------------------------------- #
log_section "5. Single Motif Extraction (-e option)"

FIRST_MOTIF_ID=$(grep "^MOTIF" "$MEME_SMALL" | head -1 | awk '{print $2}')
log_info "Target motif ID for extraction: $FIRST_MOTIF_ID"

if command -v python3 &>/dev/null; then
    python3 "$PY_M2H" -i "$MEME_SMALL" -e "$FIRST_MOTIF_ID" > "$OUT_DIR/extract_py.homer" 2>/dev/null
    N=$(grep -c "^>" "$OUT_DIR/extract_py.homer" 2>/dev/null || echo 0)
    if (( N == 1 )); then
        log_ok "Python -e extraction: exactly 1 motif"
    else
        log_fail "Python -e extraction: expected 1, got $N"
    fi
fi

if (( RUST_AVAIL )); then
    "$RUST_M2H" -i "$MEME_SMALL" -e "$FIRST_MOTIF_ID" > "$OUT_DIR/extract_rust.homer" 2>/dev/null
    N=$(grep -c "^>" "$OUT_DIR/extract_rust.homer" 2>/dev/null || echo 0)
    if (( N == 1 )); then
        log_ok "Rust   -e extraction: exactly 1 motif"
    else
        log_fail "Rust   -e extraction: expected 1, got $N"
    fi
fi

# --------------------------------------------------------------------------- #
# 6. stdin Pipe Test
# --------------------------------------------------------------------------- #
log_section "6. stdin Pipe (cat | ... -i -)"

if command -v python3 &>/dev/null; then
    N=$(cat "$MEME_SMALL" | python3 "$PY_M2H" -i - 2>/dev/null | grep -c "^>" || echo 0)
    (( N >= MEME_MOTIF_N )) && log_ok "Python stdin meme2homer → $N motifs" || log_fail "Python stdin meme2homer → $N motifs (expected $MEME_MOTIF_N)"

    N=$(cat "$HOMER_INPUT" | python3 "$PY_H2M" -i - 2>/dev/null | grep -c "^MOTIF" || echo 0)
    (( N >= HOMER_MOTIF_N )) && log_ok "Python stdin homer2meme → $N motifs" || log_fail "Python stdin homer2meme → $N motifs (expected $HOMER_MOTIF_N)"
fi

if (( RUST_AVAIL )); then
    N=$(cat "$MEME_SMALL" | "$RUST_M2H" -i - 2>/dev/null | grep -c "^>" || echo 0)
    (( N >= MEME_MOTIF_N )) && log_ok "Rust   stdin meme2homer → $N motifs" || log_fail "Rust   stdin meme2homer → $N motifs (expected $MEME_MOTIF_N)"
fi

# --------------------------------------------------------------------------- #
# 7. Gzip Compressed Input Test
# --------------------------------------------------------------------------- #
log_section "7. Gzip Compressed Input (.gz)"

GZ_MEME="$OUT_DIR/JASPAR2024_small.meme.gz"
GZ_HOMER="$OUT_DIR/homer.known.motifs.gz"
gzip -k -f -c "$MEME_SMALL"  > "$GZ_MEME"
gzip -k -f -c "$HOMER_INPUT" > "$GZ_HOMER"

if command -v python3 &>/dev/null; then
    N=$(python3 "$PY_M2H" -i "$GZ_MEME" 2>/dev/null | grep -c "^>" || echo 0)
    (( N >= MEME_MOTIF_N )) && log_ok "Python meme2homer .gz → $N motifs" || log_fail "Python meme2homer .gz → $N motifs (expected $MEME_MOTIF_N)"

    N=$(python3 "$PY_H2M" -i "$GZ_HOMER" 2>/dev/null | grep -c "^MOTIF" || echo 0)
    (( N >= HOMER_MOTIF_N )) && log_ok "Python homer2meme .gz → $N motifs" || log_fail "Python homer2meme .gz → $N motifs (expected $HOMER_MOTIF_N)"
fi

if (( RUST_AVAIL )); then
    N=$("$RUST_M2H" -i "$GZ_MEME" 2>/dev/null | grep -c "^>" || echo 0)
    (( N >= MEME_MOTIF_N )) && log_ok "Rust   meme2homer .gz → $N motifs" || log_fail "Rust   meme2homer .gz → $N motifs (expected $MEME_MOTIF_N)"

    N=$("$RUST_H2M" -i "$GZ_HOMER" 2>/dev/null | grep -c "^MOTIF" || echo 0)
    (( N >= HOMER_MOTIF_N )) && log_ok "Rust   homer2meme .gz → $N motifs" || log_fail "Rust   homer2meme .gz → $N motifs (expected $HOMER_MOTIF_N)"
fi

# --------------------------------------------------------------------------- #
# 8. Large File Performance Test
# --------------------------------------------------------------------------- #
log_section "8. Large File Performance (JASPAR2024_vertebrates.meme)"

LARGE_N=$(grep -c "^MOTIF" "$MEME_LARGE")
log_info "Large file MOTIF count: $LARGE_N, file size: $(du -h "$MEME_LARGE" | cut -f1)"

if command -v perl &>/dev/null; then
    log_info "Perl   meme2homer large file..."
    time_cmd "perl meme2homer large" perl "$PERL_M2H" -i "$MEME_LARGE" -j JASPAR2024
    perl "$PERL_M2H" -i "$MEME_LARGE" -j JASPAR2024 > "$OUT_DIR/perl_m2h_large.homer" 2>/dev/null
    check_motif_count "Perl   meme2homer large" "$OUT_DIR/perl_m2h_large.homer" "^>" "$LARGE_N"
fi

if command -v python3 &>/dev/null; then
    log_info "Python meme2homer large file..."
    time_cmd "python meme2homer large" python3 "$PY_M2H" -i "$MEME_LARGE" -j JASPAR2024
    python3 "$PY_M2H" -i "$MEME_LARGE" -j JASPAR2024 > "$OUT_DIR/py_m2h_large.homer" 2>/dev/null
    check_motif_count "Python meme2homer large" "$OUT_DIR/py_m2h_large.homer" "^>" "$LARGE_N"
fi

if (( RUST_AVAIL )); then
    log_info "Rust   meme2homer large file..."
    time_cmd "rust meme2homer large" "$RUST_M2H" -i "$MEME_LARGE" -j JASPAR2024
    "$RUST_M2H" -i "$MEME_LARGE" -j JASPAR2024 > "$OUT_DIR/rust_m2h_large.homer" 2>/dev/null
    check_motif_count "Rust   meme2homer large" "$OUT_DIR/rust_m2h_large.homer" "^>" "$LARGE_N"

    if [[ -s "$OUT_DIR/perl_m2h_large.homer" && -s "$OUT_DIR/py_m2h_large.homer" && -s "$OUT_DIR/rust_m2h_large.homer" ]]; then
        PL=$(grep -c "^>" "$OUT_DIR/perl_m2h_large.homer")
        YL=$(grep -c "^>" "$OUT_DIR/py_m2h_large.homer")
        RL=$(grep -c "^>" "$OUT_DIR/rust_m2h_large.homer")
        if (( PL == YL && YL == RL )); then
            log_ok "Large file: Perl($PL) == Python($YL) == Rust($RL)"
        else
            log_fail "Large file: results inconsistent Perl($PL) Python($YL) Rust($RL)"
        fi
    fi
fi

# --------------------------------------------------------------------------- #
# 9. Output Format Compliance Check
# --------------------------------------------------------------------------- #
log_section "9. Output Format Compliance"

if [[ -s "$OUT_DIR/py_h2m.meme" ]]; then
    grep -q "^MEME version" "$OUT_DIR/py_h2m.meme"     && log_ok "MEME header present" || log_fail "MEME header missing"
    grep -q "^ALPHABET= ACGT" "$OUT_DIR/py_h2m.meme"   && log_ok "ALPHABET field present" || log_fail "ALPHABET field missing"
    grep -q "letter-probability matrix" "$OUT_DIR/py_h2m.meme" && log_ok "letter-probability matrix field present" || log_fail "letter-probability matrix field missing"
fi

if [[ -s "$OUT_DIR/py_m2h_small.homer" ]]; then
    BAD=$(grep "^>" "$OUT_DIR/py_m2h_small.homer" | awk -F'\t' 'NF!=6' | wc -l)
    (( BAD == 0 )) && log_ok "HOMER headers are all 6-column tab-separated" || log_fail "HOMER headers have $BAD lines with incorrect column count"
fi

if [[ -s "$OUT_DIR/py_h2m.meme" ]]; then
    BAD=$(awk '/^  [0-9]/{s=$1+$2+$3+$4; if(s<0.98||s>1.02) print NR": "s}' "$OUT_DIR/py_h2m.meme" | wc -l)
    (( BAD == 0 )) && log_ok "MEME matrix row sums all within [0.98, 1.02]" || log_fail "MEME matrix has $BAD rows with abnormal probability sum"
fi

# --------------------------------------------------------------------------- #
# Summary Report
# --------------------------------------------------------------------------- #
TOTAL=$(( PASS + FAIL + SKIP ))
echo ""
echo -e "${BOLD}══════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Test Summary${RESET}"
echo -e "${BOLD}══════════════════════════════════════════${RESET}"
echo -e "  Total:  $TOTAL"
echo -e "  ${GREEN}Pass:   $PASS${RESET}"
echo -e "  ${RED}Fail:   $FAIL${RESET}"
echo -e "  ${YELLOW}Skip:   $SKIP${RESET}"
echo -e "  Output: $OUT_DIR"

if (( ${#FAILURES[@]} > 0 )); then
    echo ""
    echo -e "${RED}${BOLD}  Failures:${RESET}"
    for f in "${FAILURES[@]}"; do
        echo -e "  ${RED}✘${RESET}  $f"
    done
fi

echo ""
if (( FAIL == 0 )); then
    echo -e "${GREEN}${BOLD} All tests passed!${RESET}"
    exit 0
else
    echo -e "${RED}${BOLD} $FAIL test(s) failed, check output above${RESET}"
    exit 1
fi
