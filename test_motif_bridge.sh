#!/usr/bin/env bash
# =============================================================================
# motif-bridge 全面测试脚本
# 测试对象：Perl / Python / Rust 三种实现，meme2homer / homer2meme 双向转换
# 用法：bash test_motif_bridge.sh [仓库根目录]
# 示例：bash test_motif_bridge.sh ~/lab/motif-bridge
# =============================================================================

set -euo pipefail

# --------------------------------------------------------------------------- #
# 路径配置
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
# 颜色输出
# --------------------------------------------------------------------------- #
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

PASS=0; FAIL=0; SKIP=0
declare -a FAILURES=()

log_section() { echo -e "\n${BOLD}${CYAN}══════════════════════════════════════════${RESET}"; echo -e "${BOLD}${CYAN}  $1${RESET}"; echo -e "${BOLD}${CYAN}══════════════════════════════════════════${RESET}"; }
log_ok()   { echo -e "  ${GREEN}✔${RESET}  $1"; ((PASS++)); }
log_fail() { echo -e "  ${RED}✘${RESET}  $1"; ((FAIL++)); FAILURES+=("$1"); }
log_skip() { echo -e "  ${YELLOW}⊘${RESET}  $1 ${YELLOW}[SKIP]${RESET}"; ((SKIP++)); }
log_info() { echo -e "  ${YELLOW}ℹ${RESET}  $1"; }

# --------------------------------------------------------------------------- #
# 辅助函数
# --------------------------------------------------------------------------- #
check_output() {
    local label="$1" file="$2" min_lines="${3:-5}"
    if [[ -s "$file" ]]; then
        local lines
        lines=$(wc -l < "$file")
        if (( lines >= min_lines )); then
            log_ok "$label → $(wc -l < "$file") 行，$(wc -c < "$file") 字节"
        else
            log_fail "$label → 输出行数过少（$lines 行，期望 >=$min_lines）"
        fi
    else
        log_fail "$label → 输出文件为空或不存在"
    fi
}

check_motif_count() {
    local label="$1" file="$2" pattern="$3" min_count="${4:-1}"
    local count
    count=$(grep -c "$pattern" "$file" 2>/dev/null || echo 0)
    if (( count >= min_count )); then
        log_ok "$label → 检测到 $count 个 motif"
    else
        log_fail "$label → motif 数量不足（$count，期望 >=$min_count）"
    fi
}

check_roundtrip() {
    local label="$1" orig="$2" roundtrip="$3"
    local orig_n rt_n
    orig_n=$(grep -c "^MOTIF" "$orig" 2>/dev/null || echo 0)
    rt_n=$(grep -c "^MOTIF" "$roundtrip" 2>/dev/null || echo 0)
    if (( orig_n > 0 && rt_n == orig_n )); then
        log_ok "$label 往返一致（MOTIF 数均为 $orig_n）"
    elif (( orig_n > 0 && rt_n > 0 )); then
        log_fail "$label 往返有损（原始 $orig_n，还原 $rt_n）"
    else
        log_fail "$label 往返检测失败（orig=$orig_n rt=$rt_n）"
    fi
}

time_cmd() {
    local label="$1"; shift
    local start end elapsed
    start=$(date +%s%3N)
    "$@" > /dev/null 2>&1
    end=$(date +%s%3N)
    elapsed=$(( end - start ))
    echo -e "    ${YELLOW}耗时: ${elapsed} ms${RESET}"
}

# --------------------------------------------------------------------------- #
# 环境检查
# --------------------------------------------------------------------------- #
log_section "0. 环境检查"

[[ -d "$DATA_DIR" ]]      && log_ok "data/ 目录存在" || { log_fail "data/ 目录不存在: $DATA_DIR"; exit 1; }
[[ -f "$HOMER_INPUT" ]]   && log_ok "homer.known.motifs 存在 ($(wc -l < "$HOMER_INPUT") 行)" || log_fail "homer.known.motifs 不存在"
[[ -f "$MEME_SMALL" ]]    && log_ok "JASPAR2024_small.meme 存在 ($(wc -l < "$MEME_SMALL") 行)" || log_fail "JASPAR2024_small.meme 不存在"
[[ -f "$MEME_LARGE" ]]    && log_ok "JASPAR2024_vertebrates.meme 存在 ($(wc -l < "$MEME_LARGE") 行)" || log_fail "JASPAR2024_vertebrates.meme 不存在"

command -v perl   &>/dev/null && log_ok "Perl:   $(perl -e 'print $]')" || log_skip "Perl 未安装"
command -v python3 &>/dev/null && log_ok "Python: $(python3 --version 2>&1)" || log_skip "Python3 未安装"

if [[ -x "$RUST_M2H" && -x "$RUST_H2M" ]]; then
    log_ok "Rust 二进制已编译: $RUST_BIN"
    RUST_AVAIL=1
else
    log_info "Rust 二进制未找到，尝试编译..."
    if command -v cargo &>/dev/null; then
        echo "    cargo build --release ..."
        (cd "$REPO_ROOT/rust_scripts" && cargo build --release -q) \
            && log_ok "Rust 编译成功" && RUST_AVAIL=1 \
            || { log_fail "Rust 编译失败"; RUST_AVAIL=0; }
    else
        log_skip "cargo 未安装，跳过 Rust 测试"
        RUST_AVAIL=0
    fi
fi

mkdir -p "$OUT_DIR"

# --------------------------------------------------------------------------- #
# 1. meme2homer — 小文件
# --------------------------------------------------------------------------- #
log_section "1. meme2homer (JASPAR2024_small.meme)"

MEME_MOTIF_N=$(grep -c "^MOTIF" "$MEME_SMALL")
log_info "输入 MOTIF 数: $MEME_MOTIF_N"

if command -v perl &>/dev/null; then
    perl "$PERL_M2H" -i "$MEME_SMALL" -j JASPAR2024 > "$OUT_DIR/perl_m2h_small.homer" 2>/dev/null
    check_output  "Perl   meme2homer small" "$OUT_DIR/perl_m2h_small.homer"
    check_motif_count "Perl   meme2homer motif 数" "$OUT_DIR/perl_m2h_small.homer" "^>" "$MEME_MOTIF_N"
fi

if command -v python3 &>/dev/null; then
    python3 "$PY_M2H" -i "$MEME_SMALL" -j JASPAR2024 > "$OUT_DIR/py_m2h_small.homer" 2>/dev/null
    check_output  "Python meme2homer small" "$OUT_DIR/py_m2h_small.homer"
    check_motif_count "Python meme2homer motif 数" "$OUT_DIR/py_m2h_small.homer" "^>" "$MEME_MOTIF_N"
fi

if (( RUST_AVAIL )); then
    "$RUST_M2H" -i "$MEME_SMALL" -j JASPAR2024 > "$OUT_DIR/rust_m2h_small.homer" 2>/dev/null
    check_output  "Rust   meme2homer small" "$OUT_DIR/rust_m2h_small.homer"
    check_motif_count "Rust   meme2homer motif 数" "$OUT_DIR/rust_m2h_small.homer" "^>" "$MEME_MOTIF_N"
fi

# --------------------------------------------------------------------------- #
# 2. homer2meme — HOMER 已知 motif 库
# --------------------------------------------------------------------------- #
log_section "2. homer2meme (homer.known.motifs)"

HOMER_MOTIF_N=$(grep -c "^>" "$HOMER_INPUT")
log_info "输入 MOTIF 数: $HOMER_MOTIF_N"

if command -v perl &>/dev/null; then
    perl "$PERL_H2M" -i "$HOMER_INPUT" > "$OUT_DIR/perl_h2m.meme" 2>/dev/null
    check_output      "Perl   homer2meme" "$OUT_DIR/perl_h2m.meme"
    check_motif_count "Perl   homer2meme motif 数" "$OUT_DIR/perl_h2m.meme" "^MOTIF" "$HOMER_MOTIF_N"
fi

if command -v python3 &>/dev/null; then
    python3 "$PY_H2M" -i "$HOMER_INPUT" > "$OUT_DIR/py_h2m.meme" 2>/dev/null
    check_output      "Python homer2meme" "$OUT_DIR/py_h2m.meme"
    check_motif_count "Python homer2meme motif 数" "$OUT_DIR/py_h2m.meme" "^MOTIF" "$HOMER_MOTIF_N"
fi

if (( RUST_AVAIL )); then
    "$RUST_H2M" -i "$HOMER_INPUT" > "$OUT_DIR/rust_h2m.meme" 2>/dev/null
    check_output      "Rust   homer2meme" "$OUT_DIR/rust_h2m.meme"
    check_motif_count "Rust   homer2meme motif 数" "$OUT_DIR/rust_h2m.meme" "^MOTIF" "$HOMER_MOTIF_N"
fi

# --------------------------------------------------------------------------- #
# 3. 三实现输出一致性对比
# --------------------------------------------------------------------------- #
log_section "3. 三实现输出一致性对比"

if [[ -s "$OUT_DIR/perl_m2h_small.homer" && -s "$OUT_DIR/py_m2h_small.homer" ]]; then
    if diff -q "$OUT_DIR/perl_m2h_small.homer" "$OUT_DIR/py_m2h_small.homer" &>/dev/null; then
        log_ok "meme2homer small: Perl == Python（逐字节一致）"
    else
        DIFF_LINES=$(diff "$OUT_DIR/perl_m2h_small.homer" "$OUT_DIR/py_m2h_small.homer" | grep "^[<>]" | wc -l)
        log_fail "meme2homer small: Perl ≠ Python（差异 $DIFF_LINES 行）"
    fi
fi

if (( RUST_AVAIL )) && [[ -s "$OUT_DIR/perl_m2h_small.homer" && -s "$OUT_DIR/rust_m2h_small.homer" ]]; then
    if diff -q "$OUT_DIR/perl_m2h_small.homer" "$OUT_DIR/rust_m2h_small.homer" &>/dev/null; then
        log_ok "meme2homer small: Perl == Rust（逐字节一致）"
    else
        DIFF_LINES=$(diff "$OUT_DIR/perl_m2h_small.homer" "$OUT_DIR/rust_m2h_small.homer" | grep "^[<>]" | wc -l)
        log_fail "meme2homer small: Perl ≠ Rust（差异 $DIFF_LINES 行）"
    fi
fi

if [[ -s "$OUT_DIR/perl_h2m.meme" && -s "$OUT_DIR/py_h2m.meme" ]]; then
    P_N=$(grep -c "^MOTIF" "$OUT_DIR/perl_h2m.meme")
    Y_N=$(grep -c "^MOTIF" "$OUT_DIR/py_h2m.meme")
    if (( P_N == Y_N )); then
        log_ok "homer2meme: Perl == Python（MOTIF 数均为 $P_N）"
    else
        log_fail "homer2meme: Perl($P_N) ≠ Python($Y_N) motif 数不一致"
    fi
fi

if (( RUST_AVAIL )) && [[ -s "$OUT_DIR/perl_h2m.meme" && -s "$OUT_DIR/rust_h2m.meme" ]]; then
    P_N=$(grep -c "^MOTIF" "$OUT_DIR/perl_h2m.meme")
    R_N=$(grep -c "^MOTIF" "$OUT_DIR/rust_h2m.meme")
    if (( P_N == R_N )); then
        log_ok "homer2meme: Perl == Rust（MOTIF 数均为 $P_N）"
    else
        log_fail "homer2meme: Perl($P_N) ≠ Rust($R_N) motif 数不一致"
    fi
fi

# --------------------------------------------------------------------------- #
# 4. 往返转换（Round-trip）
# --------------------------------------------------------------------------- #
log_section "4. 往返转换 Round-trip"

if command -v python3 &>/dev/null && [[ -s "$OUT_DIR/py_m2h_small.homer" ]]; then
    python3 "$PY_H2M" -i "$OUT_DIR/py_m2h_small.homer" > "$OUT_DIR/roundtrip_py_meme.meme" 2>/dev/null
    check_roundtrip "Python meme→homer→meme" "$MEME_SMALL" "$OUT_DIR/roundtrip_py_meme.meme"
fi

if command -v python3 &>/dev/null && [[ -s "$OUT_DIR/py_h2m.meme" ]]; then
    python3 "$PY_M2H" -i "$OUT_DIR/py_h2m.meme" > "$OUT_DIR/roundtrip_py_homer.homer" 2>/dev/null
    ORIG_N=$(grep -c "^>" "$HOMER_INPUT")
    RT_N=$(grep -c "^>" "$OUT_DIR/roundtrip_py_homer.homer" 2>/dev/null || echo 0)
    if (( ORIG_N == RT_N )); then
        log_ok "Python homer→meme→homer 往返一致（$ORIG_N 个 motif）"
    else
        log_fail "Python homer→meme→homer 往返有损（orig=$ORIG_N rt=$RT_N）"
    fi
fi

if (( RUST_AVAIL )) && [[ -s "$OUT_DIR/rust_m2h_small.homer" ]]; then
    "$RUST_H2M" -i "$OUT_DIR/rust_m2h_small.homer" > "$OUT_DIR/roundtrip_rust_meme.meme" 2>/dev/null
    check_roundtrip "Rust   meme→homer→meme" "$MEME_SMALL" "$OUT_DIR/roundtrip_rust_meme.meme"
fi

# --------------------------------------------------------------------------- #
# 5. 单 motif 提取（-e 参数）
# --------------------------------------------------------------------------- #
log_section "5. 单 motif 提取 (-e 参数)"

FIRST_MOTIF_ID=$(grep "^MOTIF" "$MEME_SMALL" | head -1 | awk '{print $2}')
log_info "提取目标 motif ID: $FIRST_MOTIF_ID"

if command -v python3 &>/dev/null; then
    python3 "$PY_M2H" -i "$MEME_SMALL" -e "$FIRST_MOTIF_ID" > "$OUT_DIR/extract_py.homer" 2>/dev/null
    N=$(grep -c "^>" "$OUT_DIR/extract_py.homer" 2>/dev/null || echo 0)
    if (( N == 1 )); then
        log_ok "Python -e 提取：正好 1 个 motif"
    else
        log_fail "Python -e 提取：期望 1 个，得到 $N 个"
    fi
fi

if (( RUST_AVAIL )); then
    "$RUST_M2H" -i "$MEME_SMALL" -e "$FIRST_MOTIF_ID" > "$OUT_DIR/extract_rust.homer" 2>/dev/null
    N=$(grep -c "^>" "$OUT_DIR/extract_rust.homer" 2>/dev/null || echo 0)
    if (( N == 1 )); then
        log_ok "Rust   -e 提取：正好 1 个 motif"
    else
        log_fail "Rust   -e 提取：期望 1 个，得到 $N 个"
    fi
fi

# --------------------------------------------------------------------------- #
# 6. stdin 管道测试
# --------------------------------------------------------------------------- #
log_section "6. stdin 管道 (cat | ... -i -)"

if command -v python3 &>/dev/null; then
    N=$(cat "$MEME_SMALL" | python3 "$PY_M2H" -i - 2>/dev/null | grep -c "^>" || echo 0)
    (( N >= MEME_MOTIF_N )) && log_ok "Python stdin meme2homer → $N motif" || log_fail "Python stdin meme2homer → $N motif（期望 $MEME_MOTIF_N）"

    N=$(cat "$HOMER_INPUT" | python3 "$PY_H2M" -i - 2>/dev/null | grep -c "^MOTIF" || echo 0)
    (( N >= HOMER_MOTIF_N )) && log_ok "Python stdin homer2meme → $N motif" || log_fail "Python stdin homer2meme → $N motif（期望 $HOMER_MOTIF_N）"
fi

if (( RUST_AVAIL )); then
    N=$(cat "$MEME_SMALL" | "$RUST_M2H" -i - 2>/dev/null | grep -c "^>" || echo 0)
    (( N >= MEME_MOTIF_N )) && log_ok "Rust   stdin meme2homer → $N motif" || log_fail "Rust   stdin meme2homer → $N motif（期望 $MEME_MOTIF_N）"
fi

# --------------------------------------------------------------------------- #
# 7. gzip 压缩输入测试
# --------------------------------------------------------------------------- #
log_section "7. gzip 压缩输入 (.gz)"

GZ_MEME="$OUT_DIR/JASPAR2024_small.meme.gz"
GZ_HOMER="$OUT_DIR/homer.known.motifs.gz"
gzip -k -f -c "$MEME_SMALL"  > "$GZ_MEME"
gzip -k -f -c "$HOMER_INPUT" > "$GZ_HOMER"

if command -v python3 &>/dev/null; then
    N=$(python3 "$PY_M2H" -i "$GZ_MEME" 2>/dev/null | grep -c "^>" || echo 0)
    (( N >= MEME_MOTIF_N )) && log_ok "Python meme2homer .gz → $N motif" || log_fail "Python meme2homer .gz → $N motif（期望 $MEME_MOTIF_N）"

    N=$(python3 "$PY_H2M" -i "$GZ_HOMER" 2>/dev/null | grep -c "^MOTIF" || echo 0)
    (( N >= HOMER_MOTIF_N )) && log_ok "Python homer2meme .gz → $N motif" || log_fail "Python homer2meme .gz → $N motif（期望 $HOMER_MOTIF_N）"
fi

if (( RUST_AVAIL )); then
    N=$("$RUST_M2H" -i "$GZ_MEME" 2>/dev/null | grep -c "^>" || echo 0)
    (( N >= MEME_MOTIF_N )) && log_ok "Rust   meme2homer .gz → $N motif" || log_fail "Rust   meme2homer .gz → $N motif（期望 $MEME_MOTIF_N）"

    N=$("$RUST_H2M" -i "$GZ_HOMER" 2>/dev/null | grep -c "^MOTIF" || echo 0)
    (( N >= HOMER_MOTIF_N )) && log_ok "Rust   homer2meme .gz → $N motif" || log_fail "Rust   homer2meme .gz → $N motif（期望 $HOMER_MOTIF_N）"
fi

# --------------------------------------------------------------------------- #
# 8. 大文件性能测试
# --------------------------------------------------------------------------- #
log_section "8. 大文件性能 (JASPAR2024_vertebrates.meme)"

LARGE_N=$(grep -c "^MOTIF" "$MEME_LARGE")
log_info "大文件 MOTIF 数: $LARGE_N，文件大小: $(du -h "$MEME_LARGE" | cut -f1)"

if command -v python3 &>/dev/null; then
    log_info "Python meme2homer 大文件..."
    time_cmd "python meme2homer large" python3 "$PY_M2H" -i "$MEME_LARGE" -j JASPAR2024
    python3 "$PY_M2H" -i "$MEME_LARGE" -j JASPAR2024 > "$OUT_DIR/py_m2h_large.homer" 2>/dev/null
    check_motif_count "Python meme2homer large" "$OUT_DIR/py_m2h_large.homer" "^>" "$LARGE_N"
fi

if (( RUST_AVAIL )); then
    log_info "Rust   meme2homer 大文件..."
    time_cmd "rust meme2homer large" "$RUST_M2H" -i "$MEME_LARGE" -j JASPAR2024
    "$RUST_M2H" -i "$MEME_LARGE" -j JASPAR2024 > "$OUT_DIR/rust_m2h_large.homer" 2>/dev/null
    check_motif_count "Rust   meme2homer large" "$OUT_DIR/rust_m2h_large.homer" "^>" "$LARGE_N"

    if [[ -s "$OUT_DIR/py_m2h_large.homer" && -s "$OUT_DIR/rust_m2h_large.homer" ]]; then
        PL=$(grep -c "^>" "$OUT_DIR/py_m2h_large.homer")
        RL=$(grep -c "^>" "$OUT_DIR/rust_m2h_large.homer")
        (( PL == RL )) && log_ok "大文件: Python($PL) == Rust($RL)" || log_fail "大文件: Python($PL) ≠ Rust($RL)"
    fi
fi

# --------------------------------------------------------------------------- #
# 9. 输出格式合规性检查
# --------------------------------------------------------------------------- #
log_section "9. 输出格式合规性"

if [[ -s "$OUT_DIR/py_h2m.meme" ]]; then
    grep -q "^MEME version" "$OUT_DIR/py_h2m.meme"     && log_ok "MEME header 存在" || log_fail "MEME header 缺失"
    grep -q "^ALPHABET= ACGT" "$OUT_DIR/py_h2m.meme"   && log_ok "ALPHABET 字段存在" || log_fail "ALPHABET 字段缺失"
    grep -q "letter-probability matrix" "$OUT_DIR/py_h2m.meme" && log_ok "letter-probability matrix 字段存在" || log_fail "letter-probability matrix 字段缺失"
fi

if [[ -s "$OUT_DIR/py_m2h_small.homer" ]]; then
    BAD=$(grep "^>" "$OUT_DIR/py_m2h_small.homer" | awk -F'\t' 'NF!=6' | wc -l)
    (( BAD == 0 )) && log_ok "HOMER header 全部为 6 列 tab 分隔" || log_fail "HOMER header 有 $BAD 行列数不等于 6"
fi

if [[ -s "$OUT_DIR/py_h2m.meme" ]]; then
    BAD=$(awk '/^  [0-9]/{s=$1+$2+$3+$4; if(s<0.98||s>1.02) print NR": "s}' "$OUT_DIR/py_h2m.meme" | wc -l)
    (( BAD == 0 )) && log_ok "MEME 矩阵行求和全部在 [0.98, 1.02]" || log_fail "MEME 矩阵有 $BAD 行概率和异常"
fi

# --------------------------------------------------------------------------- #
# 汇总报告
# --------------------------------------------------------------------------- #
TOTAL=$(( PASS + FAIL + SKIP ))
echo ""
echo -e "${BOLD}══════════════════════════════════════════${RESET}"
echo -e "${BOLD}  测试结果汇总${RESET}"
echo -e "${BOLD}══════════════════════════════════════════${RESET}"
echo -e "  总计:  $TOTAL"
echo -e "  ${GREEN}通过:  $PASS${RESET}"
echo -e "  ${RED}失败:  $FAIL${RESET}"
echo -e "  ${YELLOW}跳过:  $SKIP${RESET}"
echo -e "  输出目录: $OUT_DIR"

if (( ${#FAILURES[@]} > 0 )); then
    echo ""
    echo -e "${RED}${BOLD}  失败项列表:${RESET}"
    for f in "${FAILURES[@]}"; do
        echo -e "  ${RED}✘${RESET}  $f"
    done
fi

echo ""
if (( FAIL == 0 )); then
    echo -e "${GREEN}${BOLD}  🎉 所有测试通过！${RESET}"
    exit 0
else
    echo -e "${RED}${BOLD}  ❌ 有 $FAIL 项测试失败，请检查上方输出${RESET}"
    exit 1
fi
