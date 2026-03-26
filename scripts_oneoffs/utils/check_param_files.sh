#!/bin/bash
# check_param_files.sh
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
# bash scripts/check_param_files.sh 2>&1 | tee param_report.txt

SEP="════════════════════════════════════════════════════════════════"

show() {
    echo "$SEP"
    echo "$1"
    echo "$SEP"
    if [ -f "$2" ]; then
        echo "FOUND: $2  ($(wc -l < "$2") lines)"
        head -20 "$2"
    else
        echo "NOT FOUND: $2"
    fi
    echo ""
}

# ── ZINC2 ─────────────────────────────────────────────────────────────────────
show "ZINC2 — haplotype parameters" "helpfiles/ZINC2/ZINC2_haplotype_parameters.R"
show "ZINC2 — test file (females)"  "helpfiles/ZINC2/Zinc2.test.F.txt"
show "ZINC2 — test file (males)"    "helpfiles/ZINC2/Zinc2.test.M.txt"

# ── AGING ─────────────────────────────────────────────────────────────────────
show "AGING — haplotype parameters" "helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R"
show "AGING — test file (males day 20)" "helpfiles/AGE_SY/AGE_SY20_M.test.txt"

# ── PUPATION ─────────────────────────────────────────────────────────────────
echo "$SEP"
echo "PUPATION — all files in Sarah's helpfiles/"
echo "$SEP"
ls /dfs7/adl/sruckman/XQTL/XQTL2/helpfiles/ 2>/dev/null || echo "NOT FOUND"
echo ""

for f in /dfs7/adl/sruckman/XQTL/XQTL2/helpfiles/*; do
    [ -f "$f" ] || continue
    echo "$SEP"
    echo "PUPATION — $(basename $f)  ($(wc -l < "$f") lines)"
    echo "$SEP"
    head -20 "$f"
    echo ""
done

# ── MALATHION ────────────────────────────────────────────────────────────────
echo "$SEP"
echo "MALATHION — all files in helperfiles/"
echo "$SEP"
ls /dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/helperfiles/ 2>/dev/null || echo "NOT FOUND"
echo ""

for f in /dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/helperfiles/*; do
    [ -f "$f" ] || continue
    echo "$SEP"
    echo "MALATHION — $(basename $f)  ($(wc -l < "$f") lines)"
    echo "$SEP"
    head -20 "$f"
    echo ""
done

# ── GENERIC REFERENCE ────────────────────────────────────────────────────────
show "GENERIC haplotype parameters"           "helpfiles/generic_haplotype_parameters.R"
show "GENERIC haplotype parameters (A panel)" "helpfiles/A_generic_haplotype_parameters.R"
