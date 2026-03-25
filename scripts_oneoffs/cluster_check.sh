#!/bin/bash
###############################################################################
# cluster_check.sh  (v4 — safe, no large file dumps)
#
# Run from /dfs7/adl/tdlong/fly_pool/XQTL2:
#   bash scripts/cluster_check.sh 2>&1 | tee cluster_report.txt
# Then: scp cluster_report.txt back and paste to Claude.
###############################################################################

XQTL2="/dfs7/adl/tdlong/fly_pool/XQTL2"
SARAH="/dfs7/adl/sruckman/XQTL/XQTL2"
MAL_BASE="/dfs7/adl/tdlong/fly_pool/newpipeline_Nov23"
SEP="────────────────────────────────────────────────────────────────"

safe_head() {
    # safe_head FILE NLINES LABEL
    local f=$1 n=${2:-3} lbl=${3:-}
    [ -n "$lbl" ] && echo "  [$lbl]"
    if [ -f "$f" ]; then
        local nl; nl=$(wc -l < "$f")
        echo "  file: $f  ($nl lines)"
        head -"$n" "$f"
    else
        echo "  (not found: $f)"
    fi
}

echo "###############################################################################"
echo "# CLUSTER DIAGNOSTIC REPORT v4"
echo "# Host: $(hostname)   Date: $(date)"
echo "###############################################################################"
echo ""

# ── AGING ─────────────────────────────────────────────────────────────────────
echo "$SEP"
echo "AGING"
echo "$SEP"
echo "Scan output dir contents:"
ls "$XQTL2/process/AGE_SY/AGE_SY20_M_freqs250/" 2>/dev/null || echo "  (missing or empty)"
echo ""

# ── PUPATION ──────────────────────────────────────────────────────────────────
echo "$SEP"
echo "PUPATION"
echo "$SEP"
echo "Non-RefAlt/hap files in pupal_dynamic/:"
ls "$SARAH/process/pupal_dynamic/" 2>/dev/null | grep -vE "RefAlt|R\.haps|\.out\.rds" || echo "  (none yet)"
echo ""

# ── ZINC2 ─────────────────────────────────────────────────────────────────────
echo "$SEP"
echo "ZINC2 — pseudoscan column names"
echo "$SEP"
safe_head "$XQTL2/process/ZINC2/ZINC2_M/ZINC2_M.pseudoscan.chr2L.txt" 2 "ZINC2_M pseudoscan chr2L"
echo ""
safe_head "$XQTL2/process/ZINC2/ZINC2_M/ZINC2_M.pseudoscan.txt" 2 "ZINC2_M pseudoscan concatenated"
echo ""
safe_head "$XQTL2/process/ZINC2/ZINC2_F/ZINC2_F.pseudoscan.txt" 2 "ZINC2_F pseudoscan concatenated"
echo ""

# ── MALATHION ─────────────────────────────────────────────────────────────────
echo "$SEP"
echo "MALATHION — parameter file and scan script"
echo "$SEP"
echo "R.pseudoscan.txt — line count and first 5 lines:"
safe_head "$MAL_BASE/R.pseudoscan.txt" 5
echo ""
echo "helperfiles/ listing:"
ls "$MAL_BASE/helperfiles/" 2>/dev/null || echo "  (empty or missing)"
echo ""
echo "First file in helperfiles/ (first 10 lines):"
FIRST=$(ls "$MAL_BASE/helperfiles/" 2>/dev/null | head -1)
[ -n "$FIRST" ] && safe_head "$MAL_BASE/helperfiles/$FIRST" 10 || echo "  (no files)"
echo ""
echo "haps2pseudoN_scan.R — first 80 lines:"
echo "---"
head -80 "$MAL_BASE/scripts/haps2pseudoN_scan.R" 2>/dev/null || echo "  (not found)"
echo "---"
echo ""
echo "Any existing scan .txt output in process/ (not RefAlt/haps/badloci):"
ls "$MAL_BASE/process/"*.txt 2>/dev/null | grep -vE "RefAlt|R\.haps|R\.bad" | head -10 || echo "  (none)"
echo ""

# ── SLURM ─────────────────────────────────────────────────────────────────────
echo "$SEP"
echo "SLURM — queue and recent job history"
echo "$SEP"
squeue -u "$USER"
echo ""
echo "sacct: jobs submitted today (last 20):"
sacct --starttime=$(date +%Y-%m-%d) -u "$USER" \
    --format=JobID%15,JobName%12,State%12,ExitCode,Elapsed,Submit \
    2>/dev/null | tail -25
echo ""
echo "###############################################################################"
echo "# END OF REPORT"
echo "###############################################################################"
