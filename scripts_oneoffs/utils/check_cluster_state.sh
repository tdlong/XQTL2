#!/bin/bash
# check_cluster_state.sh
# Run on HPC3 from /dfs7/adl/tdlong/fly_pool/XQTL2 to verify all inputs
# needed for the new XQTL2 pipeline reanalysis are present.
# Usage: bash check_cluster_state.sh 2>&1 | tee check_cluster_state.out.txt
# Then scp the .out.txt back to the local project.

BASE=/dfs7/adl/tdlong/fly_pool/XQTL2
CHRS="chrX chr2L chr2R chr3L chr3R"

SEP="========================================================================"

echo "Run date: $(date)"
echo "Host: $(hostname)"
echo "Base dir: $BASE"
echo ""

# ── 1. Pipeline scripts ────────────────────────────────────────────────────
echo $SEP
echo "1. PIPELINE SCRIPTS"
echo $SEP
for f in scripts/run_scan.sh scripts/run_snp_scan.sh \
          scripts/smooth_haps.sh scripts/hap_scan.sh \
          scripts/snp_scan.sh scripts/concat_scans.sh \
          scripts/REFALT2haps.sh; do
    if [ -f "$BASE/$f" ]; then
        echo "  OK  $f"
    else
        echo "  MISSING  $f"
    fi
done
echo ""

# ── 2. SNP frequency tables ────────────────────────────────────────────────
echo $SEP
echo "2. SNP FREQUENCY TABLES"
echo $SEP
for f in helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
          helpfiles/FREQ_SNPs_Bpop.cM.txt.gz \
          helpfiles/FREQ_SNPs.cM.txt.gz; do
    if [ -f "$BASE/$f" ]; then
        echo "  OK   $f  ($(du -sh $BASE/$f | cut -f1))"
    else
        echo "  MISSING  $f"
    fi
done
echo ""

# ── 3. Haplotype files — ZINC_Hanson ──────────────────────────────────────
echo $SEP
echo "3. HAPLOTYPE FILES — ZINC_Hanson (pA)"
echo $SEP
echo "  process/ZINC_Hanson/ contents:"
ls $BASE/process/ZINC_Hanson/ 2>/dev/null || echo "  (directory missing)"
echo ""
echo "  Looking for R.haps.*.out.rds files:"
for chr in $CHRS; do
    f=$BASE/process/ZINC_Hanson/R.haps.${chr}.out.rds
    if [ -f "$f" ]; then
        echo "  OK   R.haps.${chr}.out.rds  ($(du -sh $f | cut -f1))"
    else
        echo "  MISSING  R.haps.${chr}.out.rds"
    fi
done
echo ""

# ── 4. Haplotype files — ZINC2 (pB) ───────────────────────────────────────
echo $SEP
echo "4. HAPLOTYPE FILES — ZINC2 (pB)"
echo $SEP
echo "  process/ZINC2/ contents:"
ls $BASE/process/ZINC2/ 2>/dev/null || echo "  (directory missing)"
echo ""
echo "  Looking for R.haps.*.out.rds files:"
for chr in $CHRS; do
    f=$BASE/process/ZINC2/R.haps.${chr}.out.rds
    if [ -f "$f" ]; then
        echo "  OK   R.haps.${chr}.out.rds  ($(du -sh $f | cut -f1))"
    else
        echo "  MISSING  R.haps.${chr}.out.rds"
    fi
done
echo ""

# ── 5. Helpfiles — ZINC_Hanson ─────────────────────────────────────────────
echo $SEP
echo "5. HELPFILES — ZINC_Hanson (pA)"
echo $SEP
echo "  helpfiles/ZINC_Hanson/ contents:"
ls $BASE/helpfiles/ZINC_Hanson/ 2>/dev/null || echo "  (directory missing)"
echo ""
for f in helpfiles/ZINC_Hanson/hap_params.R \
          helpfiles/ZINC_Hanson/design.txt; do
    echo "  --- $f ---"
    if [ -f "$BASE/$f" ]; then
        cat $BASE/$f
    else
        echo "  (not found)"
    fi
    echo ""
done

# ── 6. Helpfiles — ZINC2 (pB) ─────────────────────────────────────────────
echo $SEP
echo "6. HELPFILES — ZINC2 (pB)"
echo $SEP
echo "  helpfiles/ZINC2/ contents:"
ls $BASE/helpfiles/ZINC2/ 2>/dev/null || echo "  (directory missing)"
echo ""
for f in helpfiles/ZINC2/hap_params.R \
          helpfiles/ZINC2/design.txt \
          helpfiles/ZINC2/design_F.txt \
          helpfiles/ZINC2/design_M.txt; do
    echo "  --- $f ---"
    if [ -f "$BASE/$f" ]; then
        cat $BASE/$f
    else
        echo "  (not found)"
    fi
    echo ""
done

# ── 7. Existing scan output ────────────────────────────────────────────────
echo $SEP
echo "7. EXISTING SCAN OUTPUT (to avoid overwriting)"
echo $SEP
echo "  process/ZINC_Hanson/:"
ls $BASE/process/ZINC_Hanson/ 2>/dev/null | grep -v "^R\." || echo "  (none)"
echo ""
echo "  process/ZINC2/:"
ls $BASE/process/ZINC2/ 2>/dev/null | grep -v "^R\." || echo "  (none)"
echo ""

# ── 8. Founder BAM list ────────────────────────────────────────────────────
echo $SEP
echo "8. FOUNDER BAM LIST"
echo $SEP
if [ -f "$BASE/helpfiles/founder.bams.txt" ]; then
    echo "  OK  helpfiles/founder.bams.txt:"
    cat $BASE/helpfiles/founder.bams.txt
else
    echo "  MISSING  helpfiles/founder.bams.txt"
    echo "  Checking data/founders/ directly:"
    ls $BASE/data/founders/*.bam 2>/dev/null | head -20 || echo "  (not found)"
fi
echo ""

echo $SEP
echo "DONE. Paste the contents of this file back to Claude."
echo $SEP
