#!/bin/bash
###############################################################################
# diag_smooth_chrX.sh — Rerun smooth_haps.R for chrX on ZINC2_F only,
#                        capturing masking log to diagnose the coverage gap.
#
# Context: v3 post-fix scan has no windows in chrX 20.56–22.6 Mb.
# smooth_haps.R is supposed to interpolate through masked gaps, but that
# region has zero coverage. This job reruns the smooth step for chrX and
# captures verbose output so we can see how many founder-windows are masked
# and whether the gap extends beyond the smoothing kernel range.
#
# Overwrites: process/ZINC2/ZINC2_F_v3/ZINC2_F_v3.smooth.chrX.rds
#             process/ZINC2/ZINC2_F_v3/ZINC2_F_v3.meansBySample.chrX.txt
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

set -e
mkdir -p logs

# ── Job 1: lightweight gap diagnostic (reads RDS, no smoothing, ~3G) ─────────
jid_gap=$(sbatch --parsable \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=0:30:00 \
    --job-name=diag_chrX_gap \
    --output=logs/diag_chrX_gap.out \
    --wrap="module load R/4.2.2 && Rscript scripts_oneoffs/ZINC2/diag_chrX_gap.R")
echo "diag_chrX_gap:    ${jid_gap}"

# ── Job 2: full smooth rerun (2x6G, captures masking + covariance log) ───────
jid_smooth=$(sbatch --parsable \
    -A tdlong_lab -p standard --cpus-per-task=2 --mem-per-cpu=6G --time=1:00:00 \
    --job-name=diag_smooth_chrX \
    --output=logs/diag_smooth_chrX.out \
    --wrap="module load R/4.2.2 && \
Rscript scripts/smooth_haps.R \
    --chr       chrX \
    --dir       process/ZINC2 \
    --outdir    ZINC2_F_v3 \
    --rfile     helpfiles/ZINC2/Zinc2.test.F.txt \
    --smooth-kb 250")
echo "diag_smooth_chrX: ${jid_smooth}"

echo ""
echo "When complete:"
echo "  cat logs/diag_chrX_gap.out"
echo "  cat logs/diag_smooth_chrX.out"
