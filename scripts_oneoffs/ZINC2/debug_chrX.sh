#!/bin/bash
###############################################################################
# debug_chrX.sh — Rerun smooth_haps for chrX, then hap_scan, then diagnostics
#
# Submits 3 chained SLURM jobs:
#   1. smooth_haps.R for chrX only
#   2. hap_scan.R for chrX only
#   3. chrX_diag.R — text diagnostics + frequency plot + Wald plot
#
# Results auto-push to XQTL2-dev so Claude can read them.
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
# Usage:    bash scripts_oneoffs/ZINC2/debug_chrX.sh
###############################################################################

set -e

# Pull pipeline code from XQTL2 (origin) — scripts/ lives there, not in XQTL2-dev
git pull origin main

DESIGN=helpfiles/ZINC2/Zinc2.test.F.txt
DIR=process/ZINC2
SCAN=ZINC2_F_v3
RESULTS=scripts_oneoffs/ZINC2/debug_chrX_results

mkdir -p $RESULTS

# 1. smooth_haps for chrX
jid_smooth=$(sbatch --parsable \
    -A tdlong_lab -p standard --cpus-per-task=2 --mem-per-cpu=6G --time=1:00:00 \
    --job-name=smooth_chrX \
    --wrap="module load R/4.2.2 && Rscript scripts/smooth_haps.R \
        --chr chrX --dir $DIR --outdir $SCAN \
        --rfile $DESIGN --smooth-kb 250")
echo "smooth: $jid_smooth"

# 2. hap_scan for chrX
jid_hap=$(sbatch --parsable --dependency=afterok:$jid_smooth \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=1:00:00 \
    --job-name=hapscan_chrX \
    --wrap="module load R/4.2.2 && Rscript scripts/hap_scan.R \
        --chr chrX --dir ${DIR}/${SCAN} --outdir ${SCAN} \
        --rfile $DESIGN")
echo "hap_scan: $jid_hap"

# 3. diagnostics (text + plots)
jid_diag=$(sbatch --parsable --dependency=afterok:$jid_hap \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=0:30:00 \
    --job-name=diag_chrX \
    --wrap="module load R/4.2.2 && Rscript scripts_oneoffs/ZINC2/debug_chrX_diag.R && \
        cd /dfs7/adl/tdlong/fly_pool/XQTL2 && \
        git add $RESULTS/ && \
        git commit -m 'debug_chrX results' && \
        git push dev HEAD:main")
echo "diag: $jid_diag"

echo ""
echo "Jobs chained: smooth($jid_smooth) -> hap_scan($jid_hap) -> diag($jid_diag)"
echo "Results will be in: $RESULTS/"
