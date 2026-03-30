#!/bin/bash
###############################################################################
# debug_chrX_v2.sh — Rerun smooth for chrX with instrumented smooth_haps,
#                     then hap_scan, then frequency check plot.
#
# One command:
#   cd /dfs7/adl/tdlong/fly_pool/XQTL2 && bash scripts_oneoffs/ZINC2/debug_chrX_v2.sh
#
# All diagnostic output goes to scripts_oneoffs/ZINC2/debug_chrX_results/
# (git-trackable, NOT in process/). Results auto-pushed to XQTL2-dev.
###############################################################################

set -e
cd /dfs7/adl/tdlong/fly_pool/XQTL2

# Pull latest code from both repos
git pull dev main --rebase 2>/dev/null || true
git pull origin main 2>/dev/null || true

DESIGN=helpfiles/ZINC2/Zinc2.test.F.txt
DIR=process/ZINC2
SCAN=ZINC2_F_v3
RESDIR=scripts_oneoffs/ZINC2/debug_chrX_results
mkdir -p ${RESDIR}

# ── Job 1: Instrumented smooth (chrX only) ────────────────────────────────────
#     Writes diag to RESDIR (git-trackable), pipeline output to process/ as usual
JID_SMOOTH=$(sbatch --parsable \
    -A tdlong_lab -p standard \
    --cpus-per-task=2 --mem-per-cpu=6G --time=4:00:00 \
    --job-name=smooth_dbg_chrX \
    --wrap="module load R/4.2.2 && \
Rscript scripts_oneoffs/ZINC2/smooth_haps_debug.R \
    --chr chrX --dir ${DIR} --outdir ${SCAN} \
    --rfile ${DESIGN} --smooth-kb 250 \
    --diagdir ${RESDIR}")
echo "smooth: ${JID_SMOOTH}"

# ── Job 2: Haplotype scan (depends on smooth) ────────────────────────────────
JID_SCAN=$(sbatch --parsable \
    --dependency=afterok:${JID_SMOOTH} \
    -A tdlong_lab -p standard \
    --cpus-per-task=1 --mem-per-cpu=3G --time=2:00:00 \
    --job-name=hscan_dbg_chrX \
    --wrap="module load R/4.2.2 && \
Rscript scripts/hap_scan.R \
    --chr chrX --dir ${DIR}/${SCAN} --outdir ${SCAN} \
    --rfile ${DESIGN}")
echo "hap_scan: ${JID_SCAN}"

# ── Job 3: Plots + push all results ───────────────────────────────────────────
JID_DIAG=$(sbatch --parsable \
    --dependency=afterany:${JID_SMOOTH},afterany:${JID_SCAN} \
    -A tdlong_lab -p standard \
    --cpus-per-task=1 --mem-per-cpu=3G --time=0:30:00 \
    --job-name=diag_dbg_chrX \
    --wrap="cd /dfs7/adl/tdlong/fly_pool/XQTL2 && \
module load R/4.2.2 && \
Rscript scripts_oneoffs/ZINC2/debug_chrX_plots.R && \
git add ${RESDIR}/ && \
git commit -m 'debug chrX v2: diag + plots' && \
git pull dev main --rebase && \
git push dev HEAD:main || echo 'WARNING: git push failed'")
echo "diag: ${JID_DIAG}"

echo ""
echo "Jobs: smooth=${JID_SMOOTH} -> hap_scan=${JID_SCAN} -> diag=${JID_DIAG}"
echo "Monitor: squeue -u \$USER"
echo "Results will be in: ${RESDIR}/"
