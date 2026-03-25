#!/bin/bash
###############################################################################
# run_chr2R_freqsmooth.sh
#
# Diagnostic run: chr2R ONLY for pB females and males, two frequency-smoothing
# bandwidths (100 kb and 250 kb), all 4 jobs submitted in parallel.
#
# Purpose: Determine whether smoothing haplotype frequency vectors before the
# Wald test reduces the noise visible in chr2R Wald profiles.
#
# Assumes Step 3 (haplotype files) already exists at process/ZINC2_v2/
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
#   scp pipeline/run_chr2R_freqsmooth.sh tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/
#   scp pipeline/XQTL2_scripts/haps2scan.freqsmooth.R tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/scripts/
#   scp pipeline/XQTL2_scripts/haps2scan.freqsmooth.sh tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/scripts/
#   scp pipeline/XQTL2_scripts/haps2scan.stable.code.R tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/scripts/
#   ssh hpc3 "cd /dfs7/adl/tdlong/fly_pool/XQTL2 && bash run_chr2R_freqsmooth.sh"
#
# When done, scp down the 4 chr2R pseudoscan files:
#   scp 'tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2_v2/ZINC2_*_freqs*/ZINC2_*_freqs*.pseudoscan.chr2R.txt' analysis/data/process/ZINC2_v2/
###############################################################################

set -e

mkdir -p process/ZINC2_v2

echo "=== chr2R frequency-smoothing diagnostic ==="
echo "Submitting 4 jobs (F+M × 100kb+250kb), all chr2R only"
echo ""

# ── 100 kb smoothing (FREQ_SMOOTH_HALF = 20 windows × 5 kb = 100 kb) ─────────
JOB1=$(sbatch --parsable \
    scripts/haps2scan.freqsmooth.sh \
    helpfiles/ZINC2/Zinc2.test.F.txt \
    "process/ZINC2_v2" \
    "ZINC2_F_freqs100" \
    20)
echo "ZINC2_F freqs100 (chr2R): Job $JOB1"

JOB2=$(sbatch --parsable \
    scripts/haps2scan.freqsmooth.sh \
    helpfiles/ZINC2/Zinc2.test.M.txt \
    "process/ZINC2_v2" \
    "ZINC2_M_freqs100" \
    20)
echo "ZINC2_M freqs100 (chr2R): Job $JOB2"

# ── 250 kb smoothing (FREQ_SMOOTH_HALF = 50 windows × 5 kb = 250 kb) ─────────
JOB3=$(sbatch --parsable \
    scripts/haps2scan.freqsmooth.sh \
    helpfiles/ZINC2/Zinc2.test.F.txt \
    "process/ZINC2_v2" \
    "ZINC2_F_freqs250" \
    50)
echo "ZINC2_F freqs250 (chr2R): Job $JOB3"

JOB4=$(sbatch --parsable \
    scripts/haps2scan.freqsmooth.sh \
    helpfiles/ZINC2/Zinc2.test.M.txt \
    "process/ZINC2_v2" \
    "ZINC2_M_freqs250" \
    50)
echo "ZINC2_M freqs250 (chr2R): Job $JOB4"

echo ""
echo "Monitor: squeue -u \$USER"
echo ""
echo "When all 4 jobs finish, download:"
echo "  scp 'tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2_v2/ZINC2_*_freqs*/ZINC2_*_freqs*.pseudoscan.chr2R.txt' analysis/data/process/ZINC2_v2/"
echo ""
echo "Then run locally:  Rscript analysis/scripts/diag_freqsmooth.R"
