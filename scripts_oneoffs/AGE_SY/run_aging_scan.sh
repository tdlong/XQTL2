#!/bin/bash
###############################################################################
# run_aging_scan.sh
#
# Aging scan: males day 20 (AGE_SY20_M)
# Haplotype estimates already exist — goes straight to haps2scan → concat.
#
# HapParams: helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R (not needed for scan)
# TestFile:  helpfiles/AGE_SY/AGE_SY20_M.test.txt
# Haps:      process/AGE_SY/R.haps.chr*.rds  (already present)
# Scan out:  process/AGE_SY/AGE_SY20_M_freqs250/chr*.txt
# Final:     process/AGE_SY/AGE_SY20_M_freqs250/AGE_SY20_M_freqs250.pseudoscan.txt
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

set -e

TESTFILE="helpfiles/AGE_SY/AGE_SY20_M.test.txt"
PROCDIR="process/AGE_SY"
SCANNAME="AGE_SY20_M_freqs250"
FREQ_HALF=50

echo "=== Aging scan: AGE_SY20_M (250 kb freq smoothing) ==="
echo "Run from: /dfs7/adl/tdlong/fly_pool/XQTL2"
echo ""

if [ ! -f "$TESTFILE" ]; then
    echo "ERROR: TESTFILE not found: $TESTFILE"; exit 1
fi

mkdir -p "$PROCDIR"

# Step 1: Haplotypes → Wald scan (haps already exist)
echo "Submitting scan array job (chromosomes 1–5)..."
JOB_SCAN=$(sbatch --parsable --array=1-5 \
    -A tdlong_lab -p standard \
    scripts/haps2scan.freqsmooth.sh \
    "$TESTFILE" "$PROCDIR" "$SCANNAME" "$FREQ_HALF")
echo "  Scan job: ${JOB_SCAN}"
echo ""

# Step 2: Concatenate chromosome scans
echo "Submitting concat job (depends on Job ${JOB_SCAN})..."
JOB_CONCAT=$(sbatch --parsable \
    --dependency=afterok:${JOB_SCAN} \
    -A tdlong_lab -p standard \
    --wrap="bash scripts/concat_Chromosome_Scans.Andreas.sh ${PROCDIR}/${SCANNAME} && echo 'Concat done'")
echo "  Concat job: ${JOB_CONCAT}"
echo ""

echo "================================================================"
echo "All jobs submitted.  Monitor: squeue -u \$USER"
echo ""
echo "Job chain:"
echo "  Scan  (array 1-5): ${JOB_SCAN}"
echo "  Concat:            ${JOB_CONCAT}"
echo ""
echo "When complete, download the PNG:"
echo "  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/${PROCDIR}/${SCANNAME}/${SCANNAME}.5panel.png analysis/figures/"
echo "================================================================"
