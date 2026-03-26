#!/bin/bash
###############################################################################
# run_pupation_scan.sh
#
# Pupation height scan: Top vs Middle comparison
# Population: DSPR pB panel (founders B1–B7, AB8)
#
# Haplotype estimates already exist in Sarah's directory.
# This script symlinks them into a new PROCDIR in tdlong's space,
# then runs haps2scan → concat.
#
# Data (Sarah Ruckman):
#   Haps:     /dfs7/adl/sruckman/XQTL/XQTL2/process/pupal_dynamic/R.haps.chr*.rds
#   TestFile: /dfs7/adl/sruckman/XQTL/XQTL2/helpfiles/TMparameter.txt
#
# Output (tdlong's XQTL2):
#   Hap symlinks: process/PUPATION_v1/R.haps.chr*.rds
#   Scan output:  process/PUPATION_v1/PUPATION_TM_freqs250/chr*.txt
#   Final:        process/PUPATION_v1/PUPATION_TM_freqs250/PUPATION_TM_freqs250.pseudoscan.txt
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

set -e

SARAH="/dfs7/adl/sruckman/XQTL/XQTL2"
TESTFILE="${SARAH}/helpfiles/TMparameter.txt"
PROCDIR="process/PUPATION_v1"
SCANNAME="PUPATION_TM_freqs250"
FREQ_HALF=50

echo "=== Pupation Top vs Middle scan (250 kb freq smoothing) ==="
echo "Run from: /dfs7/adl/tdlong/fly_pool/XQTL2"
echo ""

if [ ! -f "$TESTFILE" ]; then
    echo "ERROR: TESTFILE not found: $TESTFILE"; exit 1
fi

mkdir -p "${PROCDIR}"

# Symlink hap files from Sarah's directory (both .rds and .out.rds variants)
for chr in chrX chr2L chr2R chr3L chr3R; do
    for ext in rds out.rds; do
        SRC="${SARAH}/process/pupal_dynamic/R.haps.${chr}.${ext}"
        DST="${PROCDIR}/R.haps.${chr}.${ext}"
        if [ ! -e "${DST}" ] && [ -f "${SRC}" ]; then
            ln -s "${SRC}" "${DST}"
            echo "  Symlinked: ${DST}"
        fi
    done
done
echo ""

# Step 1: Haplotypes → Wald scan
echo "Submitting scan array job (chromosomes 1–5)..."
JOB_SCAN=$(sbatch --parsable --array=1-5 \
    -A tdlong_lab -p standard \
    scripts/haps2scan.freqsmooth.sh \
    "${TESTFILE}" "${PROCDIR}" "${SCANNAME}" "${FREQ_HALF}")
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
echo "  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/${PROCDIR}/${SCANNAME}/${SCANNAME}.5panel.cM.png analysis/figures/"
echo "================================================================"
