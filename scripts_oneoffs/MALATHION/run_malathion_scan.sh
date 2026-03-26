#!/bin/bash
###############################################################################
# run_malathion_scan.sh
#
# Malathion resistance scan — pA (A-founder panel) experiment
# Population: pA panel (founders A1–A7, AB8)
# Design: 2 sexes treated as 2 replicates (no replicates within sex)
#   Pools: c_F_1 (control female), s_F_1 (selected female)
#          c_M_1 (control male),   s_M_1 (selected male)
#
# Pipeline: REFALT2haps → haps2scan → concat
#   HapParams: helpfiles/MALATHION/malathion_haplotype_parameters.R
#   TestFile:  helpfiles/MALATHION/malathion.test.txt
#   RefAlt:    *** SEE REFALT_DIR BELOW — confirm path with find_refalt_files.sh ***
#   Hap output:  process/MALATHION_v1/R.haps.chr*.rds
#   Scan output: process/MALATHION_v1/MALATHION_freqs250/chr*.txt
#   Final: process/MALATHION_v1/MALATHION_freqs250/MALATHION_freqs250.pseudoscan.txt
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

set -e

# ── VERIFY THIS PATH with find_refalt_files.sh output ────────────────────────
REFALT_DIR="/dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/process"
# Expected files: ${REFALT_DIR}/RefAlt.chrX.txt, RefAlt.chr2L.txt, etc.

HAP_PARAMS="helpfiles/MALATHION/malathion_haplotype_parameters.R"
TESTFILE="helpfiles/MALATHION/malathion.test.txt"
PROCDIR="process/MALATHION_v1"
SCANNAME="MALATHION_freqs250"
FREQ_HALF=50

echo "=== Malathion resistance scan, pA founders (250 kb freq smoothing) ==="
echo "Run from: /dfs7/adl/tdlong/fly_pool/XQTL2"
echo ""

# Sanity checks
if [ ! -f "$HAP_PARAMS" ]; then
    echo "ERROR: HAP_PARAMS not found: $HAP_PARAMS"; exit 1
fi
if [ ! -f "$TESTFILE" ]; then
    echo "ERROR: TESTFILE not found: $TESTFILE"; exit 1
fi

# Check at least one RefAlt file exists
if [ ! -f "${REFALT_DIR}/RefAlt.chr2L.txt" ]; then
    echo "ERROR: RefAlt files not found in: ${REFALT_DIR}"
    echo "       Run scripts/find_refalt_files.sh to locate them, then update REFALT_DIR above."
    exit 1
fi

mkdir -p "${PROCDIR}"

# Symlink RefAlt files into PROCDIR
for chr in chrX chr2L chr2R chr3L chr3R; do
    SRC="${REFALT_DIR}/RefAlt.${chr}.txt"
    DST="${PROCDIR}/RefAlt.${chr}.txt"
    if [ ! -e "${DST}" ] && [ -f "${SRC}" ]; then
        ln -s "${SRC}" "${DST}"
        echo "  Symlinked: ${DST}"
    fi
done
echo ""

# Step 1: RefAlt → Haplotypes
echo "Submitting REFALT2haps array job (chromosomes 1–5)..."
JOB_HAPS=$(sbatch --parsable --array=1-5 \
    -A tdlong_lab -p standard \
    scripts/REFALT2haps.Andreas.sh \
    "${HAP_PARAMS}" \
    "${PROCDIR}")
echo "  Haps job: ${JOB_HAPS}"
echo ""

# Step 2: Haplotypes → Wald scan
echo "Submitting scan array job (depends on Job ${JOB_HAPS})..."
JOB_SCAN=$(sbatch --parsable --array=1-5 \
    --dependency=afterok:${JOB_HAPS} \
    -A tdlong_lab -p standard \
    scripts/haps2scan.freqsmooth.sh \
    "${TESTFILE}" "${PROCDIR}" "${SCANNAME}" "${FREQ_HALF}")
echo "  Scan job: ${JOB_SCAN}"
echo ""

# Step 3: Concatenate chromosome scans
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
echo "  Haps  (array 1-5): ${JOB_HAPS}"
echo "  Scan  (array 1-5): ${JOB_SCAN}"
echo "  Concat:            ${JOB_CONCAT}"
echo ""
echo "When complete, download the PNG:"
echo "  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/${PROCDIR}/${SCANNAME}/${SCANNAME}.5panel.png analysis/figures/"
echo "================================================================"
