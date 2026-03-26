#!/bin/bash
###############################################################################
# run_freqs250_all.sh
#
# Full genome run with frequency smoothing (250 kb) + covariance smoothing
# + gap filling.  Definitive pipeline for manuscript.
#
# Scans (11 total, all using haps2scan.freqsmooth with FREQ_SMOOTH_HALF=50):
#   ZINC2_F_freqs250       pB females, 15 reps
#   ZINC2_M_freqs250       pB males, 15 reps
#   ZINC_Hanson_freqs250   pA females, 12 reps
#   ZINC2_F_freqs250_N07..N13   downsampled females
#   ZINC2_M_freqs250_N07..N13   downsampled males
#
# Assumes haplotype files (R.haps.*.rds) already exist from previous Step 3.
# NOTE: ZINC2_F_freqs250 and ZINC2_M_freqs250 may already have chr2R from the
#       diagnostic run — full run adds remaining chromosomes and overwrites chr2R
#       with the gap-filled version.
#
# Outputs:
#   Tarball 1 (download):  process/main_results.tar.gz
#     — ZINC2_F_freqs250, ZINC2_M_freqs250, ZINC_Hanson_freqs250
#   Tarball 2 (stay on cluster): process/downsample_results.tar.gz
#     — all 8 downsampled scans
#   Extract (download):    process/ZINC2_v2/downsample_pseudoscans.tsv.gz
#     — all downsampled pseudoscans in one compact file
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

set -e

echo "=== Full genome run — freq smoothing 250 kb + gap filling ==="
echo ""

mkdir -p process/ZINC2_v2 process/ZINC_Hanson_v2

# ── Generate downsampled test files ──────────────────────────────────────────
HELPDIR="helpfiles/ZINC2"
DOWNSAMPLE_DIR="${HELPDIR}/downsample"
mkdir -p "$DOWNSAMPLE_DIR"

echo "Generating downsampled test files..."
for N in 7 9 11 13; do
    NPAD=$(printf "%02d" $N)
    NLINES=$((2 * N + 1))
    for SEX in F M; do
        SRC="${HELPDIR}/Zinc2.test.${SEX}.txt"
        DST="${DOWNSAMPLE_DIR}/Zinc2.test.${SEX}.N${NPAD}.txt"
        head -n "$NLINES" "$SRC" > "$DST"
    done
done
echo "  Done."
echo ""

# ── Helper: submit one scan (all 5 chr array) + concat ───────────────────────
ALL_MAIN_CONCAT=""
ALL_DS_CONCAT=""

submit_scan() {
    local TESTFILE=$1
    local PROCDIR=$2
    local SCANNAME=$3
    local FREQ_HALF=$4
    local DEST_LIST=$5   # "main" or "ds"

    JOB=$(sbatch --parsable --array=1-5 \
        scripts/haps2scan.freqsmooth.sh \
        "$TESTFILE" "$PROCDIR" "$SCANNAME" "$FREQ_HALF")
    echo "  Scan  ${SCANNAME}: Job $JOB"

    JOB_CONCAT=$(sbatch --parsable \
        --dependency=afterok:${JOB} \
        -A tdlong_lab -p standard \
        --wrap="bash scripts/concat_Chromosome_Scans.Andreas.sh ${PROCDIR}/${SCANNAME}")
    echo "  Concat ${SCANNAME}: Job $JOB_CONCAT"

    if [ "$DEST_LIST" = "main" ]; then
        ALL_MAIN_CONCAT="${ALL_MAIN_CONCAT}:${JOB_CONCAT}"
    else
        ALL_DS_CONCAT="${ALL_DS_CONCAT}:${JOB_CONCAT}"
    fi
}

# ── Main scans (3) ───────────────────────────────────────────────────────────
echo "Main scans (pB F, pB M, Hanson)..."
submit_scan "helpfiles/ZINC2/Zinc2.test.F.txt"           "process/ZINC2_v2"       "ZINC2_F_freqs250"     50 "main"
submit_scan "helpfiles/ZINC2/Zinc2.test.M.txt"           "process/ZINC2_v2"       "ZINC2_M_freqs250"     50 "main"
submit_scan "helpfiles/ZINC_Hanson/ZINC_Hanson.test.txt" "process/ZINC_Hanson_v2" "ZINC_Hanson_freqs250" 50 "main"
echo ""

# ── Downsampled scans (8) ────────────────────────────────────────────────────
echo "Downsampled scans (F+M x N07,09,11,13)..."
for N in 7 9 11 13; do
    NPAD=$(printf "%02d" $N)
    for SEX in F M; do
        submit_scan \
            "${DOWNSAMPLE_DIR}/Zinc2.test.${SEX}.N${NPAD}.txt" \
            "process/ZINC2_v2" \
            "ZINC2_${SEX}_freqs250_N${NPAD}" \
            50 "ds"
    done
done
echo ""

# ── Tarball 1: main results (download this) ───────────────────────────────────
ALL_MAIN_CONCAT="${ALL_MAIN_CONCAT#:}"

TAR1=$(sbatch --parsable \
    --dependency=afterok:${ALL_MAIN_CONCAT} \
    -A tdlong_lab -p standard \
    --wrap="tar czf process/main_results.tar.gz \
        process/ZINC2_v2/ZINC2_F_freqs250 \
        process/ZINC2_v2/ZINC2_M_freqs250 \
        process/ZINC_Hanson_v2/ZINC_Hanson_freqs250 && \
        echo 'Tarball 1 done: process/main_results.tar.gz'")
echo "Tarball 1 (main results): Job $TAR1"

# ── Extract + Tarball 2: downsampling ────────────────────────────────────────
ALL_DS_CONCAT="${ALL_DS_CONCAT#:}"

EXTRACT=$(sbatch --parsable \
    --dependency=afterok:${ALL_DS_CONCAT} \
    -A tdlong_lab -p standard \
    --wrap="module load R/4.2.2 && Rscript scripts/extract_downsampling.R && \
        echo 'Extract done: process/ZINC2_v2/downsample_pseudoscans.tsv.gz'")
echo "Downsampling extract: Job $EXTRACT"

TAR2=$(sbatch --parsable \
    --dependency=afterok:${ALL_DS_CONCAT} \
    -A tdlong_lab -p standard \
    --wrap="tar czf process/downsample_results.tar.gz \
        process/ZINC2_v2/ZINC2_F_freqs250_N07 \
        process/ZINC2_v2/ZINC2_F_freqs250_N09 \
        process/ZINC2_v2/ZINC2_F_freqs250_N11 \
        process/ZINC2_v2/ZINC2_F_freqs250_N13 \
        process/ZINC2_v2/ZINC2_M_freqs250_N07 \
        process/ZINC2_v2/ZINC2_M_freqs250_N09 \
        process/ZINC2_v2/ZINC2_M_freqs250_N11 \
        process/ZINC2_v2/ZINC2_M_freqs250_N13 && \
        echo 'Tarball 2 done: process/downsample_results.tar.gz'")
echo "Tarball 2 (downsampling, stays on cluster): Job $TAR2"

echo ""
echo "================================================================"
echo "All jobs submitted. Monitor: squeue -u \$USER"
echo ""
echo "When Job $TAR1 finishes, download tarball 1:"
echo "  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/main_results.tar.gz analysis/data/"
echo "  cd analysis/data && tar xzf main_results.tar.gz && rm main_results.tar.gz"
echo ""
echo "When Job $EXTRACT finishes, download the downsampling extract:"
echo "  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2_v2/downsample_pseudoscans.tsv.gz analysis/data/process/ZINC2_v2/"
echo "================================================================"
