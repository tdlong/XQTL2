#!/bin/bash
###############################################################################
# run_stable_all.sh — Run stabilized Wald scan on ALL datasets, one download
#
# Uses EXISTING haplotype estimates (R.haps.*.rds) — no Step 3 re-run.
# Submits all scan + concat jobs with SLURM dependencies.
# Final job tars all stable results into one downloadable file.
#
# Scans (11 total):
#   ZINC2_F_stable         pB females, 15 reps
#   ZINC2_M_stable         pB males, 15 reps
#   ZINC_Hanson_stable     pA females, 12 reps
#   ZINC2_F_stable_N07..N13  downsampled females (7,9,11,13 reps)
#   ZINC2_M_stable_N07..N13  downsampled males (7,9,11,13 reps)
#
# Each scan = 5-chromosome array job → concat → tar
# Total: 55 array tasks + 11 concat + 1 tar = 67 jobs
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
#
# DEPLOYMENT:
#   scp pipeline/run_stable_all.sh tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/
#   ssh hpc3 "cd /dfs7/adl/tdlong/fly_pool/XQTL2 && bash run_stable_all.sh"
#
# DOWNLOAD (when all jobs done):
#   scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/stable_results.tar.gz analysis/data/
###############################################################################

set -e

echo "=== Stabilized Wald scan — all datasets ==="
echo ""

# ── Verify haplotype files exist ─────────────────────────────────────────────
for DIR in process/ZINC2_v2 process/ZINC_Hanson_v2; do
    for CHR in chrX chr2L chr2R chr3L chr3R; do
        if [ ! -f "${DIR}/R.haps.${CHR}.out.rds" ]; then
            echo "ERROR: Missing ${DIR}/R.haps.${CHR}.out.rds"
            echo "Run run_ZINC2.sh / run_ZINC_Hanson.sh Step 3 first."
            exit 1
        fi
    done
done
echo "All haplotype files found."
echo ""

# ── Generate downsampled test files ──────────────────────────────────────────
HELPDIR="helpfiles/ZINC2"
DOWNSAMPLE_DIR="${HELPDIR}/downsample"
mkdir -p "$DOWNSAMPLE_DIR"

echo "Generating downsampled test files..."
for N in 7 9 11 13; do
    NPAD=$(printf "%02d" $N)
    NLINES=$((2 * N + 1))   # header + 2*N data lines (control + zinc per rep)
    for SEX in F M; do
        SRC="${HELPDIR}/Zinc2.test.${SEX}.txt"
        DST="${DOWNSAMPLE_DIR}/Zinc2.test.${SEX}.N${NPAD}.txt"
        head -n "$NLINES" "$SRC" > "$DST"
    done
done
echo "  Done."
echo ""

# ── Submit scan + concat jobs ────────────────────────────────────────────────
ALL_CONCAT_JOBS=""

submit_scan() {
    local TESTFILE=$1
    local PROCDIR=$2
    local SCANNAME=$3

    JOB=$(sbatch --parsable -A tdlong_lab -p standard \
        scripts/haps2scan.stable.sh \
        "$TESTFILE" "$PROCDIR" "$SCANNAME")
    echo "  Scan  ${SCANNAME}: Job $JOB"

    JOB_CONCAT=$(sbatch --parsable -A tdlong_lab -p standard \
        --dependency=afterok:${JOB} \
        --wrap="bash scripts/concat_Chromosome_Scans.Andreas.sh ${PROCDIR}/${SCANNAME}")
    echo "  Concat ${SCANNAME}: Job $JOB_CONCAT"

    ALL_CONCAT_JOBS="${ALL_CONCAT_JOBS}:${JOB_CONCAT}"
}

echo "Full-data scans..."
submit_scan "helpfiles/ZINC2/Zinc2.test.F.txt"             "process/ZINC2_v2"       "ZINC2_F_stable"
submit_scan "helpfiles/ZINC2/Zinc2.test.M.txt"             "process/ZINC2_v2"       "ZINC2_M_stable"
submit_scan "helpfiles/ZINC_Hanson/ZINC_Hanson.test.txt"   "process/ZINC_Hanson_v2" "ZINC_Hanson_stable"
echo ""

echo "Downsampled scans..."
for N in 7 9 11 13; do
    NPAD=$(printf "%02d" $N)
    for SEX in F M; do
        submit_scan \
            "${DOWNSAMPLE_DIR}/Zinc2.test.${SEX}.N${NPAD}.txt" \
            "process/ZINC2_v2" \
            "ZINC2_${SEX}_stable_N${NPAD}"
    done
done
echo ""

# ── Final tar job: collect all stable results ────────────────────────────────
# Remove leading colon
ALL_CONCAT_JOBS="${ALL_CONCAT_JOBS#:}"

TAR_JOB=$(sbatch --parsable -A tdlong_lab -p standard \
    --dependency=afterok:${ALL_CONCAT_JOBS} \
    --wrap="cd /dfs7/adl/tdlong/fly_pool/XQTL2 && \
tar czf process/stable_results.tar.gz \
  process/ZINC2_v2/ZINC2_F_stable \
  process/ZINC2_v2/ZINC2_M_stable \
  process/ZINC_Hanson_v2/ZINC_Hanson_stable \
  process/ZINC2_v2/ZINC2_F_stable_N07 \
  process/ZINC2_v2/ZINC2_F_stable_N09 \
  process/ZINC2_v2/ZINC2_F_stable_N11 \
  process/ZINC2_v2/ZINC2_F_stable_N13 \
  process/ZINC2_v2/ZINC2_M_stable_N07 \
  process/ZINC2_v2/ZINC2_M_stable_N09 \
  process/ZINC2_v2/ZINC2_M_stable_N11 \
  process/ZINC2_v2/ZINC2_M_stable_N13")

echo "Final tar: Job $TAR_JOB"
echo ""
echo "================================================================"
echo "All jobs submitted. Monitor with: squeue -u \$USER"
echo ""
echo "When done, download ONE file:"
echo "  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/stable_results.tar.gz analysis/data/"
echo ""
echo "Then locally:"
echo "  cd analysis/data && tar xzf stable_results.tar.gz"
echo "================================================================"
