#!/bin/bash
###############################################################################
# submit_freqs100_scans.sh
#
# Submit 100 kb smoothing scans for all 4 experiments.
# Malathion depends on RefAlt2haps job 49840701 (still running).
# All others start immediately — haplotype files already exist.
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

FREQ_HALF=20   # 20 windows × 5 kb step = 100 kb smoothing
SARAH="/dfs7/adl/sruckman/XQTL/XQTL2"

# ── MALATHION — depends on RefAlt2haps 49840701 ──────────────────────────────
JOB_MAL_SCAN=$(sbatch --parsable --array=1-5 \
    --dependency=afterok:49840701 \
    -A tdlong_lab -p standard \
    scripts/haps2scan.freqsmooth.sh \
    helpfiles/MALATHION/malathion.test.txt \
    process/MALATHION_v1 \
    MALATHION_freqs100 \
    ${FREQ_HALF})
echo "Malathion scan  (freqs100): ${JOB_MAL_SCAN}"

JOB_MAL_CONCAT=$(sbatch --parsable \
    --dependency=afterok:${JOB_MAL_SCAN} \
    -A tdlong_lab -p standard \
    --wrap="bash scripts/concat_Chromosome_Scans.Andreas.sh process/MALATHION_v1/MALATHION_freqs100 && echo done")
echo "Malathion concat (freqs100): ${JOB_MAL_CONCAT}"

# ── AGING — haps at process/AGE_SY/ ──────────────────────────────────────────
JOB_AGE_SCAN=$(sbatch --parsable --array=1-5 \
    -A tdlong_lab -p standard \
    scripts/haps2scan.freqsmooth.sh \
    helpfiles/AGE_SY/AGE_SY20_M.test.txt \
    process/AGE_SY \
    AGE_SY20_M_freqs100 \
    ${FREQ_HALF})
echo "Aging scan      (freqs100): ${JOB_AGE_SCAN}"

JOB_AGE_CONCAT=$(sbatch --parsable \
    --dependency=afterok:${JOB_AGE_SCAN} \
    -A tdlong_lab -p standard \
    --wrap="bash scripts/concat_Chromosome_Scans.Andreas.sh process/AGE_SY/AGE_SY20_M_freqs100 && echo done")
echo "Aging concat    (freqs100): ${JOB_AGE_CONCAT}"

# ── ZINC2 MALE — haps at process/ZINC2/ ──────────────────────────────────────
JOB_ZINC_M_SCAN=$(sbatch --parsable --array=1-5 \
    -A tdlong_lab -p standard \
    scripts/haps2scan.freqsmooth.sh \
    helpfiles/ZINC2/Zinc2.test.M.txt \
    process/ZINC2 \
    ZINC2_M_freqs100 \
    ${FREQ_HALF})
echo "ZINC2 M scan    (freqs100): ${JOB_ZINC_M_SCAN}"

JOB_ZINC_M_CONCAT=$(sbatch --parsable \
    --dependency=afterok:${JOB_ZINC_M_SCAN} \
    -A tdlong_lab -p standard \
    --wrap="bash scripts/concat_Chromosome_Scans.Andreas.sh process/ZINC2/ZINC2_M_freqs100 && echo done")
echo "ZINC2 M concat  (freqs100): ${JOB_ZINC_M_CONCAT}"

# ── ZINC2 FEMALE — haps at process/ZINC2/ ────────────────────────────────────
JOB_ZINC_F_SCAN=$(sbatch --parsable --array=1-5 \
    -A tdlong_lab -p standard \
    scripts/haps2scan.freqsmooth.sh \
    helpfiles/ZINC2/Zinc2.test.F.txt \
    process/ZINC2 \
    ZINC2_F_freqs100 \
    ${FREQ_HALF})
echo "ZINC2 F scan    (freqs100): ${JOB_ZINC_F_SCAN}"

JOB_ZINC_F_CONCAT=$(sbatch --parsable \
    --dependency=afterok:${JOB_ZINC_F_SCAN} \
    -A tdlong_lab -p standard \
    --wrap="bash scripts/concat_Chromosome_Scans.Andreas.sh process/ZINC2/ZINC2_F_freqs100 && echo done")
echo "ZINC2 F concat  (freqs100): ${JOB_ZINC_F_CONCAT}"

# ── PUPATION — haps symlinked at process/PUPATION_v1/ ────────────────────────
JOB_PUP_SCAN=$(sbatch --parsable --array=1-5 \
    -A tdlong_lab -p standard \
    scripts/haps2scan.freqsmooth.sh \
    "${SARAH}/helpfiles/TMparameter.txt" \
    process/PUPATION_v1 \
    PUPATION_TM_freqs100 \
    ${FREQ_HALF})
echo "Pupation scan   (freqs100): ${JOB_PUP_SCAN}"

JOB_PUP_CONCAT=$(sbatch --parsable \
    --dependency=afterok:${JOB_PUP_SCAN} \
    -A tdlong_lab -p standard \
    --wrap="bash scripts/concat_Chromosome_Scans.Andreas.sh process/PUPATION_v1/PUPATION_TM_freqs100 && echo done")
echo "Pupation concat (freqs100): ${JOB_PUP_CONCAT}"

echo ""
echo "Monitor: squeue -u \$USER"
