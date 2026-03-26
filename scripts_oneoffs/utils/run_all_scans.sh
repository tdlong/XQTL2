#!/bin/bash
###############################################################################
# run_all_scans.sh
#
# Runs haps2scan → concat for all 4 experiments (haps already exist),
# then submits slide plots once all concats finish.
# All output goes to process/XQTL_talk/ — nothing else is touched.
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
# sbatch scripts/run_all_scans.sh
###############################################################################
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --job-name=run_all
#SBATCH --output=/dfs7/adl/tdlong/fly_pool/XQTL2/slurm-run_all.out
#SBATCH --error=/dfs7/adl/tdlong/fly_pool/XQTL2/slurm-run_all.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -e

SARAH="/dfs7/adl/sruckman/XQTL/XQTL2"
BASE="process/XQTL_talk"

mkdir -p ${BASE}/ZINC2 ${BASE}/AGING ${BASE}/PUPATION ${BASE}/MALATHION

# ── Helper: submit scan + concat; print concat job ID to stdout ───────────────
# submit_scan <exp> <hap_dir> <testfile> <scanname> <freq_half>
# Status messages go to stderr so stdout captures only the job ID.
submit_scan() {
    local exp=$1
    local hap_dir=$2
    local testfile=$3
    local scanname=$4
    local freq_half=$5

    local JOB_SCAN=$(sbatch --parsable --array=1-5 \
        -A tdlong_lab -p standard \
        scripts/haps2scan.freqsmooth.sh \
        "${testfile}" "${hap_dir}" "${scanname}" "${freq_half}")
    echo "${exp} scan:   ${JOB_SCAN}" >&2

    local JOB_CONCAT=$(sbatch --parsable \
        --dependency=afterok:${JOB_SCAN} \
        -A tdlong_lab -p standard \
        --wrap="bash scripts/concat_Chromosome_Scans.Andreas.sh ${hap_dir}/${scanname} && echo '${exp} concat done'")
    echo "${exp} concat: ${JOB_CONCAT}" >&2

    echo "${JOB_CONCAT}"
}

# ── Submit all 5 scan chains ──────────────────────────────────────────────────
# FREQ_HALF=50 with step=5000  → ±250 kb smoothing
# FREQ_HALF=25 with step=10000 → ±250 kb smoothing (pupation)

ID_ZM=$(submit_scan ZINC2_M \
    "${BASE}/ZINC2" \
    helpfiles/ZINC2/Zinc2.test.M.txt \
    ZINC2_M_freqs250 50)

ID_ZF=$(submit_scan ZINC2_F \
    "${BASE}/ZINC2" \
    helpfiles/ZINC2/Zinc2.test.F.txt \
    ZINC2_F_freqs250 50)

ID_AG=$(submit_scan AGING \
    "${BASE}/AGING" \
    helpfiles/AGE_SY/AGE_SY20_M.test.txt \
    AGE_SY20_M_freqs250 50)

ID_PU=$(submit_scan PUPATION \
    "${BASE}/PUPATION" \
    "${SARAH}/helpfiles/TBparameter.txt" \
    PUPATION_TB_freqs250 25)

ID_MA=$(submit_scan MALATHION \
    "${BASE}/MALATHION" \
    helpfiles/MALATHION/malathion.test.txt \
    MALATHION_freqs250 50)

# ── Submit plots once all concats finish ──────────────────────────────────────

JOB_PLOTS=$(sbatch --parsable \
    --dependency=afterok:${ID_ZM}:${ID_ZF}:${ID_AG}:${ID_PU}:${ID_MA} \
    scripts/submit_slide_plots.sh)
echo "plots job: ${JOB_PLOTS}"

echo ""
echo "All jobs submitted. Monitor: squeue -u \$USER"
