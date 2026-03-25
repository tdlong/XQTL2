#!/bin/bash
###############################################################################
# run_all_scans_125.sh
#
# Same as run_all_scans.sh but with FREQ_SMOOTH_HALF halved → ~125 kb smoothing.
# step=5000 experiments: FREQ_SMOOTH_HALF=25 (25 × 5000 = 125 kb)
# Pupation (step=10000): FREQ_SMOOTH_HALF=13 (13 × 10000 = 130 kb)
#
# Output scan dirs: *_freqs125/
# Output PNGs:      figures/slides/*_125.png
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
# sbatch scripts/run_all_scans_125.sh
###############################################################################
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --job-name=run_125
#SBATCH --output=/dfs7/adl/tdlong/fly_pool/XQTL2/slurm-run_125.out
#SBATCH --error=/dfs7/adl/tdlong/fly_pool/XQTL2/slurm-run_125.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -e

SARAH="/dfs7/adl/sruckman/XQTL/XQTL2"
BASE="process/XQTL_talk"

# ── Helper: submit scan + concat; print concat job ID to stdout ───────────────
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

ID_ZM=$(submit_scan ZINC2_M \
    "${BASE}/ZINC2" \
    helpfiles/ZINC2/Zinc2.test.M.txt \
    ZINC2_M_freqs125 25)

ID_ZF=$(submit_scan ZINC2_F \
    "${BASE}/ZINC2" \
    helpfiles/ZINC2/Zinc2.test.F.txt \
    ZINC2_F_freqs125 25)

ID_AG=$(submit_scan AGING \
    "${BASE}/AGING" \
    helpfiles/AGE_SY/AGE_SY20_M.test.txt \
    AGE_SY20_M_freqs125 25)

ID_PU=$(submit_scan PUPATION \
    "${BASE}/PUPATION" \
    "${SARAH}/helpfiles/TBparameter.txt" \
    PUPATION_TB_freqs125 13)

ID_MA=$(submit_scan MALATHION \
    "${BASE}/MALATHION" \
    helpfiles/MALATHION/malathion.test.txt \
    MALATHION_freqs125 25)

# ── Submit plots once all concats finish ──────────────────────────────────────

JOB_PLOTS=$(sbatch --parsable \
    --dependency=afterok:${ID_ZM}:${ID_ZF}:${ID_AG}:${ID_PU}:${ID_MA} \
    -A tdlong_lab -p standard \
    --job-name=plots_125 \
    --output=/dfs7/adl/tdlong/fly_pool/XQTL2/slurm-plots_125.out \
    --error=/dfs7/adl/tdlong/fly_pool/XQTL2/slurm-plots_125.err \
    --time=00:30:00 --mem=8G \
    --wrap="module load R/4.2.2
            mkdir -p figures/slides
            Rscript configs/zinc2_125.R
            Rscript configs/aging_125.R
            Rscript configs/pupation_125.R
            Rscript configs/malathion_125.R
            echo 'Done. 125kb plots complete.'")
echo "plots job: ${JOB_PLOTS}"

echo ""
echo "All 125kb jobs submitted. Monitor: squeue -u \$USER"
