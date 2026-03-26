#!/bin/bash
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --job-name=slide_plots
#SBATCH --output=/dfs7/adl/tdlong/fly_pool/XQTL2/slurm-slide_plots.out
#SBATCH --error=/dfs7/adl/tdlong/fly_pool/XQTL2/slurm-slide_plots.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

###############################################################################
# submit_slide_plots.sh
#
# Generates slide-ready 5-panel Manhattan PNGs for any scans that are complete.
# Skips any experiment whose pseudoscan file doesn't exist yet.
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
#
# sbatch scripts/submit_slide_plots.sh
#
# Download all PNGs when done:
#   scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/figures/slides/*.png .
###############################################################################

module load R/4.2.2

OUTDIR="figures/slides"
rm -rf "${OUTDIR}"
mkdir -p "${OUTDIR}"

run_config() {
    local label="$1"
    local config="$2"
    local check_file="$3"
    if [ ! -f "${check_file}" ]; then
        echo "SKIP ${label}: ${check_file} not found"
        return
    fi
    echo "PLOT ${label} -> configs/${config}"
    Rscript "configs/${config}"
}

run_config "malathion" "malathion.R" \
    "process/XQTL_talk/MALATHION/MALATHION_freqs250/MALATHION_freqs250.pseudoscan.txt"

run_config "zinc2" "zinc2.R" \
    "process/XQTL_talk/ZINC2/ZINC2_M_freqs250/ZINC2_M_freqs250.pseudoscan.txt"

run_config "pupation" "pupation.R" \
    "process/XQTL_talk/PUPATION/PUPATION_TB_freqs250/PUPATION_TB_freqs250.pseudoscan.txt"

run_config "aging" "aging.R" \
    "process/XQTL_talk/AGING/AGE_SY20_M_freqs250/AGE_SY20_M_freqs250.pseudoscan.txt"

echo ""
echo "Done. PNGs in ${OUTDIR}/:"
ls "${OUTDIR}/"
