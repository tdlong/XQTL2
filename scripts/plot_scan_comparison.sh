#!/bin/bash
#SBATCH --job-name=plot_scan
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0:30:00

# Run from XQTL2-dev root:
#   sbatch pipeline/scripts/plot_scan_comparison.sh \
#       --unsmoothed process/ZINC2/ZINC2_F_unsmoothed/ZINC2_F_unsmoothed.scan.txt \
#       --smoothed   process/ZINC2/ZINC2_F_v4/ZINC2_F_v4.scan.txt \
#       --chr        chr2L \
#       --out        process/ZINC2/scan_comparison_chr2L.png

module load R/4.2.2

Rscript pipeline/scripts/plot_scan_comparison.R "$@"
