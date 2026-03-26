#!/bin/bash
#SBATCH --job-name=compareRefAlt
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1:00:00

# Compare original joint-called RefAlt against the new founder-sites RefAlt.
# Run after consolidate_refalt_fromsites.sh finishes for all chromosomes.
# Output goes to stdout/SLURM log — redirect or check with: cat slurm-<jobid>.out
#
# Usage: sbatch compare_refalt_zinc2.sh <original_dir> <new_dir>
#   original_dir : directory containing original RefAlt.<chr>.txt  (e.g. process/ZINC2)
#   new_dir      : directory containing new RefAlt.<chr>.txt       (e.g. process/ZINC2_fromsites)

module load R/4.2.2

original_dir=$1
new_dir=$2

Rscript scripts_oneoffs/compare_refalt_zinc2.R $original_dir $new_dir
