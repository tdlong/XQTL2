#!/bin/bash
#SBATCH --job-name=consolidateRefAlt
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-5
#SBATCH --time=1:00:00

# Consolidate per-sample RefAlt files into a single wide RefAlt table.
# One array task per chromosome. Run after all bam2bcf2REFALT_fromsites.sh jobs finish.
#
# Usage: sbatch consolidate_refalt_fromsites.sh <sample_dir> <sites_dir>
#   sample_dir : directory containing RefAlt.<samplename>.<chr>.txt files
#   sites_dir  : directory containing founder_sites.<chr>.vcf.gz

module load R/4.2.2

sample_dir=$1
sites_dir=$2

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

Rscript scripts_oneoffs/consolidate_refalt_fromsites.R \
  $sample_dir \
  $mychr \
  ${sites_dir}/founder_sites.${mychr}.vcf.gz \
  ${sample_dir}/RefAlt.${mychr}.txt
