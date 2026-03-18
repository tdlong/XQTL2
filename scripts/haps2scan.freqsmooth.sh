#!/bin/bash
#SBATCH --job-name=freqsmooth
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5        # all chromosomes (chrX=1 chr2L=2 chr2R=3 chr3L=4 chr3R=5)

module load R/4.2.2

Rfile=$1
mydir=$2
myoutdir=$3
FREQ_SMOOTH_HALF=$4

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

Rscript scripts/haps2scan.freqsmooth.R $mychr $Rfile $mydir $myoutdir $FREQ_SMOOTH_HALF
