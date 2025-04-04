#!/bin/bash
#SBATCH --job-name=pseudoscan
#SBATCH -A tdlong_lab        ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --cpus-per-task=1 
#SBATCH --array=1-5

module load R/4.2.2

Rfile=$1
mydir=$2
myoutdir=$3

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

Rscript scripts/haps2scan.Apr2025.R $mychr $Rfile $mydir $myoutdir

