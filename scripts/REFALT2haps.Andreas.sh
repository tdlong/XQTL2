#!/bin/bash
#SBATCH --job-name=RefAlt2hap
#SBATCH -A tdlong_lab        ## account to charge 
#SBATCH --mem-per-cpu=10G
#SBATCH -p highmem
#SBATCH --cpus-per-task=1 
#SBATCH --array=1-5

module load R/4.2.2

parfile=$1
mydir=$2

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

Rscript scripts/REFALT2haps.Andreas.R $mychr $parfile $mydir

