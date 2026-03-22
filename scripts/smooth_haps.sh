#!/bin/bash
#SBATCH --job-name=smooth_haps
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2          ## R is single-threaded; 2 cores to get 12G total
#SBATCH --mem-per-cpu=6G
#SBATCH --time=2:00:00
#SBATCH --array=1-5

module load R/4.2.2

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

# Parse named arguments forwarded from the pipeline driver
while [[ $# -gt 0 ]]; do
  case $1 in
    --rfile)    RFILE="$2";     shift 2 ;;
    --dir)      DIR="$2";       shift 2 ;;
    --outdir)   OUTDIR="$2";    shift 2 ;;
    --smooth-kb) SMOOTH_KB="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

Rscript scripts/smooth_haps.R \
    --chr       "${mychr}"    \
    --dir       "${DIR}"      \
    --outdir    "${OUTDIR}"   \
    --rfile     "${RFILE}"    \
    --smooth-kb "${SMOOTH_KB}"
