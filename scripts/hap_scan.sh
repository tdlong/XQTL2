#!/bin/bash
#SBATCH --job-name=hap_scan
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=4:00:00
#SBATCH --array=1-5

module load R/4.2.2

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

while [[ $# -gt 0 ]]; do
  case $1 in
    --rfile)  RFILE="$2";  shift 2 ;;
    --dir)    DIR="$2";    shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

Rscript "$(dirname $(readlink -f $0))/hap_scan.R" \
    --chr     "${mychr}"  \
    --dir     "${DIR}/${OUTDIR}" \
    --outdir  "${OUTDIR}" \
    --rfile   "${RFILE}"
