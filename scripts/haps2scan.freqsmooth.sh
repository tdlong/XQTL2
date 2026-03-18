#!/bin/bash
#SBATCH --job-name=scan_smooth
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5

module load R/4.2.2

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --rfile)          Rfile="$2";          shift 2 ;;
        --dir)            mydir="$2";          shift 2 ;;
        --outdir)         myoutdir="$2";       shift 2 ;;
        --cov-smooth-kb)  COV_SMOOTH_KB="$2";  shift 2 ;;
        --freq-smooth-kb) FREQ_SMOOTH_KB="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

[[ -z "$Rfile"         ]] && { echo "Error: --rfile required";          exit 1; }
[[ -z "$mydir"         ]] && { echo "Error: --dir required";            exit 1; }
[[ -z "$myoutdir"      ]] && { echo "Error: --outdir required";         exit 1; }
[[ -z "$COV_SMOOTH_KB" ]] && { echo "Error: --cov-smooth-kb required";  exit 1; }
[[ -z "$FREQ_SMOOTH_KB" ]] && { echo "Error: --freq-smooth-kb required"; exit 1; }

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

Rscript scripts/haps2scan.freqsmooth.R \
    $mychr $Rfile $mydir $myoutdir $COV_SMOOTH_KB $FREQ_SMOOTH_KB
