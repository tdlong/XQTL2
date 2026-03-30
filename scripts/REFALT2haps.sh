#!/bin/bash
#SBATCH --job-name=RefAlt2hap
#SBATCH -A tdlong_lab
#SBATCH --mem-per-cpu=10G
#SBATCH -p highmem
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --array=1-5

module load R/4.2.2

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --parfile) parfile="$2"; shift 2 ;;
        --dir)     mydir="$2";   shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

[[ -z "$parfile" ]] && { echo "Error: --parfile required"; exit 1; }
[[ -z "$mydir"   ]] && { echo "Error: --dir required";    exit 1; }

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

Rscript "$(dirname $(readlink -f $0))/REFALT2haps.R" $mychr $parfile $mydir
