#!/bin/bash
# run_haps.sh — submit the haplotype estimation step and print the job ID
#
# Submits REFALT2haps for all 5 chromosomes as a SLURM array job.
# Prints the job ID to stdout so callers can condition downstream jobs on it.
#
# Usage:
#   JID=$(bash pipeline/scripts/run_haps.sh \
#       --parfile helpfiles/<project>/hap_params.R \
#       --dir     process/<project>)
#
#   bash pipeline/scripts/run_scan.sh --after $JID ...

set -e

MEM_PER_CPU=10G
CPUS_PER_TASK=1
PARTITION=highmem
ACCOUNT=tdlong_lab
TIME=1-00:00:00     # matches REFALT2haps.sh's own header; large chromosomes (chr3R) exceed 4h

while [[ $# -gt 0 ]]; do
  case $1 in
    --parfile)      PARFILE="$2";       shift 2 ;;
    --dir)          DIR="$2";           shift 2 ;;
    --after)        AFTER="$2";         shift 2 ;;
    --mem-per-cpu)  MEM_PER_CPU="$2";   shift 2 ;;
    --time)         TIME="$2";          shift 2 ;;
    -p|--partition) PARTITION="$2";     shift 2 ;;
    -A|--account)   ACCOUNT="$2";       shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -z "$PARFILE" ]] && { echo "Error: --parfile required" >&2; exit 1; }
[[ -z "$DIR" ]]     && { echo "Error: --dir required" >&2; exit 1; }

DEP=""
[[ -n "$AFTER" ]] && DEP="--dependency=afterok:${AFTER}"

sbatch --parsable ${DEP} \
    -A ${ACCOUNT} -p ${PARTITION} \
    --cpus-per-task=${CPUS_PER_TASK} --mem-per-cpu=${MEM_PER_CPU} \
    --time=${TIME} --array=1-5 \
    pipeline/scripts/REFALT2haps.sh \
    --parfile "${PARFILE}" --dir "${DIR}" \
    | cut -d_ -f1
