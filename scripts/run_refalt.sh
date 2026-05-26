#!/bin/bash
# run_refalt.sh — submit the REFALT allele-count step and print the job ID
#
# Submits bam2bcf2REFALT for all 5 chromosomes as a SLURM array job.
# Prints the job ID to stdout so callers can condition downstream jobs on it.
#
# Usage:
#   JID_REFALT=$(bash pipeline/scripts/run_refalt.sh \
#       --bamlist helpfiles/<project>/bam_list.txt \
#       --dir     process/<project>)
#
#   JID_HAPS=$(bash pipeline/scripts/run_haps.sh \
#       --after   $JID_REFALT \
#       --parfile helpfiles/<project>/hap_params.R \
#       --dir     process/<project>)

set -e

MEM_PER_CPU=6G
CPUS_PER_TASK=2
PARTITION=standard
ACCOUNT=tdlong_lab

while [[ $# -gt 0 ]]; do
  case $1 in
    --bamlist)      BAMLIST="$2";      shift 2 ;;
    --dir)          DIR="$2";          shift 2 ;;
    --after)        AFTER="$2";        shift 2 ;;
    --mem-per-cpu)  MEM_PER_CPU="$2";  shift 2 ;;
    --cpus-per-task) CPUS_PER_TASK="$2"; shift 2 ;;
    -p|--partition) PARTITION="$2";    shift 2 ;;
    -A|--account)   ACCOUNT="$2";      shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -z "$BAMLIST" ]] && { echo "Error: --bamlist required" >&2; exit 1; }
[[ -z "$DIR" ]]     && { echo "Error: --dir required" >&2; exit 1; }

DEP=""
[[ -n "$AFTER" ]] && DEP="--dependency=afterok:${AFTER}"

mkdir -p "${DIR}"

sbatch --parsable ${DEP} \
    -A ${ACCOUNT} -p ${PARTITION} \
    --cpus-per-task=${CPUS_PER_TASK} --mem-per-cpu=${MEM_PER_CPU} \
    --time=5-00:00:00 --array=1-5 \
    pipeline/scripts/bam2bcf2REFALT.sh \
    "${BAMLIST}" "${DIR}" \
    | cut -d_ -f1
