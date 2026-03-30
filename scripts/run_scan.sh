#!/bin/bash
# run_scan.sh — submit the haplotype scan pipeline with one command
#
# Submits smooth_haps → hap_scan → concat with SLURM dependency chaining.
#
# Usage:
#   bash scripts/run_scan.sh \
#       --design  helpfiles/<project>/design.txt \
#       --dir     process/<project> \
#       --scan    <scan_name> \
#       --smooth         250 \
#       --mem-per-cpu    6G \
#       --cpus-per-task  2 \
#       -p               highmem \
#       --after          <jobid>

set -e

# ── Defaults ──────────────────────────────────────────────────────────────────
SMOOTH_KB=250
MEM_PER_CPU=3G
CPUS_PER_TASK=1
PARTITION=standard
ACCOUNT=tdlong_lab

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case $1 in
    --design)        DESIGN="$2";       shift 2 ;;
    --dir)           DIR="$2";          shift 2 ;;
    --scan)          SCAN="$2";         shift 2 ;;
    --smooth)        SMOOTH_KB="$2";    shift 2 ;;
    --mem-per-cpu)   MEM_PER_CPU="$2";  shift 2 ;;
    --cpus-per-task) CPUS_PER_TASK="$2"; shift 2 ;;
    -p|--partition)  PARTITION="$2";    shift 2 ;;
    -A|--account)    ACCOUNT="$2";      shift 2 ;;
    --after)         AFTER="$2";        shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# ── Validate required arguments ──────────────────────────────────────────────
missing=""
[[ -z "$DESIGN" ]] && missing="$missing --design"
[[ -z "$DIR" ]]    && missing="$missing --dir"
[[ -z "$SCAN" ]]   && missing="$missing --scan"
if [[ -n "$missing" ]]; then
    echo "Error: missing required arguments:$missing" >&2
    echo "Usage: bash scripts/run_scan.sh --design <file> --dir <dir> --scan <name> [options]" >&2
    exit 1
fi

OUTDIR=${DIR}/${SCAN}
mkdir -p "${OUTDIR}"

# ── Build dependency for smooth step ─────────────────────────────────────────
DEP_SMOOTH=""
[[ -n "$AFTER" ]] && DEP_SMOOTH="--dependency=afterok:${AFTER}"

# ── smooth haplotype frequencies ─────────────────────────────────────────────
jid_smooth=$(sbatch --parsable ${DEP_SMOOTH} \
    -A ${ACCOUNT} -p ${PARTITION} --cpus-per-task=${CPUS_PER_TASK} --mem-per-cpu=${MEM_PER_CPU} \
    --array=1-5 "$(dirname $(readlink -f $0))/smooth_haps.sh" \
    --rfile     "${DESIGN}" \
    --dir       "${DIR}" \
    --outdir    "${SCAN}" \
    --smooth-kb "${SMOOTH_KB}")
echo "smooth:   $jid_smooth"

# ── haplotype scan (Wald + H2) ───────────────────────────────────────────────
jid_hap=$(sbatch --parsable --dependency=afterok:${jid_smooth} \
    -A ${ACCOUNT} -p ${PARTITION} --cpus-per-task=${CPUS_PER_TASK} --mem-per-cpu=${MEM_PER_CPU} \
    --array=1-5 "$(dirname $(readlink -f $0))/hap_scan.sh" \
    --rfile  "${DESIGN}" \
    --dir    "${DIR}" \
    --outdir "${SCAN}")
echo "hap_scan: $jid_hap"

# ── concat chromosomes ───────────────────────────────────────────────────────
jid_concat=$(sbatch --parsable --dependency=afterok:${jid_hap} \
    -A ${ACCOUNT} -p ${PARTITION} --cpus-per-task=${CPUS_PER_TASK} --mem-per-cpu=${MEM_PER_CPU} --time=1:00:00 \
    --wrap="bash $(readlink -f $(dirname $0))/concat_scans.sh ${OUTDIR}")
echo "concat:   $jid_concat"

echo "done:     ${jid_concat}"
