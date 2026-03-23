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
#       --smooth  250 \
#       --mem     6G \
#       --after   <jobid>

set -e

# ── Defaults ──────────────────────────────────────────────────────────────────
SMOOTH_KB=250
MEM_PER_CPU=3G

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case $1 in
    --design) DESIGN="$2";    shift 2 ;;
    --dir)    DIR="$2";       shift 2 ;;
    --scan)   SCAN="$2";      shift 2 ;;
    --smooth) SMOOTH_KB="$2"; shift 2 ;;
    --mem)    MEM_PER_CPU="$2"; shift 2 ;;
    --after)  AFTER="$2";     shift 2 ;;
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
    --mem-per-cpu=${MEM_PER_CPU} \
    --array=1-5 scripts/smooth_haps.sh \
    --rfile     "${DESIGN}" \
    --dir       "${DIR}" \
    --outdir    "${SCAN}" \
    --smooth-kb "${SMOOTH_KB}")
echo "smooth:   $jid_smooth"

# ── haplotype scan (Wald + H2) ───────────────────────────────────────────────
jid_hap=$(sbatch --parsable --dependency=afterok:${jid_smooth} \
    --array=1-5 scripts/hap_scan.sh \
    --rfile  "${DESIGN}" \
    --dir    "${DIR}" \
    --outdir "${SCAN}")
echo "hap_scan: $jid_hap"

# ── concat chromosomes ───────────────────────────────────────────────────────
jid_concat=$(sbatch --parsable --dependency=afterok:${jid_hap} \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=1:00:00 \
    --wrap="bash scripts/concat_scans.sh ${OUTDIR}")
echo "concat:   $jid_concat"

echo "done:     ${jid_concat}"
