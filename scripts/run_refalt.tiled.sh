#!/bin/bash
# run_refalt.tiled.sh — PROPOSED alternative to run_refalt.sh (under validation).
#
# Submits the tiled/scatter SNP caller (bam2bcf2REFALT.tiled.sh, one genome tile
# per array task) and chains a reassembly job (reassemble_refalt.sh) that stitches
# the per-tile tables back into per-chromosome RefAlt.<chr>.txt files. The array
# size is derived from the tiling so it always matches the number of tiles.
#
# This is a candidate replacement for run_refalt.sh, NOT the validated path. Run
# it into a SEPARATE --dir and compare against a validated run with
# compare_refalt.sh before adopting it. It does not touch run_refalt.sh or
# bam2bcf2REFALT.sh.
#
# Prints the reassembly job ID to stdout (the reassembled per-chromosome tables
# are its output), so downstream jobs can --after it.
#
# Usage:
#   JID=$(bash pipeline/scripts/run_refalt.tiled.sh \
#       --bamlist helpfiles/<project>/bam_list.txt \
#       --dir     process/<project>_tiled)

set -e

MEM_PER_CPU=6G
CPUS_PER_TASK=2
PARTITION=standard
ACCOUNT=tdlong_lab
REF=pipeline/ref/dm6.fa
WINDOW=5000000      # tile core size (bp)
PAD=10000           # calling pad per side (bp)

while [[ $# -gt 0 ]]; do
  case $1 in
    --bamlist)       BAMLIST="$2";       shift 2 ;;
    --dir)           DIR="$2";           shift 2 ;;
    --after)         AFTER="$2";         shift 2 ;;
    --window)        WINDOW="$2";        shift 2 ;;
    --pad)           PAD="$2";           shift 2 ;;
    --ref)           REF="$2";           shift 2 ;;
    --mem-per-cpu)   MEM_PER_CPU="$2";   shift 2 ;;
    --cpus-per-task) CPUS_PER_TASK="$2"; shift 2 ;;
    -p|--partition)  PARTITION="$2";     shift 2 ;;
    -A|--account)    ACCOUNT="$2";       shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -z "$BAMLIST" ]] && { echo "Error: --bamlist required" >&2; exit 1; }
[[ -z "$DIR" ]]     && { echo "Error: --dir required" >&2; exit 1; }
[[ -f "${REF}.fai" ]] || { echo "Error: fai not found: ${REF}.fai" >&2; exit 1; }

DEP=""
[[ -n "$AFTER" ]] && DEP="--dependency=afterok:${AFTER}"

mkdir -p "${DIR}"

# Number of tiles = number of array tasks (single source of truth: make_tiles.sh).
NTILES=$(bash pipeline/scripts/make_tiles.sh "${REF}.fai" "$WINDOW" "$PAD" | wc -l | tr -d ' ')
[[ "$NTILES" -ge 1 ]] || { echo "Error: no tiles produced from ${REF}.fai" >&2; exit 1; }

# Scatter: one tile per array task.
JID_SCATTER=$(sbatch --parsable ${DEP} \
    -A ${ACCOUNT} -p ${PARTITION} \
    --cpus-per-task=${CPUS_PER_TASK} --mem-per-cpu=${MEM_PER_CPU} \
    --time=1-00:00:00 --array=1-${NTILES} \
    pipeline/scripts/bam2bcf2REFALT.tiled.sh \
    "${BAMLIST}" "${DIR}" "${WINDOW}" "${PAD}" \
    | cut -d_ -f1)

# Gather: reassemble per-chromosome tables once every tile succeeds.
JID_GATHER=$(sbatch --parsable --dependency=afterok:${JID_SCATTER} \
    -A ${ACCOUNT} -p ${PARTITION} \
    pipeline/scripts/reassemble_refalt.sh \
    "${DIR}" "${WINDOW}" "${PAD}")

# Downstream should wait on the reassembly (it produces RefAlt.<chr>.txt).
echo "${JID_GATHER}"
