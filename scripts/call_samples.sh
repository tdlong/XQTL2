#!/bin/bash
# call_samples.sh — count BAMs against a catalog and merge into RefAlt.<chr>.txt.
#
# PROPOSED founder-catalog caller (separate from the validated bam2bcf2REFALT.sh).
# Calling samples is SEPARATE from building the catalog: this never builds a
# catalog, it counts against an existing one you choose with --catalog (the
# population you are calling).
#
# Counts EXACTLY the BAMs in --bamlist (samples and whichever founders you want as
# columns) against --catalog, then merges into drop-in RefAlt.<chr>.txt. Writes
# only into <--dir>/Calls:
#   Calls/counts/<sample>.tsv.gz     per-sample REF/ALT counts
#   Calls/RefAlt.<chr>.txt           the deliverable
#
# No skip guard: it counts (overwrites) exactly the BAMs you list. To ADD samples,
# rerun with a --bamlist of just the new BAMs — their counts land next to the
# existing ones and everything is re-merged.
#
# Usage:
#   JID=$(bash pipeline/scripts/call_samples.sh \
#       --catalog process/<project>/Catalog \
#       --bamlist helpfiles/<project>/bam_list.txt \
#       --dir     process/<project>)

set -e

PARTITION=standard
ACCOUNT=tdlong_lab

while [[ $# -gt 0 ]]; do
  case $1 in
    --catalog)      CATDIR="$2";  shift 2 ;;
    --bamlist)      BAMLIST="$2"; shift 2 ;;
    --dir)          DIR="$2";     shift 2 ;;
    --after)        AFTER="$2";   shift 2 ;;
    -p|--partition) PARTITION="$2"; shift 2 ;;
    -A|--account)   ACCOUNT="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -z "${CATDIR:-}" ]]  && { echo "Error: --catalog (catalog dir) required" >&2; exit 1; }
[[ -z "${BAMLIST:-}" ]] && { echo "Error: --bamlist required" >&2; exit 1; }
[[ -z "${DIR:-}" ]]     && { echo "Error: --dir required" >&2; exit 1; }
[[ -f "${CATDIR}/catalog.tsv.gz" ]] || { echo "Error: no catalog.tsv.gz in ${CATDIR}; run build_catalog.sh first" >&2; exit 1; }

DEP=""
[[ -n "${AFTER:-}" ]] && DEP="--dependency=afterok:${AFTER}"

CALLSDIR="${DIR}/Calls"
mkdir -p "${CALLSDIR}"
NBAM=$(grep -cve '^[[:space:]]*$' "$BAMLIST")
[[ "$NBAM" -ge 1 ]] || { echo "Error: no BAMs in $BAMLIST" >&2; exit 1; }
echo "counting ${NBAM} BAM(s) against ${CATDIR}/catalog.tsv.gz -> ${CALLSDIR}"

# 1. Count each BAM against the catalog (array, one SAMPLE per task).
JID_COUNT=$(sbatch --parsable ${DEP} \
    -A ${ACCOUNT} -p ${PARTITION} --array=1-${NBAM} \
    pipeline/scripts/catalog_count.sh \
    "${BAMLIST}" "${CATDIR}" "${CALLSDIR}" \
    | cut -d_ -f1)

# 2. Merge every counts/*.tsv.gz into drop-in RefAlt.<chr>.txt.
JID_MERGE=$(sbatch --parsable --dependency=afterok:${JID_COUNT} \
    -A ${ACCOUNT} -p ${PARTITION} --cpus-per-task=2 --mem-per-cpu=6G --time=02:00:00 \
    --wrap="module load R/4.2.2; Rscript pipeline/scripts/catalog_merge.R ${CALLSDIR}")

echo "${JID_MERGE}"
