#!/bin/bash
# build_catalog.sh — build a founder SNP catalog (the "database").
#
# PROPOSED founder-catalog caller (separate from the validated bam2bcf2REFALT.sh).
# A catalog is defined by its FOUNDER SET; build one per population and point
# call_samples.sh at whichever you need. Building is a deliberate, explicit act,
# SEPARATE from calling samples: running this (re)builds the catalog and OVERWRITES
# whatever is in --out. There is no reuse/staleness check — run it to rebuild.
#
# Writes into --out (the catalog dir):
#   catalog.tsv.gz (+ .tbi)   the catalog (CHROM POS REF,ALT)
#   founders.bams.txt         the founder set that defined it
#   catalog.stats.txt         per-rule SNP tally (candidates, dropped-by-rule, kept)
#   work/                     per-chromosome build intermediates (removed unless --keep-work)
#
# SNP rules and the --exempt-founders default (B5:chr2L) are documented in
# catalog_build.sh and the README appendix.
#
# Usage:
#   bash pipeline/scripts/build_catalog.sh \
#       --founders helpfiles/B_founders.bams.txt --out process/<project>/Catalog
#       [--exempt-founders B5:chr2L] [--min-dp 10] [--maxaf 0.03] [--snpgap 5] [--keep-work]

set -e

PARTITION=standard
ACCOUNT=tdlong_lab
KEEP_WORK=0
PASS=""   # tunables passed through to catalog_build.sh

while [[ $# -gt 0 ]]; do
  case $1 in
    --founders)        FOUNDERS="$2";              shift 2 ;;
    --out)             CATDIR="$2";                shift 2 ;;
    --exempt-founders) PASS="$PASS --exempt-founders $2"; shift 2 ;;
    --min-dp)          PASS="$PASS --min-dp $2";   shift 2 ;;
    --maxaf)           PASS="$PASS --maxaf $2";    shift 2 ;;
    --snpgap)          PASS="$PASS --snpgap $2";   shift 2 ;;
    --keep-work)       KEEP_WORK=1;                shift 1 ;;
    -p|--partition)    PARTITION="$2";             shift 2 ;;
    -A|--account)      ACCOUNT="$2";               shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -z "${FOUNDERS:-}" ]] && { echo "Error: --founders required" >&2; exit 1; }
[[ -z "${CATDIR:-}" ]]   && { echo "Error: --out (catalog dir) required" >&2; exit 1; }
[[ -f "$FOUNDERS" ]]     || { echo "Error: founder list not found: $FOUNDERS" >&2; exit 1; }

mkdir -p "${CATDIR}"
grep -ve '^[[:space:]]*$' "$FOUNDERS" > "${CATDIR}/founders.bams.txt"    # record the founder set
echo "building catalog in ${CATDIR} from $(grep -cve '^[[:space:]]*$' "$FOUNDERS") founder BAM(s)${PASS:+ (}${PASS}${PASS:+ )}"

# 1. Per-chromosome founder calling + fixation filter (array 1-5) -> work/*.bed
JID_BUILD=$(sbatch --parsable -A ${ACCOUNT} -p ${PARTITION} --array=1-5 \
    pipeline/scripts/catalog_build.sh --founders "${FOUNDERS}" --catdir "${CATDIR}" ${PASS} \
    | cut -d_ -f1)

# 2. Assemble the catalog + sum the per-rule tally.
JID_GATHER=$(sbatch --parsable --dependency=afterok:${JID_BUILD} \
    -A ${ACCOUNT} -p ${PARTITION} \
    pipeline/scripts/catalog_gather.sh "${CATDIR}")

# 3. Remove the (large) build intermediates unless asked to keep them.
LAST=${JID_GATHER}
if [[ "${KEEP_WORK}" != "1" ]]; then
  LAST=$(sbatch --parsable --dependency=afterok:${JID_GATHER} \
      -A ${ACCOUNT} -p ${PARTITION} --time=00:10:00 \
      --wrap="rm -rf ${CATDIR}/work")
fi

echo "${LAST}"
