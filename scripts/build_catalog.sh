#!/bin/bash
# build_catalog.sh — build a founder SNP catalog (the "database").
#
# Stage 1 of the DEFAULT SNP caller (the founder-catalog caller); call_samples.sh
# is stage 2. The superseded joint QUAL caller is in scripts/legacy/.
# A catalog is defined by its FOUNDER SET; build one per population and point
# call_samples.sh at whichever you need. Building is a deliberate, explicit act,
# SEPARATE from calling samples: running this (re)builds the catalog and OVERWRITES
# whatever is in --out. There is no reuse/staleness check.
#
# Two stages, so that thresholds are a fast re-cut, not a rebuild (XQTL2 #21):
#   catalog_build + catalog_gather  -> catalog.annot.tsv.gz  (ALL candidate SNPs +
#         per-founder AD + distance to nearest founder indel) -- the slow part,
#         parallelized per chromosome, done ONCE.
#   catalog_filter                  -> catalog.tsv.gz + catalog.stats.txt  (applies
#         --min-dp/--maxaf/--snpgap/--exempt-founders). Re-run catalog_filter.sh
#         alone to re-cut thresholds in seconds, no founder recall.
#
# Writes into --out (the catalog dir):
#   catalog.annot.tsv.gz   annotated source (re-cuttable)
#   catalog.tsv.gz (+.tbi)  the -T positions file under the chosen thresholds
#   catalog.stats.txt       per-rule SNP tally
#   catalog.founders.txt    founder order of the AD columns
#   founders.bams.txt       the founder set that defined it
#   work/                   per-chromosome build intermediates (removed unless --keep-work)
#
# Usage:
#   bash pipeline/scripts/build_catalog.sh \
#       --founders pipeline/helpfiles/B_founders.bams.txt --out process/<project>/Catalog
#       [--min-dp 10] [--maxaf 0.03] [--snpgap 25] [--exempt-founders B5:chr2L] [--keep-work]

set -e

PARTITION=standard
ACCOUNT=tdlong_lab
KEEP_WORK=0
FILT=""   # threshold args passed through to catalog_filter.sh

while [[ $# -gt 0 ]]; do
  case $1 in
    --founders)        FOUNDERS="$2";               shift 2 ;;
    --out)             CATDIR="$2";                 shift 2 ;;
    --min-dp)          FILT="$FILT --min-dp $2";    shift 2 ;;
    --maxaf)           FILT="$FILT --maxaf $2";     shift 2 ;;
    --snpgap)          FILT="$FILT --snpgap $2";    shift 2 ;;
    --exempt-founders) FILT="$FILT --exempt-founders $2"; shift 2 ;;
    --keep-work)       KEEP_WORK=1;                 shift 1 ;;
    -p|--partition)    PARTITION="$2";              shift 2 ;;
    -A|--account)      ACCOUNT="$2";                shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -z "${FOUNDERS:-}" ]] && { echo "Error: --founders required" >&2; exit 1; }
[[ -z "${CATDIR:-}" ]]   && { echo "Error: --out (catalog dir) required" >&2; exit 1; }
[[ -f "$FOUNDERS" ]]     || { echo "Error: founder list not found: $FOUNDERS" >&2; exit 1; }

mkdir -p "${CATDIR}"
grep -ve '^[[:space:]]*$' "$FOUNDERS" > "${CATDIR}/founders.bams.txt"
echo "building catalog in ${CATDIR} from $(grep -cve '^[[:space:]]*$' "$FOUNDERS") founder BAM(s)"

# 1. Per-chromosome founder calling + annotation (array 1-5) -> work/annot.<chr>.txt
JID_BUILD=$(sbatch --parsable -A ${ACCOUNT} -p ${PARTITION} --array=1-5 \
    pipeline/scripts/catalog_build.sh --founders "${FOUNDERS}" --catdir "${CATDIR}" \
    | cut -d_ -f1)

# 2. Assemble the annotated catalog (catalog.annot.tsv.gz).
JID_GATHER=$(sbatch --parsable --dependency=afterok:${JID_BUILD} \
    -A ${ACCOUNT} -p ${PARTITION} \
    pipeline/scripts/catalog_gather.sh "${CATDIR}")

# 3. Apply thresholds -> catalog.tsv.gz + catalog.stats.txt.
JID_FILTER=$(sbatch --parsable --dependency=afterok:${JID_GATHER} \
    -A ${ACCOUNT} -p ${PARTITION} \
    pipeline/scripts/catalog_filter.sh --catdir "${CATDIR}" ${FILT})

# 4. Remove the (large) build intermediates unless asked to keep them.
LAST=${JID_FILTER}
if [[ "${KEEP_WORK}" != "1" ]]; then
  LAST=$(sbatch --parsable --dependency=afterok:${JID_FILTER} \
      -A ${ACCOUNT} -p ${PARTITION} --time=00:10:00 \
      --wrap="rm -rf ${CATDIR}/work")
fi

echo "${LAST}"
