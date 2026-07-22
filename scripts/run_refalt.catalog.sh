#!/bin/bash
# run_refalt.catalog.sh — PROPOSED parallel REFALT pipeline (under evaluation).
#
# Alternative to run_refalt.sh. Instead of jointly calling all BAMs and filtering
# on QUAL, it:
#   1. catalog_build.sh  (array 1-5, ONE CHROMOSOME per task) — the slow founder
#      calling, parallelized by chromosome; writes catalog.<chr>.bed pieces.
#   2. catalog_gather.sh — concatenates the pieces into one catalog.tsv.gz.
#   3. catalog_count.sh  (array 1-<#BAMs>, ONE SAMPLE per task) — counts each BAM
#      against the fixed catalog; already-counted samples skip themselves.
#   4. catalog_merge.R   — merges per-sample counts into drop-in RefAlt.<chr>.txt.
# The validated run_refalt.sh / bam2bcf2REFALT.sh are untouched.
#
# --dir is PERSISTENT project state. The catalog is built once (parallel over
# chromosomes) and reused; each sample's counts are written once. To ADD samples,
# append their BAMs to bam_list.txt and rerun the SAME command: the founders are
# not recalled and prior samples are not recounted — only the new BAMs are
# counted, then everything is re-merged.
#
# Run this into a SEPARATE --dir and compare against a validated run with
# compare_refalt_calls.R (the two callsets are deliberately NOT identical).
#
# Prints the merge job ID to stdout (its output is RefAlt.<chr>.txt).
#
# Usage:
#   JID=$(bash pipeline/scripts/run_refalt.catalog.sh \
#       --bamlist helpfiles/<project>/bam_list.txt \
#       --parfile helpfiles/<project>/hap_params.R \
#       --dir     process/<project>_catalog)

set -e

PARTITION=standard
ACCOUNT=tdlong_lab

while [[ $# -gt 0 ]]; do
  case $1 in
    --bamlist)      BAMLIST="$2";   shift 2 ;;
    --parfile)      PARFILE="$2";   shift 2 ;;
    --founders)     FOUNDERS="$2";  shift 2 ;;
    --dir)          DIR="$2";       shift 2 ;;
    --after)        AFTER="$2";     shift 2 ;;
    -p|--partition) PARTITION="$2"; shift 2 ;;
    -A|--account)   ACCOUNT="$2";   shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -z "${BAMLIST:-}" ]] && { echo "Error: --bamlist required" >&2; exit 1; }
[[ -z "${DIR:-}" ]]     && { echo "Error: --dir required" >&2; exit 1; }

DEP=""
[[ -n "${AFTER:-}" ]] && DEP="--dependency=afterok:${AFTER}"

mkdir -p "${DIR}"
NBAM=$(grep -cve '^[[:space:]]*$' "$BAMLIST")
[[ "$NBAM" -ge 1 ]] || { echo "Error: no BAMs in $BAMLIST" >&2; exit 1; }

# Build args: founders come from --parfile (+ --bamlist, resolved by SM tag) or
# from a direct --founders BAM list.
if [[ -n "${FOUNDERS:-}" ]]; then
  BUILD_FOUNDERS="--founders ${FOUNDERS}"
else
  [[ -z "${PARFILE:-}" ]] && { echo "Error: --parfile (or --founders) required" >&2; exit 1; }
  BUILD_FOUNDERS="--parfile ${PARFILE} --bamlist ${BAMLIST}"
fi

# 1. Build the founder catalog, parallelized by chromosome (array 1-5).
#    Per-chr pieces self-reuse, so a rerun to add samples costs nothing here.
JID_BUILD=$(sbatch --parsable ${DEP} -A ${ACCOUNT} -p ${PARTITION} \
    --array=1-5 \
    pipeline/scripts/catalog_build.sh ${BUILD_FOUNDERS} --dir "${DIR}" \
    | cut -d_ -f1)

# 2. Gather the per-chromosome pieces into one genome-wide catalog.tsv.gz.
JID_CAT=$(sbatch --parsable --dependency=afterok:${JID_BUILD} \
    -A ${ACCOUNT} -p ${PARTITION} \
    pipeline/scripts/catalog_gather.sh "${DIR}")

# 3. Count every BAM against the catalog (array, one SAMPLE per task).
#    Already-counted samples skip themselves, so this only does new work.
JID_COUNT=$(sbatch --parsable --dependency=afterok:${JID_CAT} \
    -A ${ACCOUNT} -p ${PARTITION} --array=1-${NBAM} \
    pipeline/scripts/catalog_count.sh \
    "${BAMLIST}" "${DIR}" \
    | cut -d_ -f1)

# 4. Merge per-sample counts into drop-in RefAlt.<chr>.txt. This is the
#    memory-heavy step (a genome-wide join over all samples), so it goes to
#    highmem like REFALT2haps. 2 x 10G = 20G; standard caps at 6G/core and
#    highmem at 10G/core. Merge memory scales with samples x sites — tune from
#    seff once a real experiment has run.
JID_MERGE=$(sbatch --parsable --dependency=afterok:${JID_COUNT} \
    -A ${ACCOUNT} -p highmem --cpus-per-task=2 --mem-per-cpu=10G --time=02:00:00 \
    --wrap="module load R/4.2.2; Rscript pipeline/scripts/catalog_merge.R ${DIR}")

echo "${JID_MERGE}"
