#!/bin/bash
# run_refalt.catalog.sh — PROPOSED parallel REFALT pipeline (under evaluation).
#
# Alternative to run_refalt.sh. Instead of jointly calling all BAMs and filtering
# on QUAL, it: (1) builds a SNP catalog from the founders once, (2) counts every
# BAM against that fixed catalog as an array (one element per sample; adding a
# sample later is one more element), (3) merges into drop-in RefAlt.<chr>.txt.
# The validated run_refalt.sh / bam2bcf2REFALT.sh are untouched.
#
# Run this into a SEPARATE --dir and compare against a validated run with
# compare_refalt_calls.R (the two callsets are deliberately NOT identical).
#
# The founder set is read from the project config (--parfile hap_params.R) and
# their BAMs resolved from --bamlist by SM tag; the same bam list is counted.
#
# --dir is PERSISTENT project state. The catalog is built once (or supplied via
# --catalog) and reused; each sample's counts are written once. To ADD samples,
# append their BAMs to bam_list.txt and rerun the SAME command: the founders are
# not recalled and prior samples are not recounted — only the new BAMs are
# counted, then everything is re-merged.
#
# Prints the merge job ID to stdout (its output is RefAlt.<chr>.txt).
#
# Usage:
#   JID=$(bash pipeline/scripts/run_refalt.catalog.sh \
#       --bamlist helpfiles/<project>/bam_list.txt \
#       --parfile helpfiles/<project>/hap_params.R \
#       --dir     process/<project>_catalog)
#
#   # standard design with a precalled (shipped) catalog — no --parfile needed:
#   #   ... --bamlist ... --catalog pipeline/helpfiles/catalog_Bpop.tsv.gz --dir ...

set -e

PARTITION=standard
ACCOUNT=tdlong_lab

while [[ $# -gt 0 ]]; do
  case $1 in
    --bamlist)      BAMLIST="$2";   shift 2 ;;
    --parfile)      PARFILE="$2";   shift 2 ;;
    --dir)          DIR="$2";       shift 2 ;;
    --catalog)      CATALOG="$2";   shift 2 ;;
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

dircat="${DIR}/catalog.tsv.gz"

# 1. Catalog. The catalog is per-project persistent state, built ONCE:
#    - --catalog <file>  : use a precalled catalog (e.g. a shipped A/B-pop one) and skip building.
#    - already in --dir  : reuse it (a rerun to add samples must NOT recall founders).
#    - otherwise         : build it from the founders named in --parfile.
JID_CAT=""
if [[ -n "${CATALOG:-}" && ! -s "$dircat" ]]; then
  [[ -s "${CATALOG}.tbi" ]] || { echo "Error: ${CATALOG}.tbi missing (bgzipped + tabixed catalog required)" >&2; exit 1; }
  cp "$CATALOG" "$dircat"; cp "${CATALOG}.tbi" "${dircat}.tbi"
  echo "using precalled catalog: $CATALOG -> $dircat (build skipped)"
elif [[ -s "$dircat" ]]; then
  echo "catalog already present in $DIR; reusing (founders not recalled)"
else
  [[ -z "${PARFILE:-}" ]] && { echo "Error: --parfile required to build a catalog (or pass --catalog)" >&2; exit 1; }
  JID_CAT=$(sbatch --parsable ${DEP} -A ${ACCOUNT} -p ${PARTITION} \
      pipeline/scripts/catalog_build.sh \
      --parfile "${PARFILE}" --bamlist "${BAMLIST}" --dir "${DIR}" \
      | cut -d_ -f1)
fi

# 2. Count every BAM against the catalog (array). Already-counted samples
#    (founders + prior samples) skip themselves, so this only does new work.
#    Depend on the build job if we submitted one; else on --after (if any).
if [[ -n "$JID_CAT" ]]; then COUNT_DEP="--dependency=afterok:${JID_CAT}"; else COUNT_DEP="${DEP}"; fi
JID_COUNT=$(sbatch --parsable ${COUNT_DEP} \
    -A ${ACCOUNT} -p ${PARTITION} --array=1-${NBAM} \
    pipeline/scripts/catalog_count.sh \
    "${BAMLIST}" "${DIR}" \
    | cut -d_ -f1)

# 3. Merge per-sample counts into drop-in RefAlt.<chr>.txt (one job).
JID_MERGE=$(sbatch --parsable --dependency=afterok:${JID_COUNT} \
    -A ${ACCOUNT} -p ${PARTITION} --cpus-per-task=1 --mem-per-cpu=16G --time=02:00:00 \
    --wrap="module load R; Rscript pipeline/scripts/catalog_merge.R ${DIR}")

echo "${JID_MERGE}"
