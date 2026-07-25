#!/bin/bash
#SBATCH --job-name=catalog_gather
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:30:00

# catalog_gather.sh — assemble the per-chromosome annotation pieces into one
# annotated catalog table.
#
# Part of the PROPOSED founder-catalog caller; a worker driven by build_catalog.sh.
# Reads <catdir>/work/annot.<chr>.txt (from the catalog_build.sh array) and writes
# <catdir>/catalog.annot.tsv.gz with a header naming the founder AD columns:
#   CHROM POS REF ALT dist_indel AD_<founder1> AD_<founder2> ...
# This is the re-cuttable source: catalog_filter.sh turns it into catalog.tsv.gz
# (the -T positions file) under whatever thresholds you choose, with no rebuild.
#
# Usage:
#   sbatch catalog_gather.sh <catalog dir>

set -euo pipefail

catdir=$1
[[ -z "${catdir:-}" ]] && { echo "Usage: catalog_gather.sh <catalog dir>" >&2; exit 1; }

module load bcftools/1.21

# Founder names in AD-column order (the sample order of any founder BCF).
shopt -s nullglob
somebcf=$(ls "${catdir}"/work/founders.calls.chr*.bcf 2>/dev/null | head -1)
[[ -n "$somebcf" ]] || { echo "Error: no work/founders.calls.chr*.bcf in $catdir" >&2; exit 1; }
mapfile -t fnames < <(bcftools query -l "$somebcf")

hdr="CHROM	POS	REF	ALT	dist_indel"
for f in "${fnames[@]}"; do hdr="${hdr}	AD_${f}"; done

pieces=( "${catdir}"/work/annot.chr*.txt )
[[ ${#pieces[@]} -ge 1 ]] || { echo "Error: no work/annot.chr*.txt in $catdir" >&2; exit 1; }

{ printf '%s\n' "$hdr"; cat "${pieces[@]}"; } | bgzip > "${catdir}/catalog.annot.tsv.gz"
printf '%s\n' "${fnames[@]}" > "${catdir}/catalog.founders.txt"

echo "annotated catalog: $(( $(zcat "${catdir}/catalog.annot.tsv.gz" | wc -l) - 1 )) candidate SNPs, ${#fnames[@]} founders -> ${catdir}/catalog.annot.tsv.gz"
