#!/bin/bash
#SBATCH --job-name=catalog_gather
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:30:00

# catalog_gather.sh — concatenate the per-chromosome catalog pieces
# (catalog.<chr>.bed, written by the catalog_build.sh array) into one
# genome-wide, tabix-indexed catalog.tsv.gz used by the per-sample count step.
#
# Runs after the build array (afterok). Reused like the pieces: if catalog.tsv.gz
# already exists it is left in place unless FORCE=1.
#
# Usage:
#   sbatch catalog_gather.sh <output_dir>

set -euo pipefail

output=$1
[[ -z "${output:-}" ]] && { echo "Usage: catalog_gather.sh <output_dir>" >&2; exit 1; }

cat="${output}/catalog.tsv.gz"
if [[ -s "$cat" && "${FORCE:-0}" != "1" ]]; then
  echo "catalog already assembled ($cat); reusing. FORCE=1 to rebuild."
  exit 0
fi

module load bcftools/1.21

shopt -s nullglob
pieces=( "${output}"/catalog.chr*.bed )
[[ ${#pieces[@]} -ge 1 ]] || { echo "Error: no catalog.chr*.bed pieces in $output" >&2; exit 1; }

# pieces are CHROM POS0 POS REF ALT. The catalog is consumed by both
# `bcftools mpileup -T` and `bcftools call -m -C alleles -T` (catalog_count.sh);
# `-C alleles` requires 3 columns CHROM<TAB>POS<TAB>REF,ALT (REF/ALT comma-joined),
# and mpileup -T accepts that form too — so emit REF,ALT joined, not as two columns.
cat "${pieces[@]}" \
  | awk -v OFS='\t' '{print $1, $3, $4","$5}' \
  | sort -k1,1 -k2,2n \
  | bgzip > "$cat"
tabix -s1 -b2 -e2 "$cat"

echo "catalog: $(zcat "$cat" | wc -l) sites from ${#pieces[@]} chromosome piece(s) -> $cat"
