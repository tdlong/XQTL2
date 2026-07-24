#!/bin/bash
#SBATCH --job-name=catalog_gather
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:30:00

# catalog_gather.sh — assemble the per-chromosome catalog pieces into one catalog.
#
# Part of the PROPOSED founder-catalog caller; a worker driven by build_catalog.sh.
# Reads <catdir>/work/catalog.<chr>.bed (written by the catalog_build.sh array) and
# writes <catdir>/catalog.tsv.gz, plus <catdir>/catalog.stats.txt summing the
# per-chromosome per-rule tallies. No reuse guard: running it rebuilds the catalog.
#
# The catalog is CHROM<TAB>POS<TAB>REF,ALT (REF/ALT comma-joined) — the form both
# `bcftools mpileup -T` and `bcftools call -m -C alleles -T` accept (catalog_count.sh).
#
# Usage:
#   sbatch catalog_gather.sh <catalog dir>

set -euo pipefail

catdir=$1
[[ -z "${catdir:-}" ]] && { echo "Usage: catalog_gather.sh <catalog dir>" >&2; exit 1; }

out="${catdir}/catalog.tsv.gz"
module load bcftools/1.21

shopt -s nullglob
pieces=( "${catdir}"/work/catalog.chr*.bed )
[[ ${#pieces[@]} -ge 1 ]] || { echo "Error: no work/catalog.chr*.bed in $catdir" >&2; exit 1; }

# pieces are CHROM POS0 POS REF ALT -> CHROM POS REF,ALT
cat "${pieces[@]}" \
  | awk -v OFS='\t' '{print $1, $3, $4","$5}' \
  | sort -k1,1 -k2,2n \
  | bgzip > "$out"
tabix -s1 -b2 -e2 "$out"

# Sum the per-rule tallies across chromosomes into one catalog.stats.txt.
awk '
  { for (i = 1; i <= NF; i++) { split($i, kv, "="); t[kv[1]] += kv[2] } }
  END {
    printf "catalog SNP rules (summed over chromosomes):\n"
    printf "  candidate biallelic SNPs : %d\n", t["candidates"]
    printf "  dropped, founder depth   : %d\n", t["depth_fail"]
    printf "  dropped, not near-fixed  : %d\n", t["notfixed_fail"]
    printf "  dropped, not segregating : %d\n", t["notseg_fail"]
    printf "  KEPT (catalog)           : %d\n", t["kept"]
  }' "${catdir}"/work/stats.chr*.txt > "${catdir}/catalog.stats.txt"

echo "catalog: $(zcat "$out" | wc -l) sites -> $out"
cat "${catdir}/catalog.stats.txt"
