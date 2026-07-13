#!/bin/bash
#SBATCH --job-name=reassemble
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00

# reassemble_refalt.sh — stitch per-tile RefAlt tables into per-chromosome ones.
#
# The scatter job (bam2bcf2REFALT.sh) writes RefAlt.tile_<idx>.txt, each holding
# only its tile's core positions. Cores tile a chromosome with no overlap and in
# genomic order, so per chromosome we take the header from the first tile and
# append the bodies (tail -n +2) of the rest. The result, RefAlt.<chr>.txt, is
# byte-compatible with the old one-file-per-chromosome output, so everything
# downstream (REFALT2haps, ...) is unchanged.
#
# Run as a dependent job (afterok) on the scatter array; see run_refalt.sh.
#
# Usage:
#   sbatch reassemble_refalt.sh <output_dir> [window_bp] [pad_bp]

set -e

ref="pipeline/ref/dm6.fa"
output=$1
W=${2:-5000000}
PAD=${3:-10000}

tiles=$(bash pipeline/scripts/make_tiles.sh "${ref}.fai" "$W" "$PAD")

# unique chromosomes, in first-seen (genomic) order
for chr in $(echo "$tiles" | awk '!seen[$2]++ {print $2}'); do
  out="${output}/RefAlt.${chr}.txt"
  idxs=$(echo "$tiles" | awk -v c="$chr" '$2 == c {print $1}' | sort -n)

  # verify every tile for this chr is present before touching the output
  missing=""
  for idx in $idxs; do
    [[ -s "${output}/RefAlt.tile_${idx}.txt" ]] || missing="$missing $idx"
  done
  [[ -n "$missing" ]] && { echo "Error: $chr missing tile(s):$missing" >&2; exit 1; }

  first=$(echo "$idxs" | head -1)
  cat "${output}/RefAlt.tile_${first}.txt" > "$out"
  for idx in $(echo "$idxs" | tail -n +2); do
    tail -n +2 "${output}/RefAlt.tile_${idx}.txt" >> "$out"
  done
  echo "wrote $out ($(echo "$idxs" | wc -l | tr -d ' ') tiles)"
done

echo "Reassembly complete. Per-tile RefAlt.tile_*.txt left in place; remove once verified."
