#!/bin/bash
#SBATCH --job-name=callSNPs_tiled
#SBATCH -A tdlong_lab        ## account to charge
#SBATCH -p standard          ## partition/queue name
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1-00:00:00

# bam2bcf2REFALT.tiled.sh — PROPOSED alternative to bam2bcf2REFALT.sh.
#
# This is NOT the validated caller. It is a candidate replacement under
# validation: it must produce byte-identical RefAlt.<chr>.txt tables to
# bam2bcf2REFALT.sh (verify with compare_refalt.sh) before it replaces it.
#
# bcftools mpileup/call is single-threaded, so the validated one-array-task-per-
# chromosome layout scales badly as datasets grow. Here the genome is tiled into
# ~5 Mb windows (see make_tiles.sh) and each array task calls one tile, cutting
# per-task wall-clock by roughly the number of tiles per chromosome.
#
# mpileup runs on a PADDED region so positions at the tile boundary get full
# read support and realignment context; only positions inside the tile CORE are
# written out. Cores tile each chromosome with no overlap, so the per-tile RefAlt
# tables reassemble by plain concatenation (see reassemble_refalt.sh).
#
# The per-tile BCF is a throwaway intermediate (nothing downstream reads it) so
# it is deleted after the RefAlt table is built.
#
# Usage (normally via run_refalt.tiled.sh, which sizes --array to the tiling):
#   sbatch --array=1-N bam2bcf2REFALT.tiled.sh <bam_list> <output_dir> [window_bp] [pad_bp]

module load bcftools/1.21

ref="pipeline/ref/dm6.fa"
# passed from command line
bams=$1
output=$2
W=${3:-5000000}     # tile core size (bp)
PAD=${4:-10000}     # padding added to each side for calling (bp)

# Select this task's tile. make_tiles.sh is deterministic, so this row matches
# the one the wrapper used to size the array.
row=$(bash pipeline/scripts/make_tiles.sh "${ref}.fai" "$W" "$PAD" \
        | awk -v i="$SLURM_ARRAY_TASK_ID" '$1 == i')
[[ -z "$row" ]] && { echo "No tile for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
read -r idx chr cstart cend pstart pend <<< "$row"

reg="${chr}:${pstart}-${pend}"          # PADDED region: call here
tmpbcf="${output}/calls.tile_${idx}.bcf"
out="${output}/RefAlt.tile_${idx}.txt"

# Call on the padded region. -r (regions) uses the BAM index to fetch only reads
# overlapping the tile — far less I/O than the validated caller's -t (targets),
# which streamed the whole chromosome.
bcftools mpileup -I -d 1000 -r "$reg" -a "FORMAT/AD,FORMAT/DP" -f "$ref" -b "$bams" \
  | bcftools call -mv -Ob > "$tmpbcf"

# Header (identical across all tiles: same sample order from the same bam list).
echo -ne "CHROM\tPOS" > "$out"
bcftools query -l "$tmpbcf" | awk '{printf("\tREF_%s\tALT_%s",$1,$1)}' >> "$out"
echo -ne "\n" >> "$out"

# Body, filtered to the tile CORE [cstart,cend] so tiles never overlap.
bcftools view -m2 -M2 -v snps -i 'QUAL>59' "$tmpbcf" \
  | bcftools query -f'%CHROM %POS [ %AD{0} %AD{1}] [%GT]\n' \
  | grep -v '\.' | awk 'NF-=1' \
  | awk -v s="$cstart" -v e="$cend" '$2 >= s && $2 <= e' >> "$out"

# The BCF is not consumed downstream — drop it.
rm -f "$tmpbcf"
