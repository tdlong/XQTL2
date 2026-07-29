#!/bin/bash
# make_tiles.sh — tile a reference into padded windows for scatter SNP calling.
#
# Reads a .fai index and emits one tab-separated line per tile:
#   idx  chr  core_start  core_end  pad_start  pad_end
#
# All coordinates are 1-based inclusive. SNPs are CALLED on the padded region
# [pad_start,pad_end] (so boundary positions get full read/realignment context)
# but only KEPT if core_start <= POS <= core_end. Cores tile the chromosome
# with no gaps and no overlap, so per-tile outputs reassemble by plain
# concatenation — no dedup, no merge.
#
# This script is the single source of truth for the tiling: the submission
# wrapper calls it to size the array (--array=1-N), and each array task calls
# it to select its own row by SLURM_ARRAY_TASK_ID. They therefore agree by
# construction.
#
# Usage:
#   make_tiles.sh <ref.fai> <window_bp> <pad_bp> [chr1,chr2,...]
#
# The optional chromosome list restricts (and orders) the tiled sequences.
# Default matches the callable euchromatin used elsewhere in the pipeline.

set -e

FAI=$1
W=$2
PAD=$3
CHRS=${4:-chrX,chr2L,chr2R,chr3L,chr3R}

[[ -z "$FAI" || -z "$W" || -z "$PAD" ]] && {
  echo "Usage: make_tiles.sh <ref.fai> <window_bp> <pad_bp> [chr1,chr2,...]" >&2
  exit 1
}
[[ -f "$FAI" ]] || { echo "Error: fai not found: $FAI" >&2; exit 1; }

awk -v W="$W" -v PAD="$PAD" -v CHRS="$CHRS" '
BEGIN { n = split(CHRS, want, ","); OFS = "\t"; idx = 0 }
{ len[$1] = $2 }                       # column 1 = name, column 2 = length
END {
  for (i = 1; i <= n; i++) {
    c = want[i]; L = len[c]
    if (L == "") continue              # requested chr absent from this fai
    for (s = 1; s <= L; s += W) {
      e = s + W - 1;      if (e > L) e = L
      ps = s - PAD;       if (ps < 1) ps = 1
      pe = e + PAD;       if (pe > L) pe = L
      idx++
      print idx, c, s, e, ps, pe
    }
  }
}' "$FAI"
