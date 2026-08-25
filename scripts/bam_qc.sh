#!/bin/bash
###############################################################################
# bam_qc.sh — per-sample alignment QC, one samtools flagstat per BAM
#
# Companion to refalt_qc.R. RefAlt tells you a sample is thin; this tells you
# whether that is because it had few reads or because it had plenty that did not
# map — contamination, the wrong reference, adapter dimer. Different diagnosis,
# different fix.
#
# Joins to refalt_qc.txt on `sample` (the BAM basename without .bam), which is
# the same name call_samples.sh uses for the RefAlt columns.
#
# has_unmapped is reported because pct_mapped is only meaningful if unmapped
# reads are still present. A BAM that arrived via --skip-fq2bam may have had them
# stripped upstream, in which case pct_mapped reads a meaningless 100.0. This
# pipeline's own fq2bam.sh keeps them (no -F filter), so for BAMs it produced the
# number is real.
#
# Note: no duplicate metric. These libraries are not deduplicated by design — at
# 50-200x, reads sharing a start site occur legitimately, especially from Nextera
# where the transposase has insertion-site preference, so marking them as
# duplicates would delete real molecules and distort the very allele frequencies
# the pipeline measures.
#
# Usage:
#   bash pipeline/scripts/bam_qc.sh --bamlist helpfiles/<project>/bam_list.txt \
#       [--out process/<project>/Calls/bam_qc.txt]
###############################################################################
set -euo pipefail

BAMLIST=""
OUT=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --bamlist) BAMLIST="$2"; shift 2 ;;
    --out)     OUT="$2";     shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done
[[ -n "$BAMLIST" ]] || { echo "Usage: bam_qc.sh --bamlist <file> [--out <file>]" >&2; exit 1; }
[[ -f "$BAMLIST" ]] || { echo "Error: no such file: $BAMLIST" >&2; exit 1; }
[[ -n "$OUT" ]] || OUT="bam_qc.txt"

command -v samtools >/dev/null || { echo "Error: samtools not on PATH (module load samtools)" >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"
printf 'sample\tn_reads\tpct_mapped\tpct_proper_pair\thas_unmapped\n' > "$OUT"

n=0
while read -r bam; do
  [[ -n "$bam" ]] || continue
  if [[ ! -f "$bam" ]]; then
    echo "  warning: missing, skipped: $bam" >&2
    continue
  fi
  s=$(basename "$bam"); s=${s%.bam}

  # One pass; parse the fixed flagstat lines.
  fs=$(samtools flagstat "$bam")
  total=$(awk '/in total/        {print $1; exit}' <<< "$fs")
  mapped=$(awk '/ mapped \(/     {print $1; exit}' <<< "$fs")
  proper=$(awk '/properly paired/{print $1; exit}' <<< "$fs")

  if [[ "$total" -gt 0 ]]; then
    unmapped=$(( total - mapped ))
    pm=$(awk -v a="$mapped" -v b="$total" 'BEGIN{printf "%.2f", 100*a/b}')
    pp=$(awk -v a="$proper" -v b="$total" 'BEGIN{printf "%.2f", 100*a/b}')
    hu=$([[ "$unmapped" -gt 0 ]] && echo yes || echo no)
  else
    pm="NA"; pp="NA"; hu="NA"
  fi

  printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$total" "$pm" "$pp" "$hu" >> "$OUT"
  n=$((n + 1))
done < "$BAMLIST"

echo "Wrote $OUT ($n sample(s))"

# has_unmapped=no means the percentages above describe a pre-filtered BAM and
# cannot be compared with one this pipeline aligned.
if awk -F'\t' 'NR>1 && $5=="no"{c++} END{exit !(c>0)}' "$OUT"; then
  echo ""
  echo "Note: some BAMs contain no unmapped reads, so their pct_mapped is not a"
  echo "real mapping rate — those reads were removed before this pipeline saw them:"
  awk -F'\t' 'NR>1 && $5=="no"{print "  " $1}' "$OUT"
fi
