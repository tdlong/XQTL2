#!/bin/bash
# compare_refalt.sh — check that two sets of RefAlt.<chr>.txt tables agree.
#
# Intended to validate a change to the SNP caller (e.g. the tiled/scatter
# bam2bcf2REFALT.sh vs the old one-task-per-chromosome version): run both into
# separate --dir directories, then diff. For a correct change the per-chromosome
# tables are byte-identical, so `diff` is empty.
#
# For each chromosome it reports one of:
#   IDENTICAL   — files match exactly (including line order)
#   REORDERED   — same set of lines, different order (sorted diff is empty)
#   DIFFERS     — content differs; prints line counts and the first few diffs
#   MISSING     — a file is absent in one directory
#
# Exits non-zero if any chromosome is DIFFERS or MISSING (REORDERED is a warning
# — investigate, but it is not a content mismatch).
#
# Usage:
#   compare_refalt.sh <dirA> <dirB> [chr1,chr2,...]

set -u

DIRA=$1
DIRB=$2
CHRS=${3:-chrX,chr2L,chr2R,chr3L,chr3R}

[[ -z "${DIRA:-}" || -z "${DIRB:-}" ]] && {
  echo "Usage: compare_refalt.sh <dirA> <dirB> [chr1,chr2,...]" >&2
  exit 2
}

status=0
IFS=',' read -r -a chrs <<< "$CHRS"

for chr in "${chrs[@]}"; do
  a="${DIRA}/RefAlt.${chr}.txt"
  b="${DIRB}/RefAlt.${chr}.txt"

  if [[ ! -f "$a" || ! -f "$b" ]]; then
    printf '%-8s MISSING   (%s%s)\n' "$chr" \
      "$([[ -f "$a" ]] || echo "A absent ")" "$([[ -f "$b" ]] || echo "B absent")"
    status=1
    continue
  fi

  na=$(wc -l < "$a" | tr -d ' ')
  nb=$(wc -l < "$b" | tr -d ' ')

  if diff -q "$a" "$b" >/dev/null; then
    printf '%-8s IDENTICAL (%s lines)\n' "$chr" "$na"
  elif diff -q <(sort "$a") <(sort "$b") >/dev/null; then
    printf '%-8s REORDERED (A=%s B=%s lines; same set, different order)\n' "$chr" "$na" "$nb"
  else
    printf '%-8s DIFFERS   (A=%s B=%s lines)\n' "$chr" "$na" "$nb"
    echo "  first differences (< A  > B):"
    diff "$a" "$b" | grep -E '^[<>]' | head -6 | sed 's/^/    /'
    # lines present in only one side (content, order-independent)
    only=$(comm -3 <(sort "$a") <(sort "$b") | grep -c .)
    echo "  lines unique to one side (sorted): $only"
    status=1
  fi
done

exit $status
