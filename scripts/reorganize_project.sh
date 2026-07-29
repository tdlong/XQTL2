#!/bin/bash
# reorganize_project.sh — migrate a project dir to the Calls/Haps/Scans layout.
#
# The pipeline output is organized into stage subfolders:
#   <project>/Calls/    RefAlt.<chr>.txt, per-sample counts, (legacy) calls.<chr>.bcf
#   <project>/Haps/     R.haps.<chr>.rds, R.haps.<chr>.out.rds
#   <project>/Scans/    <scan_name>/ ...
#   <project>/Catalog/  (founder-catalog caller only)
#
# Older projects have these files flat in <project>/. This moves them into the
# subfolders in place — no recalling, no recomputing. Idempotent: run it again and
# anything already in place is left alone. Print what it moves; move nothing silently.
#
# Usage:
#   bash pipeline/scripts/reorganize_project.sh process/<project>

set -euo pipefail

proj=${1:-}
[[ -n "$proj" ]]  || { echo "Usage: reorganize_project.sh process/<project>" >&2; exit 1; }
[[ -d "$proj" ]]  || { echo "Error: not a directory: $proj" >&2; exit 1; }

cd "$proj"
mkdir -p Calls Haps Scans
shopt -s nullglob
moved=0

move_to () {  # move_to <dest> <glob...>
  local dest=$1; shift
  local f
  for f in "$@"; do
    [[ -e "$f" ]] || continue
    mv "$f" "$dest"/
    echo "  $f -> $dest/"
    moved=$((moved + 1))
  done
}

# Calls: ref/alt tables, per-sample counts, and the legacy QUAL-caller BCFs.
move_to Calls RefAlt.*.txt calls.*.bcf calls.*.bcf.csi calls.*.bcf.tbi
[[ -d counts ]] && { mv counts Calls/; echo "  counts/ -> Calls/"; moved=$((moved + 1)); }

# Haps: haplotype RDS.
move_to Haps R.haps.*.rds

# Scans: every remaining top-level directory is a scan output (skip the stage dirs).
for d in */; do
  d=${d%/}
  case "$d" in
    Calls|Haps|Scans|Catalog) ;;
    *) mv "$d" Scans/; echo "  $d/ -> Scans/"; moved=$((moved + 1)) ;;
  esac
done

if [[ "$moved" -eq 0 ]]; then
  echo "$proj: nothing to move — already in Calls/Haps/Scans layout."
else
  echo "$proj: reorganized ($moved item(s) moved)."
fi
