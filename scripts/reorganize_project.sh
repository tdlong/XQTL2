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
# The whole plan is checked before anything moves, so the project either migrates
# or it does not — it is never left half-way (XQTL2 #30). Three cases per item:
#
#   destination missing              -> move
#   destination is the same file     -> already migrated, skip (this is the
#                                       top-level RefAlt.<chr>.txt symlink into
#                                       Calls/ that flat-path tooling leaves behind)
#   destination exists and differs   -> refuse the run and report every conflict,
#                                       because moving would destroy the file in
#                                       place. --force overrides and overwrites.
#
# That last case is not hypothetical: a project called with the legacy caller has
# top-level RefAlt.*, and re-calling it with call_samples.sh writes Calls/RefAlt.*.
# A bare mv would replace the new tables with the old ones and look like success.
#
# Usage:
#   bash pipeline/scripts/reorganize_project.sh process/<project> [--force]

set -euo pipefail

proj=""
force=false
while [[ $# -gt 0 ]]; do
  case $1 in
    --force) force=true; shift ;;
    -*) echo "Unknown argument: $1" >&2; exit 1 ;;
    *) [[ -z "$proj" ]] || { echo "Error: only one project directory" >&2; exit 1; }
       proj=$1; shift ;;
  esac
done
[[ -n "$proj" ]]  || { echo "Usage: reorganize_project.sh process/<project> [--force]" >&2; exit 1; }
[[ -d "$proj" ]]  || { echo "Error: not a directory: $proj" >&2; exit 1; }

cd "$proj"
shopt -s nullglob

# ── Plan ────────────────────────────────────────────────────────────────────
plan_src=(); plan_dst=()      # to move
skip_src=()                   # already migrated (same file)
conf_src=(); conf_dst=()      # destination exists and differs

consider () {  # consider <dest> <item...>
  local dest=$1; shift
  local f base target
  for f in "$@"; do
    [[ -e "$f" || -L "$f" ]] || continue
    base=$(basename "$f")
    target="$dest/$base"
    if [[ ! -e "$target" && ! -L "$target" ]]; then
      plan_src+=("$f"); plan_dst+=("$dest")
    elif [[ "$f" -ef "$target" ]]; then
      # -ef compares device+inode after following symlinks, so this catches the
      # top-level symlink whose target IS the destination — mv would abort on it.
      skip_src+=("$f")
    else
      conf_src+=("$f"); conf_dst+=("$target")
    fi
  done
}

consider Calls RefAlt.*.txt calls.*.bcf calls.*.bcf.csi calls.*.bcf.tbi
[[ -d counts ]] && consider Calls counts

consider Haps R.haps.*.rds

# Scans: every remaining top-level directory is a scan output (skip the stage dirs).
for d in */; do
  d=${d%/}
  case "$d" in
    Calls|Haps|Scans|Catalog|counts) ;;
    *) consider Scans "$d" ;;
  esac
done

# ── Refuse on conflict, before touching anything ────────────────────────────
if [[ ${#conf_src[@]} -gt 0 && "$force" == false ]]; then
  echo "Error: $proj cannot be migrated — ${#conf_src[@]} destination(s) already exist and differ:" >&2
  for i in "${!conf_src[@]}"; do
    echo "  ${conf_src[$i]}  ->  ${conf_dst[$i]}  (exists, different file)" >&2
  done
  echo "" >&2
  echo "Nothing has been moved. This usually means the project was re-called with" >&2
  echo "call_samples.sh (which writes Calls/) while the legacy flat files are still" >&2
  echo "present. The files in Calls/ are normally the ones you want; delete or move" >&2
  echo "the top-level copies aside, then re-run." >&2
  echo "Use --force to overwrite the destinations with the top-level copies." >&2
  exit 1
fi

# ── Execute ─────────────────────────────────────────────────────────────────
mkdir -p Calls Haps Scans

moved=0
for i in "${!plan_src[@]}"; do
  mv "${plan_src[$i]}" "${plan_dst[$i]}"/
  echo "  ${plan_src[$i]} -> ${plan_dst[$i]}/"
  moved=$((moved + 1))
done

if [[ "$force" == true && ${#conf_src[@]} -gt 0 ]]; then
  for i in "${!conf_src[@]}"; do
    mv -f "${conf_src[$i]}" "${conf_dst[$i]}"
    echo "  ${conf_src[$i]} -> ${conf_dst[$i]}  (overwritten, --force)"
    moved=$((moved + 1))
  done
fi

for f in "${skip_src[@]:-}"; do
  [[ -n "$f" ]] || continue
  echo "  $f already points into its stage folder — left alone"
done

if [[ "$moved" -eq 0 ]]; then
  echo "$proj: nothing to move — already in Calls/Haps/Scans layout."
else
  echo "$proj: reorganized ($moved item(s) moved)."
fi
