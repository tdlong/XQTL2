#!/bin/bash
# show_project_layout.sh
# Summarize what pipeline outputs exist for a given project.
# Run from the XQTL2 project root on the cluster.
#
# Usage:
#   bash scripts/show_project_layout.sh <project_name>
#
# Example:
#   bash scripts/show_project_layout.sh ZINC2

project=$1
[[ -z "$project" ]] && { echo "Usage: $0 <project_name>"; exit 1; }

dir="process/$project"
[[ ! -d "$dir" ]] && { echo "Directory not found: $dir"; exit 1; }

chrs=(chrX chr2L chr2R chr3L chr3R)

echo "========================================"
echo "Project: $project"
echo "Root:    $dir"
echo "========================================"

# ── Step 3: REFALT counts ────────────────────────────────────────────────────
echo ""
echo "--- REFALT counts (Step 3) ---"
for chr in "${chrs[@]}"; do
    f="$dir/RefAlt.$chr.txt"
    if [[ -f "$f" ]]; then
        size=$(du -sh "$f" | cut -f1)
        echo "  [OK]  $f  ($size)"
    else
        echo "  [--]  $f"
    fi
done

# ── Step 4: Haplotypes ───────────────────────────────────────────────────────
echo ""
echo "--- Haplotype files (Step 4) ---"
for chr in "${chrs[@]}"; do
    f_rds="$dir/R.haps.$chr.rds"
    f_out="$dir/R.haps.$chr.out.rds"
    if [[ -f "$f_rds" ]]; then
        size=$(du -sh "$f_rds" | cut -f1)
        echo "  [OK]  $f_rds  ($size)"
    else
        echo "  [--]  $f_rds"
    fi
    if [[ -f "$f_out" ]]; then
        size=$(du -sh "$f_out" | cut -f1)
        echo "  [OK]  $f_out  ($size)"
    else
        echo "  [--]  $f_out"
    fi
done

# ── Step 5/6: Scan subdirectories ────────────────────────────────────────────
echo ""
echo "--- Scan subdirectories (Steps 5-6) ---"
found_scans=0
for scan_dir in "$dir"/*/; do
    [[ ! -d "$scan_dir" ]] && continue
    scan_name=$(basename "$scan_dir")
    found_scans=1
    echo ""
    echo "  Scan: $scan_name"

    # Per-chromosome scan files (any naming)
    chr_scan_count=$(ls "$scan_dir"/${scan_name}.*scan*.chr*.txt 2>/dev/null | wc -l)
    chr_means_count=$(ls "$scan_dir"/${scan_name}.meansBySample.chr*.txt 2>/dev/null | wc -l)
    echo "    per-chr scan files : $chr_scan_count / 5"
    echo "    per-chr means files: $chr_means_count / 5"

    # Concatenated scan (any naming: .scan.txt or .pseudoscan.txt)
    concat_scan=$(ls "$scan_dir"/${scan_name}.*scan.txt 2>/dev/null | grep -v chr | head -1)
    if [[ -n "$concat_scan" ]]; then
        size=$(du -sh "$concat_scan" | cut -f1)
        echo "    [OK]  $(basename $concat_scan)  ($size)"
    else
        echo "    [--]  concatenated scan"
    fi

    # meansBySample concatenated
    f="$scan_dir/${scan_name}.meansBySample.txt"
    if [[ -f "$f" ]]; then
        size=$(du -sh "$f" | cut -f1)
        echo "    [OK]  ${scan_name}.meansBySample.txt  ($size)"
    else
        echo "    [--]  ${scan_name}.meansBySample.txt"
    fi

    # Figures
    fig_count=$(ls "$scan_dir"/*.png 2>/dev/null | wc -l)
    echo "    figures: $fig_count png file(s)"
done

[[ $found_scans -eq 0 ]] && echo "  (none found)"

echo ""
echo "========================================"
