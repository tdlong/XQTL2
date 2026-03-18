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
    echo "  Scan: $scan_name  ($scan_dir)"

    # Per-chromosome scan files
    chr_count=0
    for chr in "${chrs[@]}"; do
        f="$scan_dir/${scan_name}.scan.$chr.txt"
        [[ -f "$f" ]] && chr_count=$((chr_count + 1))
    done
    echo "    per-chr scan files : $chr_count / 5"

    chr_count=0
    for chr in "${chrs[@]}"; do
        f="$scan_dir/${scan_name}.meansBySample.$chr.txt"
        [[ -f "$f" ]] && chr_count=$((chr_count + 1))
    done
    echo "    per-chr means files: $chr_count / 5"

    # Concatenated outputs
    for suffix in "scan.txt" "meansBySample.txt" "5panel.Mb.png" "5panel.cM.png" "Manhattan.png"; do
        f="$scan_dir/${scan_name}.$suffix"
        if [[ -f "$f" ]]; then
            size=$(du -sh "$f" | cut -f1)
            echo "    [OK]  ${scan_name}.$suffix  ($size)"
        else
            echo "    [--]  ${scan_name}.$suffix"
        fi
    done
done

[[ $found_scans -eq 0 ]] && echo "  (none found)"

echo ""
echo "========================================"
