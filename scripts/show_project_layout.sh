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

    # Per-chromosome files — check both current (.scan.) and legacy (.pseudoscan.) names
    scan_count=0; pseudo_count=0
    means_count=0; means_pseudo_count=0
    for chr in "${chrs[@]}"; do
        [[ -f "$scan_dir/${scan_name}.scan.$chr.txt"       ]] && scan_count=$((scan_count+1))
        [[ -f "$scan_dir/${scan_name}.pseudoscan.$chr.txt" ]] && pseudo_count=$((pseudo_count+1))
        [[ -f "$scan_dir/${scan_name}.meansBySample.$chr.txt"       ]] && means_count=$((means_count+1))
        [[ -f "$scan_dir/${scan_name}.meansBySample.pseudo.$chr.txt" ]] && means_pseudo_count=$((means_pseudo_count+1))
    done
    if [[ $scan_count -gt 0 ]]; then
        echo "    per-chr scan files (.scan.)       : $scan_count / 5"
    fi
    if [[ $pseudo_count -gt 0 ]]; then
        echo "    per-chr scan files (.pseudoscan.) : $pseudo_count / 5  [legacy naming]"
    fi
    if [[ $scan_count -eq 0 && $pseudo_count -eq 0 ]]; then
        echo "    per-chr scan files : 0 / 5  (deleted after concat, or not yet run)"
    fi
    echo "    per-chr means files: $means_count / 5"

    # Concatenated outputs — check both naming conventions
    for base_suffix in "scan.txt" "meansBySample.txt" "5panel.Mb.png" "5panel.cM.png" "Manhattan.png"; do
        f_new="$scan_dir/${scan_name}.$base_suffix"
        f_old="$scan_dir/${scan_name}.pseudoscan.${base_suffix#scan.}"  # only relevant for scan.txt

        if [[ -f "$f_new" ]]; then
            size=$(du -sh "$f_new" | cut -f1)
            echo "    [OK]  ${scan_name}.$base_suffix  ($size)"
        else
            # Check legacy name only for scan.txt
            f_legacy="$scan_dir/${scan_name}.pseudoscan.txt"
            if [[ "$base_suffix" == "scan.txt" && -f "$f_legacy" ]]; then
                size=$(du -sh "$f_legacy" | cut -f1)
                echo "    [OK]  ${scan_name}.pseudoscan.txt  ($size)  [legacy name — rename to .scan.txt to use with current scripts]"
            else
                echo "    [--]  ${scan_name}.$base_suffix"
            fi
        fi
    done
done

[[ $found_scans -eq 0 ]] && echo "  (none found)"

echo ""
echo "========================================"
