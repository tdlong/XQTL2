#!/bin/bash
# fix_move_to_oneoffs.sh
# One-time fix: the prior script incorrectly moved files to analysis/.
# This moves them to scripts_oneoffs/ where they belong.
# Safe: only moves files that exist, does not delete anything.

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "Working in: $ROOT"
mkdir -p scripts_oneoffs

FILES=(
    aging.R
    malathion.R
    zinc2.R
    pupation.R
    run_all_scans.sh
    run_all_scans_125.sh
    run_pupation_scan.sh
    run_malathion_scan.sh
    run_aging_scan.sh
    submit_freqs100_scans.sh
    submit_slide_plots.sh
    submit_plots.sh
    plot_pseudoscan.R
    plot_single_scan.R
    plot_zinc_overlay.R
    plot_for_slides.R
    make_scan_plot.R
    diagnose_scans.sh
    cluster_check.sh
    check_param_files.sh
    extract_downsampling.R
    diag_freqsmooth.R
    find_refalt_files.sh
    fix_refalt_b3852.sh
    merge_b3852.sh
    fq2bam_illumina.sh
    malathion_genes.txt
    diagnose_scans.out
    refalt_report.txt
    param_report.txt
    cluster_report.txt
    run_freqs250_all.sh
    run_chr2R_freqsmooth.sh
    run_stable_all.sh
)

for f in "${FILES[@]}"; do
    src="analysis/$f"
    if [ -f "$src" ]; then
        echo "  mv $src -> scripts_oneoffs/$f"
        mv "$src" "scripts_oneoffs/$f"
    fi
done

echo ""
echo "Done. Remaining contents of analysis/:"
ls analysis/ 2>/dev/null || echo "  (empty)"
