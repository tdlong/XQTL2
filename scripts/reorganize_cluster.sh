#!/bin/bash
# reorganize_cluster.sh
# Run ONCE from the XQTL2 project root to move experiment-specific files
# into scripts_oneoffs/ and clean up the root directory.
# Safe: only moves files, does not delete anything.

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "Working in: $ROOT"
mkdir -p scripts_oneoffs

# -----------------------------------------------------------------------
# Move experiment-specific scripts out of scripts/ -> scripts_oneoffs/
# (freqsmooth and stable are kept in scripts/ — they will be tracked)
# -----------------------------------------------------------------------
ONEOFF_SCRIPTS=(
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
)

for f in "${ONEOFF_SCRIPTS[@]}"; do
    src="scripts/$f"
    if [ -f "$src" ]; then
        echo "  mv $src -> scripts_oneoffs/$f"
        mv "$src" "scripts_oneoffs/$f"
    fi
done

# -----------------------------------------------------------------------
# Move root-level clutter -> scripts_oneoffs/
# -----------------------------------------------------------------------
ROOT_CLUTTER=(
    diagnose_scans.out
    refalt_report.txt
    param_report.txt
    cluster_report.txt
    run_freqs250_all.sh
    run_chr2R_freqsmooth.sh
    run_stable_all.sh
)

for f in "${ROOT_CLUTTER[@]}"; do
    if [ -f "$f" ]; then
        echo "  mv $f -> scripts_oneoffs/$f"
        mv "$f" "scripts_oneoffs/$f"
    fi
done

# -----------------------------------------------------------------------
# Large files: just report their presence (do not move — user decides)
# -----------------------------------------------------------------------
echo ""
echo "Large files in root (not moved — handle manually):"
for f in founders_bam_files.tar AGE_SY.tar .RData .Rhistory; do
    if [ -f "$f" ]; then
        ls -lh "$f"
    fi
done

echo ""
echo "Done. Review scripts_oneoffs/ and then run: git status"
