#!/bin/bash
# dump_configs.sh
# Prints the contents of all config/plot scripts referenced in run_all_scans_125.sh
# Run from the XQTL2 project root on the cluster.

FILES=(
    configs/zinc2_125.R
    configs/aging_125.R
    configs/pupation_125.R
    configs/malathion_125.R
    scripts_oneoffs/plot_pseudoscan.R
    scripts/plot_pseudoscan.R
)

for f in "${FILES[@]}"; do
    if [ -f "$f" ]; then
        echo "########################################"
        echo "# FILE: $f"
        echo "########################################"
        cat "$f"
        echo ""
    else
        echo "# NOT FOUND: $f"
        echo ""
    fi
done
