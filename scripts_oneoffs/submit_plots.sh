#!/bin/bash
###############################################################################
# submit_plots.sh
#
# Submit plot jobs for all completed scans.
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
#
#   bash scripts/submit_plots.sh
#
# Downloads (after jobs finish):
#   scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/figures/*.png .
###############################################################################

set -e
XQTL2="/dfs7/adl/tdlong/fly_pool/XQTL2"
MAL_BASE="/dfs7/adl/tdlong/fly_pool/newpipeline_Nov23"
SARAH="/dfs7/adl/sruckman/XQTL/XQTL2"
FLYMAP="$XQTL2/helpfiles/flymap.r6.txt"
FIGDIR="$XQTL2/figures"
mkdir -p "$FIGDIR"

echo "=== Submitting plot jobs ==="
echo ""

# ── ZINC2 overlay (M + F, candidate gene annotation) ─────────────────────────
# Malathion gene positions (FlyBase r6):
#   Ace: chr3R 9069721 / Cyp6g1: chr2R 9064616 / Mdr65: chr2R 9163757
GENE_FILE="$XQTL2/helpfiles/malathion_genes.txt"
if [ ! -f "$GENE_FILE" ]; then
    printf "gene_name\tchr\tpos_bp\nAce\tchr3R\t9069721\nCyp6g1\tchr2R\t9064616\nMdr65\tchr2R\t9163757\n" \
        > "$GENE_FILE"
    echo "Created: $GENE_FILE"
fi

JOB_ZINC=$(sbatch --parsable \
    -A tdlong_lab -p standard \
    --wrap="module load R/4.2.2 && \
            Rscript $XQTL2/scripts/plot_zinc_overlay.R \
                $XQTL2/process/ZINC2/ZINC2_M \
                $XQTL2/process/ZINC2/ZINC2_F \
                $FIGDIR/zinc2_overlay \
                $FLYMAP \
                $GENE_FILE && \
            echo 'ZINC2 overlay done'")
echo "  ZINC2 overlay: job $JOB_ZINC  →  $FIGDIR/zinc2_overlay.png"

# ── Pupation TopMid scan ──────────────────────────────────────────────────────
JOB_PUP=$(sbatch --parsable \
    -A tdlong_lab -p standard \
    --wrap="module load R/4.2.2 && \
            Rscript $XQTL2/scripts/plot_single_scan.R \
                $SARAH/process/pupal_dynamic/TopMid \
                $FIGDIR/pupation_TopMid \
                $FLYMAP \
                'Pupation: Top vs Middle' && \
            echo 'Pupation plot done'")
echo "  Pupation TopMid: job $JOB_PUP  →  $FIGDIR/pupation_TopMid.png"

# ── Malathion scan ────────────────────────────────────────────────────────────
JOB_MAL=$(sbatch --parsable \
    -A tdlong_lab -p standard \
    --wrap="module load R/4.2.2 && \
            Rscript $XQTL2/scripts/plot_single_scan.R \
                $MAL_BASE/process \
                $FIGDIR/malathion \
                $FLYMAP \
                'Malathion resistance' \
                $GENE_FILE && \
            echo 'Malathion plot done'")
echo "  Malathion:       job $JOB_MAL  →  $FIGDIR/malathion.png"

echo ""
echo "Monitor: squeue -u \$USER"
echo ""
echo "When all jobs finish, download all PNGs:"
echo "  scp 'tdlong@hpc3.rcic.uci.edu:$FIGDIR/*.png' analysis/figures/"
echo ""
echo "NOTE: Aging scan failed — check log before submitting aging plot."
echo "  cat slurm-49687503_1.out   (or ls slurm-*.out | tail -5)"
