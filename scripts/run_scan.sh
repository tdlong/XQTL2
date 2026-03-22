#!/bin/bash
# run_scan.sh — submit the full scan pipeline with one command
#
# Submits smooth_haps → (hap_scan || snp_scan) → concat → figure
# with proper SLURM dependency chaining.
#
# Usage:
#   bash scripts/run_scan.sh \
#       --design    helpfiles/<project>/design.txt \
#       --dir       process/<project> \
#       --scan      <scan_name> \
#       --smooth    250 \
#       --snp-table helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
#       --founders  A1,A2,A3,A4,A5,A6,A7,AB8 \
#       --figure    helpfiles/<project>/figure.R \
#       --after     <jobid>

set -e

# ── Defaults ──────────────────────────────────────────────────────────────────
SMOOTH_KB=250

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case $1 in
    --design)    DESIGN="$2";    shift 2 ;;
    --dir)       DIR="$2";       shift 2 ;;
    --scan)      SCAN="$2";      shift 2 ;;
    --smooth)    SMOOTH_KB="$2"; shift 2 ;;
    --snp-table) SNP_TABLE="$2"; shift 2 ;;
    --founders)  FOUNDERS="$2";  shift 2 ;;
    --figure)    FIGURE="$2";    shift 2 ;;
    --after)     AFTER="$2";     shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# ── Validate required arguments ──────────────────────────────────────────────
missing=""
[[ -z "$DESIGN" ]] && missing="$missing --design"
[[ -z "$DIR" ]]    && missing="$missing --dir"
[[ -z "$SCAN" ]]   && missing="$missing --scan"
if [[ -n "$missing" ]]; then
    echo "Error: missing required arguments:$missing" >&2
    echo "Usage: bash scripts/run_scan.sh --design <file> --dir <dir> --scan <name> [options]" >&2
    exit 1
fi

# SNP scan requires both --snp-table and --founders
if [[ -n "$SNP_TABLE" && -z "$FOUNDERS" ]] || [[ -z "$SNP_TABLE" && -n "$FOUNDERS" ]]; then
    echo "Error: --snp-table and --founders must be used together" >&2
    exit 1
fi

OUTDIR=${DIR}/${SCAN}
mkdir -p "${OUTDIR}"

# ── Build dependency for smooth step ─────────────────────────────────────────
DEP_SMOOTH=""
[[ -n "$AFTER" ]] && DEP_SMOOTH="--dependency=afterok:${AFTER}"

# ── 5a: smooth haplotype frequencies ─────────────────────────────────────────
jid_smooth=$(sbatch --parsable ${DEP_SMOOTH} \
    --array=1-5 scripts/smooth_haps.sh \
    --rfile     "${DESIGN}" \
    --dir       "${DIR}" \
    --outdir    "${SCAN}" \
    --smooth-kb "${SMOOTH_KB}")
echo "smooth:     $jid_smooth"

# ── 5b: haplotype scan (Wald + H2) ──────────────────────────────────────────
jid_hap=$(sbatch --parsable --dependency=afterok:${jid_smooth} \
    --array=1-5 scripts/hap_scan.sh \
    --rfile  "${DESIGN}" \
    --dir    "${DIR}" \
    --outdir "${SCAN}")
echo "hap_scan:   $jid_hap"

# ── 6a: concat haplotype scan ───────────────────────────────────────────────
jid_concat=$(sbatch --parsable --dependency=afterok:${jid_hap} \
    -A tdlong_lab -p standard --mem=10G --time=1:00:00 \
    --wrap="bash scripts/concat_scans.sh ${OUTDIR}")
echo "concat:     $jid_concat"

# Track all final jobs for figure dependency
jid_final="${jid_concat}"

# ── 5c + 6b: SNP scan (only if --snp-table provided) ────────────────────────
if [[ -n "$SNP_TABLE" ]]; then
    jid_snp=$(sbatch --parsable --dependency=afterok:${jid_smooth} \
        --array=1-5 scripts/snp_scan.sh \
        --rfile     "${DESIGN}" \
        --dir       "${DIR}" \
        --outdir    "${SCAN}" \
        --snp-table "${SNP_TABLE}" \
        --founders  "${FOUNDERS}")
    echo "snp_scan:   $jid_snp"

    jid_snp_concat=$(sbatch --parsable --dependency=afterok:${jid_snp} \
        -A tdlong_lab -p standard --mem=10G --time=1:00:00 \
        --wrap="bash scripts/concat_scans.sh --snp ${OUTDIR}")
    echo "snp_concat: $jid_snp_concat"

    jid_final="${jid_concat}:${jid_snp_concat}"
fi

# ── 7: figure (only if --figure provided) ────────────────────────────────────
if [[ -n "$FIGURE" ]]; then
    jid_fig=$(sbatch --parsable --dependency=afterok:${jid_final} \
        -A tdlong_lab -p standard --mem=10G --time=1:00:00 \
        --wrap="module load R/4.2.2 && Rscript ${FIGURE}")
    echo "figure:     $jid_fig"
    jid_final="${jid_fig}"
fi

echo "done:       ${jid_final}"
