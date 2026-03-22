#!/bin/bash
# run_snp_scan.sh — submit the SNP scan pipeline with one command
#
# Imputes per-SNP allele frequencies from smoothed haplotype estimates
# and runs a Wald test at every SNP. Requires that smooth_haps has
# already been run (via run_scan.sh or standalone).
#
# Usage:
#   bash scripts/run_snp_scan.sh \
#       --design    helpfiles/<project>/design.txt \
#       --dir       process/<project> \
#       --scan      <scan_name> \
#       --snp-table helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
#       --founders  A1,A2,A3,A4,A5,A6,A7,AB8 \
#       --after     <jobid>

set -e

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case $1 in
    --design)    DESIGN="$2";    shift 2 ;;
    --dir)       DIR="$2";       shift 2 ;;
    --scan)      SCAN="$2";      shift 2 ;;
    --snp-table) SNP_TABLE="$2"; shift 2 ;;
    --founders)  FOUNDERS="$2";  shift 2 ;;
    --after)     AFTER="$2";     shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# ── Validate required arguments ──────────────────────────────────────────────
missing=""
[[ -z "$DESIGN" ]]    && missing="$missing --design"
[[ -z "$DIR" ]]       && missing="$missing --dir"
[[ -z "$SCAN" ]]      && missing="$missing --scan"
[[ -z "$SNP_TABLE" ]] && missing="$missing --snp-table"
[[ -z "$FOUNDERS" ]]  && missing="$missing --founders"
if [[ -n "$missing" ]]; then
    echo "Error: missing required arguments:$missing" >&2
    echo "Usage: bash scripts/run_snp_scan.sh --design <file> --dir <dir> --scan <name> --snp-table <file> --founders <list> [options]" >&2
    exit 1
fi

OUTDIR=${DIR}/${SCAN}

# ── Build dependency ─────────────────────────────────────────────────────────
DEP=""
[[ -n "$AFTER" ]] && DEP="--dependency=afterok:${AFTER}"

# ── SNP scan ─────────────────────────────────────────────────────────────────
jid_snp=$(sbatch --parsable ${DEP} \
    --array=1-5 scripts/snp_scan.sh \
    --rfile     "${DESIGN}" \
    --dir       "${DIR}" \
    --outdir    "${SCAN}" \
    --snp-table "${SNP_TABLE}" \
    --founders  "${FOUNDERS}")
echo "snp_scan:   $jid_snp"

# ── concat SNP scan chromosomes ──────────────────────────────────────────────
jid_concat=$(sbatch --parsable --dependency=afterok:${jid_snp} \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=0:10:00 \
    --wrap="bash scripts/concat_scans.sh --snp ${OUTDIR}")
echo "snp_concat: $jid_concat"

echo "done:       ${jid_concat}"
