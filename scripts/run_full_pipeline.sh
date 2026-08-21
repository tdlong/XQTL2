#!/bin/bash
###############################################################################
# run_full_pipeline.sh — end-to-end XQTL pipeline with SLURM dependencies
#
# Chains Steps 2–6: fq2bam → call_samples → REFALT2haps → scan → figures.
# Each step waits for the previous one to finish successfully before starting.
#
# REFALT comes from the founder-catalog caller (call_samples.sh). The catalog
# itself is NOT built here: build_catalog.sh overwrites its --out with no
# staleness check, so building is a deliberate one-time act you run yourself.
# Build it once, then point every run at it with --catalog:
#
#   bash pipeline/scripts/build_catalog.sh \
#       --founders pipeline/helpfiles/A_founders.bams.txt \
#       --out      process/<project>/Catalog
#
# Usage:
#   bash pipeline/scripts/run_full_pipeline.sh \
#       --project   <project> \
#       --barcodes  helpfiles/<project>/<project>.barcodes.txt \
#       --rawdir    data/raw/<project> \
#       --bamdir    data/bam/<project> \
#       --catalog   process/<project>/Catalog \
#       --parfile   helpfiles/<project>/hap_params.R \
#       --design    helpfiles/<project>/design.txt \
#       --scan      <scan_name> \
#       --founders  A \
#       --snp-scan
#
# --catalog defaults to process/<project>/Catalog.
#
# --snp-scan adds the SNP scan (Step 5b). It needs no extra inputs: founder
# states are read from RefAlt, and founder names from --parfile.
#
# Skip early steps if you already have BAMs or haplotypes:
#       --skip-fq2bam         start at Step 3 (REFALT)
#       --skip-refalt         skip REFALT creation; re-run haplotypes + scan
#       --skip-haps           start at Step 5 (scan) — haplotypes must exist
#       --after <jid>         make hap scan wait for this SLURM job (use with --skip-haps
#                             to share one haplotype job across multiple scans)
#
# All SLURM flags (--mem-per-cpu, -p, -A, etc.) are passed through.
###############################################################################
set -e

# ── Defaults ────────────────────────────────────────────────────────────────
SMOOTH_KB=250
MEM_PER_CPU=3G
CPUS_PER_TASK=1
PARTITION=standard
ACCOUNT=tdlong_lab
SKIP_FQ2BAM=false
SKIP_REFALT=false
SKIP_HAPS=false
SNP_SCAN=false
AFTER_HAPS=""

# ── Parse arguments ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case $1 in
    --project)       PROJECT="$2";       shift 2 ;;
    --barcodes)      BARCODES="$2";      shift 2 ;;
    --rawdir)        RAWDIR="$2";        shift 2 ;;
    --bamdir)        BAMDIR="$2";        shift 2 ;;
    --parfile)       PARFILE="$2";       shift 2 ;;
    --design)        DESIGN="$2";        shift 2 ;;
    --catalog)       CATALOG="$2";       shift 2 ;;
    --scan)          SCAN="$2";          shift 2 ;;
    --smooth)        SMOOTH_KB="$2";     shift 2 ;;
    --founders)      FOUNDERS="$2";      shift 2 ;;  # A or B
    --snp-scan)      SNP_SCAN=true;      shift ;;
    --mem-per-cpu)   MEM_PER_CPU="$2";   shift 2 ;;
    --cpus-per-task) CPUS_PER_TASK="$2"; shift 2 ;;
    -p|--partition)  PARTITION="$2";     shift 2 ;;
    -A|--account)    ACCOUNT="$2";       shift 2 ;;
    --skip-fq2bam)   SKIP_FQ2BAM=true;  shift ;;
    --skip-refalt)   SKIP_REFALT=true;   shift ;;
    --skip-haps)     SKIP_HAPS=true;     shift ;;
    --after)         AFTER_HAPS="$2";    shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# ── Validate ────────────────────────────────────────────────────────────────
missing=""
[[ -z "$PROJECT" ]] && missing="$missing --project"
[[ -z "$DESIGN" ]]  && missing="$missing --design"
[[ -z "$PARFILE" ]] && missing="$missing --parfile"
[[ -z "$SCAN" ]]    && missing="$missing --scan"
if [[ -n "$missing" ]]; then
    echo "Error: missing required arguments:$missing" >&2
    exit 1
fi

PROCESSDIR="process/${PROJECT}"
mkdir -p "${PROCESSDIR}"

SLURM_COMMON="-A ${ACCOUNT}"

# Needed by both the REFALT step and the SNP-table step below, so default it here
# rather than inside the REFALT branch (which --skip-refalt skips).
[[ -z "${CATALOG:-}" ]] && CATALOG="${PROCESSDIR}/Catalog"

echo "=== XQTL full pipeline: ${PROJECT} / ${SCAN} ==="
echo ""

# ── Step 2: fq2bam ──────────────────────────────────────────────────────────
AFTER_BAM=""
if [[ "$SKIP_FQ2BAM" == false ]]; then
    [[ -z "$BARCODES" ]] && { echo "Error: --barcodes required unless --skip-fq2bam" >&2; exit 1; }
    [[ -z "$RAWDIR" ]]   && { echo "Error: --rawdir required unless --skip-fq2bam" >&2; exit 1; }
    [[ -z "$BAMDIR" ]]   && BAMDIR="data/bam/${PROJECT}"
    mkdir -p "${BAMDIR}"

    NN=$(wc -l < "${BARCODES}")
    jid_bam=$(sbatch --parsable ${SLURM_COMMON} \
        --array=1-${NN} "pipeline/scripts/fq2bam.sh" \
        "${BARCODES}" "${RAWDIR}" "${BAMDIR}")
    echo "fq2bam:     ${jid_bam}  (${NN} samples)"
    AFTER_BAM="--dependency=afterok:${jid_bam}"
fi

# ── Step 3: REFALT — count samples against the founder catalog ──────────────
AFTER_REFALT=""
if [[ "$SKIP_REFALT" == false ]]; then
    # Build bam_list if it doesn't exist yet
    BAMLIST="helpfiles/${PROJECT}/bam_list.txt"
    if [[ ! -f "$BAMLIST" ]]; then
        [[ -z "$BAMDIR" ]]    && BAMDIR="data/bam/${PROJECT}"
        [[ -z "$FOUNDERS" ]]  && { echo "Error: --founders (A or B) required to build bam_list" >&2; exit 1; }
        mkdir -p "helpfiles/${PROJECT}"
        ls "${BAMDIR}"/*.bam > "${BAMLIST}"
        cat "pipeline/helpfiles/${FOUNDERS}_founders.bams.txt" >> "${BAMLIST}"
        echo "Built ${BAMLIST} — review before results are final"
    fi

    # The catalog itself is never built here — see the header. call_samples.sh
    # errors out if it is missing.

    # call_samples.sh takes a bare job ID, not an sbatch --dependency string.
    AFTER_BAM_JID=""
    [[ -n "$AFTER_BAM" ]] && AFTER_BAM_JID=$(echo "$AFTER_BAM" | grep -o '[0-9]*')

    # It prints progress then the merge job ID last.
    jid_refalt=$(bash pipeline/scripts/call_samples.sh \
        --catalog "${CATALOG}" \
        --bamlist "${BAMLIST}" \
        --dir     "${PROCESSDIR}" \
        ${AFTER_BAM_JID:+--after ${AFTER_BAM_JID}} \
        -p "${PARTITION}" -A "${ACCOUNT}" | tail -1)
    echo "REFALT:     ${jid_refalt}  (catalog: ${CATALOG})"
    AFTER_REFALT="--dependency=afterok:${jid_refalt}"
fi

# ── Step 4: REFALT2haps ─────────────────────────────────────────────────────
if [[ "$SKIP_HAPS" == false ]]; then
    AFTER_HAPS_FLAG=""
    [[ -n "$AFTER_REFALT" ]] && AFTER_HAPS_FLAG="--after $(echo $AFTER_REFALT | grep -o '[0-9]*')"
    # --skip-refalt with an externally submitted REFALT job (e.g. call_samples.sh run
    # by hand to count only newly added samples): chain haplotypes on it via --after.
    [[ -z "$AFTER_HAPS_FLAG" && -n "$AFTER_HAPS" ]] && AFTER_HAPS_FLAG="--after ${AFTER_HAPS}"
    jid_haps=$(bash pipeline/scripts/run_haps.sh \
        --parfile "${PARFILE}" --dir "${PROCESSDIR}" \
        ${AFTER_HAPS_FLAG} \
        -A "${ACCOUNT}")
    echo "haplotypes: ${jid_haps}"
else
    jid_haps="${AFTER_HAPS}"
    echo "haplotypes: skipped (using ${jid_haps})"
fi

# ── Step 5a: haplotype scan ─────────────────────────────────────────────────
scan_out=$(bash pipeline/scripts/run_scan.sh \
    --design "${DESIGN}" --dir "${PROCESSDIR}" --scan "${SCAN}" \
    --smooth "${SMOOTH_KB}" --mem-per-cpu "${MEM_PER_CPU}" \
    --cpus-per-task "${CPUS_PER_TASK}" -p "${PARTITION}" -A "${ACCOUNT}" \
    --after "${jid_haps}")
jid_hap=$(echo "$scan_out" | grep "^done:" | awk '{print $2}')
echo "$scan_out" | sed 's/^/  /'

# ── Step 5b: SNP scan (--snp-scan) ──────────────────────────────────────────
# Founder states come from RefAlt, which the haplotypes were also fit from, so
# there is nothing extra to build or pass (XQTL2 #35).
jid_snp=""
if [[ "$SNP_SCAN" == true ]]; then
    snp_out=$(bash pipeline/scripts/run_snp_scan.sh \
        --design "${DESIGN}" --dir "${PROCESSDIR}" --scan "${SCAN}" \
        --parfile "${PARFILE}" \
        --mem-per-cpu "${MEM_PER_CPU}" --cpus-per-task "${CPUS_PER_TASK}" \
        -p "${PARTITION}" -A "${ACCOUNT}" \
        --after "${jid_haps}")
    jid_snp=$(echo "$snp_out" | grep "^done:" | awk '{print $2}')
    echo "$snp_out" | sed 's/^/  /'
fi

# ── Step 6: figures + tarball ───────────────────────────────────────────────
SCAN_DIR="${PROCESSDIR}/Scans/${SCAN}"   # stage layout (see reorganize_project.sh)
FIG_DEPS="--dependency=afterok:${jid_hap}"
[[ -n "$jid_snp" ]] && FIG_DEPS="--dependency=afterok:${jid_hap}:${jid_snp}"

# Build the figure commands
FIG_CMD="module load R/4.2.2"
FIG_CMD="${FIG_CMD} && Rscript pipeline/scripts/plot_5panel.R --scan ${SCAN_DIR}/${SCAN}.scan.txt --out ${SCAN_DIR}/${SCAN}.wald.png --format powerpoint --threshold 10"
FIG_CMD="${FIG_CMD} && Rscript pipeline/scripts/plot_manhattan.R --scan ${SCAN_DIR}/${SCAN}.scan.txt --out ${SCAN_DIR}/${SCAN}.manhattan.png --format powerpoint --threshold 10"
FIG_CMD="${FIG_CMD} && Rscript pipeline/scripts/plot_H2_overlay.R --scan ${SCAN_DIR}/${SCAN}.scan.txt --out ${SCAN_DIR}/${SCAN}.H2.png --format powerpoint"

if [[ "$SNP_SCAN" == true ]]; then
    FIG_CMD="${FIG_CMD} && Rscript pipeline/scripts/plot_freqsmooth_snp.R --scan ${SCAN_DIR}/${SCAN}.snp_scan.txt --out ${SCAN_DIR}/${SCAN}.snp.wald.png --format powerpoint --threshold 10"
fi

FIG_CMD="${FIG_CMD} && cd ${SCAN_DIR} && tar -czf ${SCAN}.tar.gz *.txt *.png"

jid_fig=$(sbatch --parsable ${FIG_DEPS} \
    ${SLURM_COMMON} --time=1:00:00 \
    --wrap="${FIG_CMD}")
echo "figures:    ${jid_fig}"

echo ""
echo "=== All jobs submitted. Final job: ${jid_fig} ==="
echo "=== Results will be in: ${SCAN_DIR}/${SCAN}.tar.gz ==="
