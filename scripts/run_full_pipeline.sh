#!/bin/bash
###############################################################################
# run_full_pipeline.sh — end-to-end XQTL pipeline with SLURM dependencies
#
# Chains Steps 2–6: fq2bam → bam2bcf2REFALT → REFALT2haps → scan → figures.
# Each step waits for the previous one to finish successfully before starting.
#
# Usage:
#   bash pipeline/scripts/run_full_pipeline.sh \
#       --project   <project> \
#       --barcodes  helpfiles/<project>/<project>.barcodes.txt \
#       --rawdir    data/raw/<project> \
#       --bamdir    data/bam/<project> \
#       --parfile   helpfiles/<project>/hap_params.R \
#       --design    helpfiles/<project>/design.txt \
#       --scan      <scan_name> \
#       --founders  A \
#       --snp-table pipeline/helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
#       --founder-list A1,A2,A3,A4,A5,A6,A7,AB8
#
# Skip early steps if you already have BAMs or haplotypes:
#       --skip-fq2bam         start at Step 3 (REFALT)
#       --skip-refalt         start at Step 5 (scan) — haplotypes must exist
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

# ── Parse arguments ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case $1 in
    --project)       PROJECT="$2";       shift 2 ;;
    --barcodes)      BARCODES="$2";      shift 2 ;;
    --rawdir)        RAWDIR="$2";        shift 2 ;;
    --bamdir)        BAMDIR="$2";        shift 2 ;;
    --parfile)       PARFILE="$2";       shift 2 ;;
    --design)        DESIGN="$2";        shift 2 ;;
    --scan)          SCAN="$2";          shift 2 ;;
    --smooth)        SMOOTH_KB="$2";     shift 2 ;;
    --founders)      FOUNDERS="$2";      shift 2 ;;  # A or B
    --snp-table)     SNP_TABLE="$2";     shift 2 ;;
    --founder-list)  FOUNDER_LIST="$2";  shift 2 ;;
    --mem-per-cpu)   MEM_PER_CPU="$2";   shift 2 ;;
    --cpus-per-task) CPUS_PER_TASK="$2"; shift 2 ;;
    -p|--partition)  PARTITION="$2";     shift 2 ;;
    -A|--account)    ACCOUNT="$2";       shift 2 ;;
    --skip-fq2bam)   SKIP_FQ2BAM=true;  shift ;;
    --skip-refalt)   SKIP_REFALT=true;   shift ;;
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

PIPELINE_DIR="$(dirname "$(readlink -f "$0")")"
PROCESSDIR="process/${PROJECT}"
mkdir -p "${PROCESSDIR}"

SLURM_COMMON="-A ${ACCOUNT} -p ${PARTITION} --cpus-per-task=${CPUS_PER_TASK} --mem-per-cpu=${MEM_PER_CPU}"

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
        --array=1-${NN} "${PIPELINE_DIR}/fq2bam.sh" \
        "${BARCODES}" "${RAWDIR}" "${BAMDIR}")
    echo "fq2bam:     ${jid_bam}  (${NN} samples)"
    AFTER_BAM="--dependency=afterok:${jid_bam}"
fi

# ── Step 3: bam2bcf2REFALT ──────────────────────────────────────────────────
AFTER_REFALT=""
if [[ "$SKIP_REFALT" == false ]]; then
    # Build bam_list if it doesn't exist yet
    BAMLIST="helpfiles/${PROJECT}/bam_list.txt"
    if [[ ! -f "$BAMLIST" ]]; then
        [[ -z "$BAMDIR" ]]    && BAMDIR="data/bam/${PROJECT}"
        [[ -z "$FOUNDERS" ]]  && { echo "Error: --founders (A or B) required to build bam_list" >&2; exit 1; }
        mkdir -p "helpfiles/${PROJECT}"
        ls "${BAMDIR}"/*.bam > "${BAMLIST}"
        cat "${PIPELINE_DIR}/../helpfiles/${FOUNDERS}_founders.bams.txt" >> "${BAMLIST}"
        echo "Built ${BAMLIST} — review before results are final"
    fi

    jid_refalt=$(sbatch --parsable ${AFTER_BAM} ${SLURM_COMMON} \
        --array=1-5 "${PIPELINE_DIR}/bam2bcf2REFALT.sh" \
        "${BAMLIST}" "${PROCESSDIR}")
    echo "REFALT:     ${jid_refalt}"
    AFTER_REFALT="--dependency=afterok:${jid_refalt}"
fi

# ── Step 4: REFALT2haps ─────────────────────────────────────────────────────
jid_haps=$(sbatch --parsable ${AFTER_REFALT} ${SLURM_COMMON} \
    --array=1-5 "${PIPELINE_DIR}/REFALT2haps.sh" \
    --parfile "${PARFILE}" --dir "${PROCESSDIR}")
echo "haplotypes: ${jid_haps}"

# ── Step 5a: haplotype scan ─────────────────────────────────────────────────
scan_out=$(bash "${PIPELINE_DIR}/run_scan.sh" \
    --design "${DESIGN}" --dir "${PROCESSDIR}" --scan "${SCAN}" \
    --smooth "${SMOOTH_KB}" --mem-per-cpu "${MEM_PER_CPU}" \
    --cpus-per-task "${CPUS_PER_TASK}" -p "${PARTITION}" -A "${ACCOUNT}" \
    --after "${jid_haps}")
jid_hap=$(echo "$scan_out" | grep "^done:" | awk '{print $2}')
echo "$scan_out" | sed 's/^/  /'

# ── Step 5b: SNP scan (if snp-table provided) ──────────────────────────────
jid_snp=""
if [[ -n "$SNP_TABLE" && -n "$FOUNDER_LIST" ]]; then
    snp_out=$(bash "${PIPELINE_DIR}/run_snp_scan.sh" \
        --design "${DESIGN}" --dir "${PROCESSDIR}" --scan "${SCAN}" \
        --snp-table "${SNP_TABLE}" --founders "${FOUNDER_LIST}" \
        --mem-per-cpu "${MEM_PER_CPU}" --cpus-per-task "${CPUS_PER_TASK}" \
        -p "${PARTITION}" -A "${ACCOUNT}" \
        --after "${jid_haps}")
    jid_snp=$(echo "$snp_out" | grep "^done:" | awk '{print $2}')
    echo "$snp_out" | sed 's/^/  /'
fi

# ── Step 6: figures + tarball ───────────────────────────────────────────────
SCAN_DIR="${PROCESSDIR}/${SCAN}"
FIG_DEPS="--dependency=afterok:${jid_hap}"
[[ -n "$jid_snp" ]] && FIG_DEPS="--dependency=afterok:${jid_hap}:${jid_snp}"

# Build the figure commands
FIG_CMD="module load R/4.2.2"
FIG_CMD="${FIG_CMD} && Rscript ${PIPELINE_DIR}/plot_5panel.R --scan ${SCAN_DIR}/${SCAN}.scan.txt --out ${SCAN_DIR}/${SCAN}.wald.png --format powerpoint --threshold 10"
FIG_CMD="${FIG_CMD} && Rscript ${PIPELINE_DIR}/plot_manhattan.R --scan ${SCAN_DIR}/${SCAN}.scan.txt --out ${SCAN_DIR}/${SCAN}.manhattan.png --format powerpoint --threshold 10"
FIG_CMD="${FIG_CMD} && Rscript ${PIPELINE_DIR}/plot_H2_overlay.R --scan ${SCAN_DIR}/${SCAN}.scan.txt --out ${SCAN_DIR}/${SCAN}.H2.png --format powerpoint"

if [[ -n "$SNP_TABLE" ]]; then
    FIG_CMD="${FIG_CMD} && Rscript ${PIPELINE_DIR}/plot_freqsmooth_snp.R --scan ${SCAN_DIR}/${SCAN}.snp_scan.txt --out ${SCAN_DIR}/${SCAN}.snp.wald.png --format powerpoint --threshold 10"
fi

FIG_CMD="${FIG_CMD} && cd ${SCAN_DIR} && tar -czf ${SCAN}.tar.gz *.txt *.png"

jid_fig=$(sbatch --parsable ${FIG_DEPS} \
    ${SLURM_COMMON} --time=1:00:00 \
    --wrap="${FIG_CMD}")
echo "figures:    ${jid_fig}"

echo ""
echo "=== All jobs submitted. Final job: ${jid_fig} ==="
echo "=== Results will be in: ${SCAN_DIR}/${SCAN}.tar.gz ==="
