#!/bin/bash
#SBATCH --job-name=smooth_r2_diag
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1:00:00

# Run from your project directory:
#   sbatch pipeline/scripts/temp/smooth_r2_diag.sh \
#       --hapsdir  process/ZINC2 \
#       --smoothdir process/ZINC2/ZINC2_F_v3 \
#       --scan     ZINC2_F_v3 \
#       --rfile    helpfiles/ZINC2/design.txt \
#       --out      process/ZINC2/ZINC2_F_v3/smooth_r2_diag.txt

while [[ $# -gt 0 ]]; do
  case $1 in
    --hapsdir)   HAPSDIR="$2";   shift 2 ;;
    --smoothdir) SMOOTHDIR="$2"; shift 2 ;;
    --scan)      SCAN="$2";      shift 2 ;;
    --rfile)     RFILE="$2";     shift 2 ;;
    --out)       OUTFILE="$2";   shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -z "$OUTFILE" ]] && OUTFILE="${SMOOTHDIR}/smooth_r2_diag.txt"

module load R/4.2.2

Rscript "$(dirname $(readlink -f $0))/smooth_r2_diag.R" \
    --hapsdir   "${HAPSDIR}"   \
    --smoothdir "${SMOOTHDIR}" \
    --scan      "${SCAN}"      \
    --rfile     "${RFILE}"     \
    | tee "${OUTFILE}"

echo "Output written to ${OUTFILE}"
