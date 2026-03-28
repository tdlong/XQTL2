#!/bin/bash
###############################################################################
# debug_chrX_v2.sh — Rerun smooth for chrX with instrumented smooth_haps,
#                     then hap_scan, then frequency check plot.
#
# One command:
#   cd /dfs7/adl/tdlong/fly_pool/XQTL2 && bash scripts_oneoffs/ZINC2/debug_chrX_v2.sh
#
# All diagnostic output goes to scripts_oneoffs/ZINC2/debug_chrX_results/
# (git-trackable, NOT in process/). Results auto-pushed to XQTL2-dev.
###############################################################################

set -e
cd /dfs7/adl/tdlong/fly_pool/XQTL2

# Pull latest code from both repos
git pull dev main --rebase 2>/dev/null || true
git pull origin main 2>/dev/null || true

DESIGN=helpfiles/ZINC2/Zinc2.test.F.txt
DIR=process/ZINC2
SCAN=ZINC2_F_v3
RESDIR=scripts_oneoffs/ZINC2/debug_chrX_results
mkdir -p ${RESDIR}

# ── Job 1: Instrumented smooth (chrX only) ────────────────────────────────────
#     Writes diag to RESDIR (git-trackable), pipeline output to process/ as usual
JID_SMOOTH=$(sbatch --parsable \
    -A tdlong_lab -p standard \
    --cpus-per-task=2 --mem-per-cpu=6G --time=4:00:00 \
    --job-name=smooth_dbg_chrX \
    --wrap="module load R/4.2.2 && \
Rscript scripts_oneoffs/ZINC2/smooth_haps_debug.R \
    --chr chrX --dir ${DIR} --outdir ${SCAN} \
    --rfile ${DESIGN} --smooth-kb 250 \
    --diagdir ${RESDIR}")
echo "smooth: ${JID_SMOOTH}"

# ── Job 2: Haplotype scan (depends on smooth) ────────────────────────────────
JID_SCAN=$(sbatch --parsable \
    --dependency=afterok:${JID_SMOOTH} \
    -A tdlong_lab -p standard \
    --cpus-per-task=1 --mem-per-cpu=3G --time=2:00:00 \
    --job-name=hscan_dbg_chrX \
    --wrap="module load R/4.2.2 && \
Rscript scripts/hap_scan.R \
    --chr chrX --dir ${DIR}/${SCAN} --outdir ${SCAN} \
    --rfile ${DESIGN}")
echo "hap_scan: ${JID_SCAN}"

# ── Job 3: Frequency plot + push all results ──────────────────────────────────
JID_DIAG=$(sbatch --parsable \
    --dependency=afterany:${JID_SMOOTH},afterany:${JID_SCAN} \
    -A tdlong_lab -p standard \
    --cpus-per-task=1 --mem-per-cpu=3G --time=0:30:00 \
    --job-name=diag_dbg_chrX \
    --wrap="cd /dfs7/adl/tdlong/fly_pool/XQTL2 && \
module load R/4.2.2 && \
Rscript -e '
suppressPackageStartupMessages(library(tidyverse))

MEANS  <- \"${DIR}/${SCAN}/${SCAN}.meansBySample.chrX.txt\"
OUTPNG <- \"${RESDIR}/chrX_freq.png\"

df <- read.table(MEANS, header=TRUE) %>%
  as_tibble() %>%
  mutate(pos_mb = pos / 1e6,
         TRT    = factor(TRT, levels=c(\"C\",\"Z\")))

cat(sprintf(\"meansBySample: %d rows, %d neg, min=%.6f\\n\",
    nrow(df), sum(df\$freq < 0, na.rm=TRUE), min(df\$freq, na.rm=TRUE)))

p <- ggplot(df, aes(x=pos_mb, y=freq, group=interaction(TRT, REP),
                    colour=TRT, alpha=TRT)) +
  geom_line(linewidth=0.3) +
  scale_colour_manual(values=c(C=\"grey60\", Z=\"#CC3333\"), name=\"Treatment\") +
  scale_alpha_manual(values=c(C=0.5, Z=0.7), guide=\"none\") +
  facet_wrap(~founder, ncol=2, scales=\"free_y\") +
  labs(x=\"Position (Mb)\", y=\"Smoothed frequency\",
       title=\"ZINC2_F_v3: smoothed founder frequencies on chrX\",
       subtitle=\"All replicates; C=grey, Z=red\") +
  theme_bw(base_size=10) +
  theme(panel.grid.minor=element_blank(),
        strip.text=element_text(face=\"bold\"))

ggsave(OUTPNG, p, width=10, height=10, dpi=150)
cat(\"Written:\", OUTPNG, \"\\n\")
' && \
git add ${RESDIR}/ && \
git commit -m 'debug chrX v2: diag + freq plot' && \
git pull dev main --rebase && \
git push dev HEAD:main || echo 'WARNING: git push failed'")
echo "diag: ${JID_DIAG}"

echo ""
echo "Jobs: smooth=${JID_SMOOTH} -> hap_scan=${JID_SCAN} -> diag=${JID_DIAG}"
echo "Monitor: squeue -u \$USER"
echo "Results will be in: ${RESDIR}/"
