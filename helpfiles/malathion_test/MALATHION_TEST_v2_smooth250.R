#!/usr/bin/env Rscript
# MALATHION_TEST_v2_smooth250.R
#
# Publication figure for the malathion freqsmooth pipeline — 250 kb smoothing.
# Produces three scan figures + a meansBySample QC plot, then bundles all
# outputs into a single tar for download.
#
# Run from XQTL2 project root:
#   Rscript helpfiles/malathion_test/MALATHION_TEST_v2_smooth250.R

library(tidyverse)

SCAN_DIR  <- "process/malathion_test/MALATHION_TEST_v2_smooth250"
SCAN      <- "MALATHION_TEST_v2_smooth250"
SCAN_FILE <- file.path(SCAN_DIR, paste0(SCAN, ".scan.txt"))
OUT_TAR   <- file.path(SCAN_DIR, paste0(SCAN, ".hap.tar.gz"))

# ── Figure 1: Haplotype Wald -log10(p) ────────────────────────────────────────
SCAN_FILES   <- SCAN_FILE
SCAN_LABELS  <- NULL
SCAN_COLOURS <- c("#1F78B4")
THRESHOLD    <- 10
PEAKS        <- NULL
GENES        <- NULL
OUT_FILE     <- file.path(SCAN_DIR, paste0(SCAN, ".wald.png"))
FORMAT       <- "powerpoint"
source("scripts/plot_pseudoscan.R")

# ── Figure 2: Falconer H² ─────────────────────────────────────────────────────
SCAN_FILES   <- SCAN_FILE
SCAN_LABELS  <- NULL
SCAN_COLOURS <- c("#E31A1C")
YVAR         <- "Falc_H2"
YLAB         <- "Falconer H\u00b2"
THRESHOLD    <- NULL
OUT_FILE     <- file.path(SCAN_DIR, paste0(SCAN, ".falcH2.png"))
FORMAT       <- "powerpoint"
source("scripts/plot_freqsmooth_H2.R")

# ── Figure 3: Cutler H² ───────────────────────────────────────────────────────
SCAN_FILES   <- SCAN_FILE
SCAN_LABELS  <- NULL
SCAN_COLOURS <- c("#33A02C")
YVAR         <- "Cutl_H2"
YLAB         <- "Cutler H\u00b2"
THRESHOLD    <- NULL
OUT_FILE     <- file.path(SCAN_DIR, paste0(SCAN, ".cutlH2.png"))
FORMAT       <- "powerpoint"
source("scripts/plot_freqsmooth_H2.R")

# ── Figure 4: meansBySample QC — all founders, chr3L ─────────────────────────
means_file <- file.path(SCAN_DIR, paste0(SCAN, ".meansBySample.txt"))
qc_file    <- file.path(SCAN_DIR, paste0(SCAN, ".means_qc_chr3L.png"))

means_df <- read.table(means_file, header = TRUE) %>%
  filter(chr == "chr3L") %>%
  mutate(
    pos_mb = pos / 1e6,
    Pool   = ifelse(TRT == "C", "Control", "Selected"),
    Rep    = factor(REP)
  )

p_qc <- ggplot(means_df, aes(x = pos_mb, y = freq, colour = Rep, linetype = Pool)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  facet_wrap(~ founder, ncol = 4) +
  scale_linetype_manual(values = c(Control = "dashed", Selected = "solid")) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    title    = paste(SCAN, "— smoothed founder frequencies, chr3L"),
    x        = "Position (Mb)",
    y        = "Founder frequency",
    colour   = "Replicate",
    linetype = "Treatment"
  ) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", strip.background = element_blank())

ggsave(qc_file, p_qc, width = 10, height = 6, dpi = 150)
cat("Saved:", qc_file, "\n")

# ── Bundle for download ────────────────────────────────────────────────────────
cat("Writing", OUT_TAR, "\n")
system(paste(
  "tar -czf", OUT_TAR, "-C", SCAN_DIR,
  paste0(SCAN, ".wald.png"),
  paste0(SCAN, ".falcH2.png"),
  paste0(SCAN, ".cutlH2.png"),
  paste0(SCAN, ".means_qc_chr3L.png"),
  paste0(SCAN, ".scan.txt"),
  paste0(SCAN, ".meansBySample.txt")
))
cat("Done.\n")
