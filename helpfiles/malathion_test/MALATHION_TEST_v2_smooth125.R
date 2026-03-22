#!/usr/bin/env Rscript
# MALATHION_TEST_v2_smooth125.R
#
# Publication figure for the malathion freqsmooth pipeline test.
# Produces three 5-panel (per-chromosome) figures using the standard
# plot_pseudoscan.R and plot_freqsmooth_H2.R engines, then bundles
# all outputs into a single tar for download.
#
# Run from XQTL2 project root:
#   Rscript helpfiles/malathion_test/MALATHION_TEST_v2_smooth125.R

SCAN_DIR  <- "process/malathion_test/MALATHION_TEST_v2_smooth125"
SCAN      <- "MALATHION_TEST_v2_smooth125"
SCAN_FILE <- file.path(SCAN_DIR, paste0(SCAN, ".scan.txt"))
OUT_TAR   <- file.path(SCAN_DIR, paste0(SCAN, ".hap.tar.gz"))

# ── Figure 1: Haplotype Wald -log10(p) ────────────────────────────────────────
SCAN_FILES   <- SCAN_FILE
SCAN_LABELS  <- NULL
SCAN_COLOURS <- c("#1F78B4")
THRESHOLD    <- 10
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

# ── Bundle for download ────────────────────────────────────────────────────────
means_txt <- file.path(SCAN_DIR, paste0(SCAN, ".meansBySample.txt"))
cat("Writing", OUT_TAR, "\n")
system(paste(
  "tar -czf", OUT_TAR, "-C", SCAN_DIR,
  paste0(SCAN, ".wald.png"),
  paste0(SCAN, ".falcH2.png"),
  paste0(SCAN, ".cutlH2.png"),
  paste0(SCAN, ".scan.txt"),
  if (file.exists(means_txt)) paste0(SCAN, ".meansBySample.txt") else ""
))
cat("Done.\n")
