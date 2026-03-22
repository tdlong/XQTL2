#!/usr/bin/env Rscript
# MALATHION_TEST_v2_smooth250_snp.R
#
# SNP scan figures for malathion — 250 kb smoothing.
# Produces Wald + combined H² figures, then bundles into tar.
#
# Run from XQTL2 project root:
#   Rscript helpfiles/malathion_test/MALATHION_TEST_v2_smooth250_snp.R

SCAN_DIR  <- "process/malathion_test/MALATHION_TEST_v2_smooth250"
SCAN      <- "MALATHION_TEST_v2_smooth250"
SCAN_FILE <- file.path(SCAN_DIR, paste0(SCAN, ".snp_scan.txt"))
OUT_TAR   <- file.path(SCAN_DIR, paste0(SCAN, ".snp.tar.gz"))

# ── Figure 1: SNP Wald -log10(p) ──────────────────────────────────────────────
SCAN_FILES   <- SCAN_FILE
SCAN_LABELS  <- NULL
SCAN_COLOURS <- c("#1F78B4")
YVAR         <- "Wald_log10p"
YLAB         <- "-log10(p) [SNP Wald]"
THRESHOLD    <- 10
OUT_FILE     <- file.path(SCAN_DIR, paste0(SCAN, ".snp.wald.png"))
FORMAT       <- "powerpoint"
source("scripts/plot_freqsmooth_snp.R")

# ── Bundle for download ────────────────────────────────────────────────────────
cat("Writing", OUT_TAR, "\n")
system(paste(
  "tar -czf", OUT_TAR, "-C", SCAN_DIR,
  paste0(SCAN, ".snp.wald.png"),
  paste0(SCAN, ".snp_scan.txt")
))
cat("Done.\n")
