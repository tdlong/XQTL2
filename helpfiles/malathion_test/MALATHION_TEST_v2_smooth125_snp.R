#!/usr/bin/env Rscript
# MALATHION_TEST_v2_smooth125_snp.R
#
# SNP scan figures for malathion — 125 kb smoothing.
#
# Run from XQTL2 project root:
#   Rscript helpfiles/malathion_test/MALATHION_TEST_v2_smooth125_snp.R

SCAN_DIR  <- "process/malathion_test/MALATHION_TEST_v2_smooth125"
SCAN      <- "MALATHION_TEST_v2_smooth125"
OUT_TAR   <- file.path(SCAN_DIR, paste0(SCAN, ".snp.tar.gz"))

# ── Figure 1: SNP Wald -log10(p) ────────────────────────────────────────────
system(paste(
  "Rscript scripts/plot_freqsmooth_snp.R",
  "--scan",   file.path(SCAN_DIR, paste0(SCAN, ".snp_scan.txt")),
  "--out",    file.path(SCAN_DIR, paste0(SCAN, ".snp.wald.png")),
  "--format", "powerpoint",
  "--threshold", "10"
))

# ── Bundle for download ─────────────────────────────────────────────────────
cat("Writing", OUT_TAR, "\n")
system(paste(
  "tar -czf", OUT_TAR, "-C", SCAN_DIR,
  paste0(SCAN, ".snp.wald.png"),
  paste0(SCAN, ".snp_scan.txt")
))
cat("Done.\n")
