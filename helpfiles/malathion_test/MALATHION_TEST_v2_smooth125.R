#!/usr/bin/env Rscript
# MALATHION_TEST_v2_smooth125.R
#
# Publication figures for malathion haplotype scan — 125 kb smoothing.
#
# Run from XQTL2 project root:
#   Rscript helpfiles/malathion_test/MALATHION_TEST_v2_smooth125.R

SCAN_DIR  <- "process/malathion_test/MALATHION_TEST_v2_smooth125"
SCAN      <- "MALATHION_TEST_v2_smooth125"
OUT_TAR   <- file.path(SCAN_DIR, paste0(SCAN, ".hap.tar.gz"))

# ── Figure 1: Haplotype Wald -log10(p) ──────────────────────────────────────
system(paste(
  "Rscript scripts/plot_pseudoscan.R",
  "--scan",   file.path(SCAN_DIR, paste0(SCAN, ".scan.txt")),
  "--out",    file.path(SCAN_DIR, paste0(SCAN, ".wald.png")),
  "--format", "powerpoint",
  "--threshold", "10"
))

# ── Figure 2: Falconer + Cutler H² overlaid ─────────────────────────────────
system(paste(
  "Rscript scripts/plot_H2_overlay.R",
  "--scan",   file.path(SCAN_DIR, paste0(SCAN, ".scan.txt")),
  "--out",    file.path(SCAN_DIR, paste0(SCAN, ".H2.png")),
  "--format", "powerpoint"
))

# ── Bundle for download ─────────────────────────────────────────────────────
cat("Writing", OUT_TAR, "\n")
system(paste(
  "tar -czf", OUT_TAR, "-C", SCAN_DIR,
  paste0(SCAN, ".wald.png"),
  paste0(SCAN, ".H2.png"),
  paste0(SCAN, ".scan.txt"),
  paste0(SCAN, ".meansBySample.txt")
))
cat("Done.\n")
