#!/usr/bin/env Rscript
# MALATHION_TEST_v2_smooth125_snp.R
#
# Publication figure for the malathion SNP scan.
# Produces three 5-panel (per-chromosome) figures using plot_freqsmooth_snp.R,
# then bundles all outputs into a single tar for download.
#
# Run from XQTL2 project root:
#   Rscript helpfiles/malathion_test/MALATHION_TEST_v2_smooth125_snp.R

SCAN_DIR  <- "process/malathion_test/MALATHION_TEST_v2_smooth125"
SCAN      <- "MALATHION_TEST_v2_smooth125"
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

# ── Figure 2: Falconer H² ─────────────────────────────────────────────────────
SCAN_FILES   <- SCAN_FILE
SCAN_LABELS  <- NULL
SCAN_COLOURS <- c("#E31A1C")
YVAR         <- "Falc_H2"
YLAB         <- "Falconer H\u00b2"
THRESHOLD    <- NULL
OUT_FILE     <- file.path(SCAN_DIR, paste0(SCAN, ".snp.falcH2.png"))
FORMAT       <- "powerpoint"
source("scripts/plot_freqsmooth_snp.R")

# ── Figure 3: Cutler H² ───────────────────────────────────────────────────────
SCAN_FILES   <- SCAN_FILE
SCAN_LABELS  <- NULL
SCAN_COLOURS <- c("#33A02C")
YVAR         <- "Cutl_H2"
YLAB         <- "Cutler H\u00b2"
THRESHOLD    <- NULL
OUT_FILE     <- file.path(SCAN_DIR, paste0(SCAN, ".snp.cutlH2.png"))
FORMAT       <- "powerpoint"
source("scripts/plot_freqsmooth_snp.R")

# ── Bundle for download ────────────────────────────────────────────────────────
cat("Writing", OUT_TAR, "\n")
system(paste(
  "tar -czf", OUT_TAR, "-C", SCAN_DIR,
  paste0(SCAN, ".snp.wald.png"),
  paste0(SCAN, ".snp.falcH2.png"),
  paste0(SCAN, ".snp.cutlH2.png"),
  paste0(SCAN, ".snp_scan.txt")
))
cat("Done.\n")
