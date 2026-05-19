#!/usr/bin/env Rscript
###############################################################################
# smooth_r2_diag.R  —  R² between raw and smoothed haplotype frequencies,
#                       per pool per founder, across all chromosomes.
#
# Mean R² is the correction factor for converting smoothed Wald statistics
# to calibrated p-values:
#
#   corrected_log10p = -log10(pchisq(tstat / mean_R2, df, lower.tail=FALSE))
#
# Usage:
#   Rscript scripts/temp/smooth_r2_diag.R \
#       --hapsdir  process/PROJ \
#       --smoothdir process/PROJ/SCAN \
#       --scan     SCAN_NAME \
#       --rfile    helpfiles/PROJ/design.txt
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

args   <- commandArgs(trailingOnly = TRUE)
parsed <- list()
i <- 1L
while (i <= length(args)) {
  switch(args[i],
    "--hapsdir"   = { parsed$hapsdir   <- args[i+1]; i <- i+2L },
    "--smoothdir" = { parsed$smoothdir <- args[i+1]; i <- i+2L },
    "--scan"      = { parsed$scan      <- args[i+1]; i <- i+2L },
    "--rfile"     = { parsed$rfile     <- args[i+1]; i <- i+2L },
    stop(paste("Unknown argument:", args[i]))
  )
}

design.df <- read.table(parsed$rfile, header = TRUE)
chrs      <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

r2_all <- vector("list", length(chrs))

for (ci in seq_along(chrs)) {
  mychr   <- chrs[ci]
  sexlink <- if (mychr == "chrX") 0.75 else 1.0

  haps_file   <- file.path(parsed$hapsdir,
                            paste0("R.haps.", mychr, ".out.rds"))
  smooth_file <- file.path(parsed$smoothdir,
                            paste0(parsed$scan, ".smooth.", mychr, ".rds"))

  if (!file.exists(haps_file))   { cat("Missing:", haps_file,   "\n"); next }
  if (!file.exists(smooth_file)) { cat("Missing:", smooth_file, "\n"); next }

  cat(sprintf("Processing %s ...\n", mychr))

  xx1 <- readRDS(haps_file)

  freq_raw <- xx1 %>%
    select(CHROM, pos, sample, Haps, Names, Groups) %>%
    unnest(c(sample, Haps, Names, Groups)) %>%
    unnest(c(Haps, Names, Groups)) %>%
    rename(pool = sample, freq_raw = Haps, founder = Names, group = Groups) %>%
    left_join(design.df, by = c("pool" = "bam")) %>%
    filter(!is.na(TRT)) %>%
    group_by(CHROM, pos, TRT, REP, founder) %>%
    summarise(freq_raw = mean(freq_raw, na.rm = TRUE), .groups = "drop")

  sm          <- readRDS(smooth_file)
  freq_smooth <- sm$freq %>%
    select(CHROM, pos, TRT, REP, founder, freq_smooth = freq)

  joined <- inner_join(freq_raw, freq_smooth,
                       by = c("CHROM", "pos", "TRT", "REP", "founder"))

  r2_all[[ci]] <- joined %>%
    group_by(chr = mychr, TRT, REP, founder) %>%
    summarise(R2    = cor(freq_raw, freq_smooth, use = "complete.obs")^2,
              slope = coef(lm(freq_smooth ~ freq_raw))[2],
              .groups = "drop")
}

r2_table <- bind_rows(r2_all)

# ── Per-chromosome summary ────────────────────────────────────────────────────
chr_summary <- r2_table %>%
  group_by(chr) %>%
  summarise(mean_R2 = mean(R2, na.rm = TRUE),
            sd_R2   = sd(R2,   na.rm = TRUE),
            correction_factor = 1 / mean(R2, na.rm = TRUE),
            .groups = "drop")

cat("\n── R² summary by chromosome ─────────────────────────────────────────\n")
print(chr_summary %>% mutate(across(where(is.numeric), ~round(.x, 3))), n = Inf)

# ── Overall ───────────────────────────────────────────────────────────────────
overall_r2 <- mean(r2_table$R2, na.rm = TRUE)
cat(sprintf("\nOverall mean R²:   %.3f\n", overall_r2))
cat(sprintf("Correction factor: %.3f  (= 1 / mean R²)\n", 1 / overall_r2))
cat(sprintf("\nIn hap_scan.R, replace:\n"))
cat(sprintf("  -log10(pchisq(tstat, df, lower.tail=FALSE))\n"))
cat(sprintf("with:\n"))
cat(sprintf("  -log10(pchisq(tstat / %.3f, df, lower.tail=FALSE))\n", overall_r2))

# ── Full table ────────────────────────────────────────────────────────────────
cat("\n── R² per chromosome × TRT × REP × founder ──────────────────────────\n")
print(r2_table %>% mutate(across(where(is.numeric), ~round(.x, 3))), n = Inf)
