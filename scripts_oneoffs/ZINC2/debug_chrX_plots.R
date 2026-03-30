#!/usr/bin/env Rscript
###############################################################################
# debug_chrX_plots.R — Generate diagnostic plots for chrX debug
#
# Produces:
#   chrX_freq.png       — all TRT x REP overlaid (C=grey, Z=red)
#   chrX_freq_single.png — single replicate (C REP=1) to judge smoothness
#   chrX_wald.png       — Wald scores from hap_scan (if scan file exists)
#
# Usage:
#   Rscript scripts_oneoffs/ZINC2/debug_chrX_plots.R
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

DIR     <- "process/ZINC2"
SCAN    <- "ZINC2_F_v3"
RESDIR  <- "scripts_oneoffs/ZINC2/debug_chrX_results"
MEANS   <- file.path(DIR, SCAN, paste0(SCAN, ".meansBySample.chrX.txt"))
SCANFILE <- file.path(DIR, SCAN, paste0(SCAN, ".scan.chrX.txt"))

# ── Frequency plots ──────────────────────────────────────────────────────────
cat("Reading", MEANS, "\n")
df <- read.table(MEANS, header = TRUE) %>%
  as_tibble() %>%
  mutate(pos_mb = pos / 1e6,
         TRT    = factor(TRT, levels = c("C", "Z")))

cat(sprintf("  %d rows | %d neg | min=%.6f\n",
    nrow(df), sum(df$freq < 0, na.rm = TRUE), min(df$freq, na.rm = TRUE)))

# Plot 1: all replicates
p1 <- ggplot(df, aes(x = pos_mb, y = freq, group = interaction(TRT, REP),
                      colour = TRT, alpha = TRT)) +
  geom_line(linewidth = 0.3) +
  scale_colour_manual(values = c(C = "grey60", Z = "#CC3333"), name = "Treatment") +
  scale_alpha_manual(values = c(C = 0.5, Z = 0.7), guide = "none") +
  facet_wrap(~founder, ncol = 2, scales = "free_y") +
  labs(x = "Position (Mb)", y = "Smoothed frequency",
       title = "ZINC2_F_v3: smoothed founder frequencies on chrX",
       subtitle = "All replicates; C=grey, Z=red") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(RESDIR, "chrX_freq.png"), p1, width = 10, height = 10, dpi = 150)
cat("Written: chrX_freq.png\n")

# Plot 2: single replicate for smoothness assessment
df_single <- df %>% filter(TRT == "C", REP == 1)

p2 <- ggplot(df_single, aes(x = pos_mb, y = freq)) +
  geom_line(linewidth = 0.4, colour = "black") +
  facet_wrap(~founder, ncol = 2, scales = "free_y") +
  labs(x = "Position (Mb)", y = "Smoothed frequency",
       title = "ZINC2_F_v3: chrX founder frequencies (C, REP 1)",
       subtitle = "Single replicate to assess smoothness") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(RESDIR, "chrX_freq_single.png"), p2, width = 10, height = 10, dpi = 150)
cat("Written: chrX_freq_single.png\n")

# ── Wald plot ─────────────────────────────────────────────────────────────────
if (file.exists(SCANFILE)) {
  cat("Reading", SCANFILE, "\n")
  sc <- read.table(SCANFILE, header = TRUE) %>%
    as_tibble() %>%
    mutate(pos_mb = pos / 1e6)

  p3 <- ggplot(sc, aes(x = pos_mb, y = Wald_log10p)) +
    geom_line(linewidth = 0.3) +
    labs(x = "Position (Mb)", y = "-log10(p)",
         title = "ZINC2_F_v3: Wald scores on chrX") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(RESDIR, "chrX_wald.png"), p3, width = 10, height = 5, dpi = 150)
  cat("Written: chrX_wald.png\n")
} else {
  cat("No scan file found at", SCANFILE, "— skipping Wald plot\n")
}

cat("Done.\n")
