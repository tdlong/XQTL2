#!/usr/bin/env Rscript
###############################################################################
# debug_chrX_diag.R — Diagnostics for chrX gap-fill fix
#
# Produces:
#   1. Text dump of B2/B4 freqs at base of chrX, all TRT x REP
#   2. Smoothness stats per founder x TRT x REP
#   3. Plot: smoothed founder frequencies at base of chrX (all reps, C vs Z)
#   4. Plot: Wald scores (Manhattan) for chrX
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

DIR     <- "process/ZINC2/ZINC2_F_v3"
MEANS   <- file.path(DIR, "ZINC2_F_v3.meansBySample.chrX.txt")
SCAN    <- file.path(DIR, "ZINC2_F_v3.scan.chrX.txt")
OUTDIR  <- "scripts_oneoffs/ZINC2/debug_chrX_results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUTTEXT <- file.path(OUTDIR, "diag.txt")
sink(OUTTEXT, split = TRUE)
options(width = 200)

# ── Read means ───────────────────────────────────────────────────────────────
cat("Reading", MEANS, "\n")
df <- read.table(MEANS, header = TRUE) %>% as_tibble()
cat(sprintf("  %d rows | %d positions | founders: %s\n",
            nrow(df), n_distinct(df$pos),
            paste(sort(unique(df$founder)), collapse = ", ")))
cat(sprintf("  TRTs: %s | REPs: %s\n",
            paste(sort(unique(df$TRT)), collapse = ", "),
            paste(sort(unique(df$REP)), collapse = ", ")))

# ── B2 and B4 at base of chrX, all TRT x REP ────────────────────────────────
cat("\n=== B2 and B4 frequencies, base of chrX (pos > 20 Mb), all TRT x REP ===\n")
base_df <- df %>%
  filter(founder %in% c("B2", "B4"), pos > 20e6) %>%
  arrange(founder, TRT, REP, pos) %>%
  mutate(pos_mb = pos / 1e6)

for (f in c("B2", "B4")) {
  for (trt in sort(unique(base_df$TRT))) {
    for (rep in sort(unique(base_df$REP))) {
      tmp <- base_df %>% filter(founder == f, TRT == trt, REP == rep)
      cat(sprintf("\n--- %s  TRT=%s  REP=%d  (%d windows) ---\n", f, trt, rep, nrow(tmp)))
      step <- max(1L, nrow(tmp) %/% 25L)
      idx <- seq(1, nrow(tmp), by = step)
      for (i in idx) {
        cat(sprintf("  pos=%6.2fMb  freq=%.4f\n", tmp$pos_mb[i], tmp$freq[i]))
      }
    }
  }
}

# ── Wide format at base ─────────────────────────────────────────────────────
cat("\n=== All founders wide, base of chrX (pos > 20 Mb) ===\n")
for (trt in sort(unique(df$TRT))) {
  for (rep in c(1, 2)) {
    cat(sprintf("\n--- TRT=%s  REP=%d ---\n", trt, rep))
    wide <- df %>%
      filter(pos > 20e6, REP == rep, TRT == trt) %>%
      select(pos, founder, freq) %>%
      pivot_wider(names_from = founder, values_from = freq) %>%
      arrange(pos)
    step2 <- max(1L, nrow(wide) %/% 25L)
    idx2 <- seq(1, nrow(wide), by = step2)
    print(wide[idx2, ], n = 100)
  }
}

# ── Smoothness check ────────────────────────────────────────────────────────
cat("\n=== Smoothness check: mean |Δfreq| between consecutive windows ===\n")
smooth_check <- df %>%
  arrange(founder, TRT, REP, pos) %>%
  group_by(founder, TRT, REP) %>%
  summarize(
    n_windows     = n(),
    n_na          = sum(is.na(freq)),
    mean_freq     = mean(freq, na.rm = TRUE),
    mean_abs_diff = mean(abs(diff(freq)), na.rm = TRUE),
    max_abs_diff  = max(abs(diff(freq)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_diff))
print(as.data.frame(smooth_check), row.names = FALSE)

cat("\nDone with text diagnostics.\n")
sink()

# ── Plot 1: founder frequencies at base of chrX ─────────────────────────────
cat("Generating frequency plot...\n")
plot_df <- df %>%
  filter(pos > 18e6) %>%
  mutate(pos_mb = pos / 1e6, TRT = factor(TRT, levels = c("C", "Z")))

p1 <- ggplot(plot_df, aes(x = pos_mb, y = freq,
                           group = interaction(TRT, REP),
                           colour = TRT, alpha = TRT)) +
  geom_line(linewidth = 0.3) +
  scale_colour_manual(values = c(C = "grey60", Z = "#CC3333"), name = "Treatment") +
  scale_alpha_manual(values = c(C = 0.5, Z = 0.7), guide = "none") +
  facet_wrap(~founder, ncol = 2, scales = "free_y") +
  labs(x = "Position (Mb)", y = "Smoothed frequency",
       title = "ZINC2_F_v3: smoothed founder frequencies, chrX > 18Mb",
       subtitle = "All replicates; C=grey, Z=red") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(OUTDIR, "freq_base_chrX.png"), p1, width = 10, height = 10, dpi = 150)

# ── Plot 2: Wald Manhattan for chrX ─────────────────────────────────────────
cat("Generating Wald plot...\n")
if (file.exists(SCAN)) {
  scan_df <- read.table(SCAN, header = TRUE) %>%
    as_tibble() %>%
    mutate(pos_mb = pos / 1e6)

  p2 <- ggplot(scan_df, aes(x = pos_mb, y = Wald_log10p)) +
    geom_point(size = 0.3, alpha = 0.5) +
    geom_hline(yintercept = c(5, 10), linetype = "dashed", colour = "red", linewidth = 0.3) +
    labs(x = "Position (Mb)", y = "-log10(p)",
         title = "ZINC2_F_v3: Wald test, chrX",
         subtitle = sprintf("%d windows", nrow(scan_df))) +
    theme_bw(base_size = 10)

  ggsave(file.path(OUTDIR, "wald_chrX.png"), p2, width = 10, height = 4, dpi = 150)

  # Also zoom into the base
  p3 <- p2 %+% filter(scan_df, pos_mb > 18) +
    labs(title = "ZINC2_F_v3: Wald test, chrX > 18Mb (base)")
  ggsave(file.path(OUTDIR, "wald_base_chrX.png"), p3, width = 10, height = 4, dpi = 150)
} else {
  cat("WARNING: scan file not found:", SCAN, "\n")
}

cat("All results in:", OUTDIR, "\n")
