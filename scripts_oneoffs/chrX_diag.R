#!/usr/bin/env Rscript
###############################################################################
# chrX_diag.R — Diagnostic: dump smoothed founder freqs for B2 and B4
#               at the base of chrX (>20Mb), all TRT x REP combos.
#               Check smoothness across all founders.
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

MEANS <- "process/ZINC2/ZINC2_F_v3/ZINC2_F_v3.meansBySample.chrX.txt"
OUTFILE <- "scripts_oneoffs/chrX_diag_results.txt"
sink(OUTFILE, split = TRUE)

cat("Reading", MEANS, "\n")
df <- read.table(MEANS, header=TRUE) %>% as_tibble()

cat(sprintf("  %d rows | %d positions | founders: %s\n",
            nrow(df), n_distinct(df$pos),
            paste(sort(unique(df$founder)), collapse=", ")))
cat(sprintf("  TRTs: %s | REPs: %s\n",
            paste(sort(unique(df$TRT)), collapse=", "),
            paste(sort(unique(df$REP)), collapse=", ")))

options(width = 200)

# ── B2 and B4 at base of chrX, ALL TRT x REP combos ─────────────────────────
cat("\n=== B2 and B4 frequencies, base of chrX (pos > 20 Mb), all TRT x REP ===\n")
base_df <- df %>%
  filter(founder %in% c("B2","B4"), pos > 20e6) %>%
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

# ── Wide format: all founders, all TRT x REP at base ────────────────────────
cat("\n=== All founders wide, base of chrX (pos > 20 Mb) ===\n")
for (trt in sort(unique(df$TRT))) {
  for (rep in sort(unique(df$REP))) {
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

# ── Smoothness check: all TRT x REP x founder ───────────────────────────────
cat("\n=== Smoothness check: mean |Δfreq| between consecutive windows ===\n")
cat("  (smaller = smoother; computed over ALL positions)\n\n")
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

cat("\nDone.\n")
sink()
cat("Results written to:", OUTFILE, "\n")
system(sprintf("cd /dfs7/adl/tdlong/fly_pool/XQTL2 && git add %s && git commit -m 'chrX_diag results' && git push dev HEAD 2>&1", OUTFILE))
