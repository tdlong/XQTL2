#!/usr/bin/env Rscript
###############################################################################
# chrX_diag.R — Diagnostic: dump smoothed founder freqs for B2 and B4
#               at the base of chrX, plus check smoothness across all founders.
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

MEANS <- "process/ZINC2/ZINC2_F_v3/ZINC2_F_v3.meansBySample.chrX.txt"
OUTFILE <- "output/chrX_diag.txt"
sink(OUTFILE, split = TRUE)   # write to file AND stdout

cat("Reading", MEANS, "\n")
df <- read.table(MEANS, header=TRUE) %>% as_tibble()

cat(sprintf("  %d rows | %d positions | founders: %s\n",
            nrow(df), n_distinct(df$pos),
            paste(sort(unique(df$founder)), collapse=", ")))

# ── Focus on base of chrX (first 3 Mb) for B2 and B4 ────────────────────────
cat("\n=== B2 and B4 frequencies, base of chrX (pos < 3 Mb), REP=1, TRT=C ===\n")
base_df <- df %>%
  filter(founder %in% c("B2","B4"), pos < 3e6, REP == 1, TRT == "C") %>%
  arrange(founder, pos) %>%
  select(founder, pos, freq) %>%
  mutate(pos_kb = pos / 1000)

for (f in c("B2", "B4")) {
  cat(sprintf("\n--- %s ---\n", f))
  tmp <- base_df %>% filter(founder == f)
  cat(sprintf("  %d windows\n", nrow(tmp)))
  # Print every ~5th position to keep it readable
  step <- max(1L, nrow(tmp) %/% 40L)
  idx <- seq(1, nrow(tmp), by = step)
  for (i in idx) {
    cat(sprintf("  pos=%7.1fkb  freq=%.4f\n", tmp$pos_kb[i], tmp$freq[i]))
  }
}

# ── Same region, all founders, single rep — wide format ──────────────────────
cat("\n=== All founders, base of chrX (pos < 3 Mb), REP=1, TRT=C ===\n")
wide <- df %>%
  filter(pos < 3e6, REP == 1, TRT == "C") %>%
  select(pos, founder, freq) %>%
  pivot_wider(names_from = founder, values_from = freq) %>%
  arrange(pos)

# Print every ~10th row
step2 <- max(1L, nrow(wide) %/% 30L)
idx2 <- seq(1, nrow(wide), by = step2)
cat(sprintf("  %d total windows, showing every %dth\n\n", nrow(wide), step2))
options(width = 200)
print(wide[idx2, ], n = 100)

# ── Smoothness check: mean absolute first-difference per founder ─────────────
cat("\n=== Smoothness check: mean |Δfreq| between consecutive windows ===\n")
cat("  (smaller = smoother; computed over all positions, REP=1, TRT=C)\n\n")
smooth_check <- df %>%
  filter(REP == 1, TRT == "C") %>%
  arrange(founder, pos) %>%
  group_by(founder) %>%
  summarize(
    n_windows   = n(),
    n_na        = sum(is.na(freq)),
    mean_freq   = mean(freq, na.rm = TRUE),
    mean_abs_diff = mean(abs(diff(freq)), na.rm = TRUE),
    max_abs_diff  = max(abs(diff(freq)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_diff))

print(as.data.frame(smooth_check), row.names = FALSE)

cat("\nDone.\n")
sink()
cat("Results written to:", OUTFILE, "\n")
# Commit results back so they can be read remotely
system(sprintf("cd /dfs7/adl/tdlong/fly_pool/XQTL2 && git add %s && git commit -m 'chrX_diag results' && git push dev dev 2>&1", OUTFILE))
