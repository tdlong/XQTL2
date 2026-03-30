#!/usr/bin/env Rscript
###############################################################################
# smooth_haps_debug.R — Instrumented copy of scripts/smooth_haps.R
#
# Identical logic to smooth_haps.R, but writes diagnostic output to a log file
# after every pipeline step so we can see exactly where negatives appear.
#
# Diagnostics written to: <outdir>/<scan>.smooth_diag.<chr>.txt
#
# Usage: same as smooth_haps.R
#   Rscript scripts_oneoffs/ZINC2/smooth_haps_debug.R \
#       --chr chrX --dir process/ZINC2 --outdir ZINC2_F_v3 \
#       --rfile helpfiles/ZINC2/Zinc2.test.F.txt --smooth-kb 250
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args   <- commandArgs(trailingOnly = TRUE)
parsed <- list(smooth_kb = 125L)
i <- 1L
while (i <= length(args)) {
  switch(args[i],
    "--chr"       = { parsed$chr       <- args[i+1]; i <- i+2L },
    "--dir"       = { parsed$dir       <- args[i+1]; i <- i+2L },
    "--outdir"    = { parsed$outdir    <- args[i+1]; i <- i+2L },
    "--rfile"     = { parsed$rfile     <- args[i+1]; i <- i+2L },
    "--smooth-kb" = { parsed$smooth_kb <- as.integer(args[i+1]); i <- i+2L },
    "--diagdir"   = { parsed$diagdir   <- args[i+1]; i <- i+2L },
    stop(paste("Unknown argument:", args[i]))
  )
}

mychr     <- parsed$chr
smooth_kb <- parsed$smooth_kb
design.df <- read.table(parsed$rfile, header = TRUE)

dirout        <- file.path(parsed$dir, parsed$outdir)
fileout_rds   <- file.path(dirout, paste0(parsed$outdir, ".smooth.", mychr, ".rds"))
fileout_means <- file.path(dirout, paste0(parsed$outdir, ".meansBySample.", mychr, ".txt"))

# Diagnostics go to --diagdir (trackable by git), not process/ (gitignored)
diagdir <- if (!is.null(parsed$diagdir)) parsed$diagdir else dirout
dir.create(diagdir, showWarnings = FALSE, recursive = TRUE)
fileout_diag  <- file.path(diagdir, paste0(parsed$outdir, ".smooth_diag.", mychr, ".txt"))

dir.create(dirout, showWarnings = FALSE, recursive = TRUE)

# ── Diagnostic helper ─────────────────────────────────────────────────────────
diag_con <- file(fileout_diag, open = "wt")

diag <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")                    # stdout (slurm log)
  writeLines(msg, con = diag_con)   # persistent file
}

report_freq <- function(label, x) {
  diag(sprintf("  [%s] n=%d  NA=%d  NaN=%d  neg=%d  min=%.6f  max=%.6f",
    label, length(x),
    sum(is.na(x)), sum(is.nan(x)),
    sum(x < 0, na.rm = TRUE),
    min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
}

# ── Running mean (edge-aware, O(n), NA-safe) ──────────────────────────────────
running_mean <- function(x, h) {
  n  <- length(x); if (n == 0L) return(numeric(0))
  ok <- !is.na(x)
  xc <- replace(x, !ok, 0)
  cs <- c(0, cumsum(xc)); cn <- c(0, cumsum(as.numeric(ok)))
  lo <- pmax(1L, seq_len(n) - h); hi <- pmin(n, seq_len(n) + h)
  tot <- cs[hi+1L] - cs[lo]; cnt <- cn[hi+1L] - cn[lo]
  ifelse(cnt > 0L, tot/cnt, NA_real_)
}

# ── Gap filler (mean-anchored interpolation) ─────────────────────────────────
fill_gaps <- function(x, h) {
  n  <- length(x)
  if (n == 0L || !anyNA(x)) return(x)

  ok  <- !is.na(x)
  if (!any(ok)) return(x)

  rle_na  <- rle(!ok)
  ends    <- cumsum(rle_na$lengths)
  starts  <- ends - rle_na$lengths + 1L

  for (g in which(rle_na$values)) {
    ga <- starts[g]
    gb <- ends[g]

    left_idx <- which(ok & seq_along(x) < ga)
    if (length(left_idx) > 0L) {
      anchor_L <- mean(x[tail(left_idx, h)])
    } else {
      anchor_L <- NULL
    }

    right_idx <- which(ok & seq_along(x) > gb)
    if (length(right_idx) > 0L) {
      anchor_R <- mean(x[head(right_idx, h)])
    } else {
      anchor_R <- NULL
    }

    gap_len <- gb - ga + 1L
    if (!is.null(anchor_L) && !is.null(anchor_R)) {
      x[ga:gb] <- seq(anchor_L, anchor_R, length.out = gap_len)
    } else if (!is.null(anchor_L)) {
      x[ga:gb] <- anchor_L
    } else if (!is.null(anchor_R)) {
      x[ga:gb] <- anchor_R
    }
  }
  x
}

# ── Load ──────────────────────────────────────────────────────────────────────
filein <- file.path(parsed$dir, paste0("R.haps.", mychr, ".out.rds"))
diag(paste("Reading", filein))
xx1 <- readRDS(filein)

step_bp     <- as.integer(median(diff(xx1$pos), na.rm = TRUE))
smooth_half <- round(smooth_kb * 1000L / step_bp)
sexlink     <- if (mychr == "chrX") 0.75 else 1.0

diag(sprintf("  %d windows | step %d bp | smooth_half %d windows (+/-%d kb)",
            nrow(xx1), step_bp, smooth_half, smooth_kb))

options(dplyr.summarise.inform = FALSE)

# ── STEP A: Unnest ────────────────────────────────────────────────────────────
diag("Step A: unnest")

freq_raw <- xx1 %>%
  select(CHROM, pos, sample, Haps, Names, Groups) %>%
  unnest(c(sample, Haps, Names, Groups)) %>%
  unnest(c(Haps, Names, Groups)) %>%
  rename(pool = sample, freq = Haps, founder = Names, group = Groups) %>%
  left_join(design.df, by = c("pool" = "bam")) %>%
  filter(!is.na(TRT)) %>%
  mutate(Num = sexlink * Num)

report_freq("after unnest", freq_raw$freq)

# ── STEP B: Mask ──────────────────────────────────────────────────────────────
diag("Step B: mask unresolvable founders")

freq_raw <- freq_raw %>%
  group_by(CHROM, pos, pool, group) %>%
  mutate(group_size = n()) %>%
  ungroup() %>%
  mutate(freq = if_else(group_size > 1L, NA_real_, freq))

n_masked <- sum(is.na(freq_raw$freq))
n_total  <- nrow(freq_raw)
diag(sprintf("  Masked %d / %d founder-window estimates (%.1f%%) as unresolvable",
            n_masked, n_total, 100 * n_masked / n_total))
report_freq("after mask", freq_raw$freq)

# ── STEP C: Summarize (pool mean) ────────────────────────────────────────────
diag("Step C: summarize (collapse pools per TRT/REP/founder)")

freq_summ <- freq_raw %>%
  group_by(CHROM, pos, TRT, REP, founder) %>%
  summarize(freq = mean(freq, na.rm = TRUE),
            Num  = mean(Num,  na.rm = TRUE), .groups = "drop") %>%
  arrange(CHROM, pos)

report_freq("after summarize", freq_summ$freq)

# Check for NaN (all-NA groups produce NaN from mean(na.rm=TRUE))
n_nan <- sum(is.nan(freq_summ$freq))
if (n_nan > 0) {
  diag(sprintf("  WARNING: %d NaN values from all-NA groups (converting to NA)", n_nan))
  freq_summ$freq[is.nan(freq_summ$freq)] <- NA_real_
  report_freq("after NaN->NA fix", freq_summ$freq)
}

# Clamp: frequencies can't be negative (lsei solver floor is 0.0003,
# but numerical issues occasionally produce negatives)
n_neg_pre <- sum(freq_summ$freq < 0, na.rm = TRUE)
if (n_neg_pre > 0) {
  diag(sprintf("  Clamping %d negative frequencies to 0.0003", n_neg_pre))
  freq_summ$freq <- pmax(freq_summ$freq, 0.0003, na.rm = TRUE)
  report_freq("after clamp", freq_summ$freq)
}

# ── Dump B2/B4 at base of chrX for comparison with standalone replay ─────────
diag("")
diag("=== B2/B4 detail at base of chrX (21-22.5 Mb), REP=1, TRT=C ===")
diag("--- After summarize (before fill_gaps) ---")
spot <- freq_summ %>%
  filter(founder %in% c("B2", "B4"), pos > 21e6, pos < 22.5e6,
         REP == 1, TRT == "C") %>%
  arrange(founder, pos)
for (r in seq_len(nrow(spot))) {
  diag(sprintf("  %s pos=%d freq=%.6f", spot$founder[r], spot$pos[r], spot$freq[r]))
}

# ── STEP D: Fill gaps ────────────────────────────────────────────────────────
diag("")
diag("Step D: fill_gaps")

freq_filled <- freq_summ %>%
  group_by(TRT, REP, founder) %>%
  mutate(freq = fill_gaps(freq, smooth_half)) %>%
  ungroup()

report_freq("after fill_gaps", freq_filled$freq)

diag("--- B2/B4 at base after fill_gaps ---")
spot2 <- freq_filled %>%
  filter(founder %in% c("B2", "B4"), pos > 21e6, pos < 22.5e6,
         REP == 1, TRT == "C") %>%
  arrange(founder, pos)
for (r in seq_len(nrow(spot2))) {
  diag(sprintf("  %s pos=%d freq=%.6f", spot2$founder[r], spot2$pos[r], spot2$freq[r]))
}

# ── STEP E: Running mean ────────────────────────────────────────────────────
diag("")
diag("Step E: running_mean")

freq_smoothed <- freq_filled %>%
  group_by(TRT, REP, founder) %>%
  mutate(freq = running_mean(freq, smooth_half)) %>%
  ungroup()

report_freq("after running_mean", freq_smoothed$freq)

diag("--- B2/B4 at base after running_mean ---")
spot3 <- freq_smoothed %>%
  filter(founder %in% c("B2", "B4"), pos > 21e6, pos < 22.5e6,
         REP == 1, TRT == "C") %>%
  arrange(founder, pos)
for (r in seq_len(nrow(spot3))) {
  diag(sprintf("  %s pos=%d freq=%.6f", spot3$founder[r], spot3$pos[r], spot3$freq[r]))
}

# ── Per-founder negative summary ──────────────────────────────────────────────
diag("")
diag("=== Per-founder negative count in final output ===")
neg_by_founder <- freq_smoothed %>%
  filter(freq < 0) %>%
  count(founder) %>%
  arrange(desc(n))
if (nrow(neg_by_founder) > 0) {
  for (r in seq_len(nrow(neg_by_founder))) {
    diag(sprintf("  %s: %d negatives", neg_by_founder$founder[r], neg_by_founder$n[r]))
  }
} else {
  diag("  NONE")
}

# ── Smooth reconstruction covariance matrices ─────────────────────────────────
diag("")
diag("Smoothing covariance matrices...")

err_unnested <- xx1 %>%
  select(CHROM, pos, sample, Err) %>%
  unnest(c(sample, Err)) %>%
  rename(pool = sample) %>%
  left_join(design.df %>% select(bam, TRT, REP), by = c("pool" = "bam")) %>%
  filter(!is.na(TRT))

nF_cov  <- nrow(as.matrix(err_unnested$Err[[1]]))
nF2     <- nF_cov^2L
fi_tmpl <- rep(seq_len(nF_cov), nF_cov)
fj_tmpl <- rep(seq_len(nF_cov), each = nF_cov)
v_mat   <- do.call(rbind, lapply(err_unnested$Err,
             function(m) as.vector(as.matrix(m))))

err_smoothed <- err_unnested %>%
  select(-Err) %>%
  tidyr::uncount(nF2) %>%
  mutate(fi = rep(fi_tmpl, nrow(err_unnested)),
         fj = rep(fj_tmpl, nrow(err_unnested)),
         v  = as.vector(t(v_mat))) %>%
  group_by(CHROM, pos, TRT, REP, fi, fj) %>%
  summarize(v = mean(v, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHROM, pos) %>%
  group_by(TRT, REP, fi, fj) %>%
  mutate(v = running_mean(v, smooth_half)) %>%
  ungroup()

# ── Save ──────────────────────────────────────────────────────────────────────
founder_names <- sort(unique(freq_smoothed$founder))
nrepl         <- n_distinct(freq_smoothed$REP)

diag(paste("Writing smoothed RDS:", fileout_rds))
saveRDS(
  list(freq          = freq_smoothed,
       err           = err_smoothed,
       founder_names = founder_names,
       nrepl         = nrepl),
  fileout_rds
)

diag(paste("Writing meansBySample:", fileout_means))
freq_smoothed %>%
  select(chr = CHROM, pos, TRT, REP, founder, freq) %>%
  filter(!is.na(freq)) %>%
  write.table(fileout_means)

diag(paste("Diagnostics written to:", fileout_diag))
diag("Done.")
close(diag_con)
