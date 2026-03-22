#!/usr/bin/env Rscript
###############################################################################
# smooth_haps.R  —  Step 1 of the freqsmooth pipeline
#
# Reads R.haps.<chr>.out.rds. Unnests haplotype frequencies and reconstruction
# covariance matrices to long format, applies a running mean across windows
# within each (TRT, REP, founder) / (TRT, REP, i, j) group, then saves:
#   <outdir>/<scan>.smooth.<chr>.rds       -- smoothed data for steps 2 & 3
#   <outdir>/<scan>.meansBySample.<chr>.txt -- smoothed per-founder frequencies
#
# Usage:
#   Rscript scripts/smooth_haps.R \
#       --chr chrX --dir process/proj --outdir SCAN_NAME \
#       --rfile helpfiles/proj/design.txt --smooth-kb 125
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
    stop(paste("Unknown argument:", args[i]))
  )
}

mychr     <- parsed$chr
smooth_kb <- parsed$smooth_kb
design.df <- read.table(parsed$rfile, header = TRUE)

dirout        <- file.path(parsed$dir, parsed$outdir)
fileout_rds   <- file.path(dirout, paste0(parsed$outdir, ".smooth.", mychr, ".rds"))
fileout_means <- file.path(dirout, paste0(parsed$outdir, ".meansBySample.", mychr, ".txt"))
dir.create(dirout, showWarnings = FALSE, recursive = TRUE)

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

# ── Load ──────────────────────────────────────────────────────────────────────
filein <- file.path(parsed$dir, paste0("R.haps.", mychr, ".out.rds"))
cat("Reading", filein, "\n")
xx1 <- readRDS(filein)

step_bp     <- as.integer(median(diff(xx1$pos), na.rm = TRUE))
smooth_half <- round(smooth_kb * 1000L / step_bp)
sexlink     <- if (mychr == "chrX") 0.75 else 1.0

cat(sprintf("  %d windows | step %d bp | smooth_half %d windows (+/-%d kb)\n",
            nrow(xx1), step_bp, smooth_half, smooth_kb))

options(dplyr.summarise.inform = FALSE)

# ── Smooth haplotype frequencies ──────────────────────────────────────────────
# Unnest: one row per (window, pool, founder)
# Average within (window, TRT, REP, founder) to collapse any technical reps
# Then group_by(TRT, REP, founder) and mutate running_mean across windows
cat("Smoothing frequencies...\n")

freq_smoothed <- xx1 %>%
  select(CHROM, pos, sample, Haps, Names) %>%
  unnest(c(sample, Haps, Names)) %>%
  unnest(c(Haps, Names)) %>%
  rename(pool = sample, freq = Haps, founder = Names) %>%
  left_join(design.df, by = c("pool" = "bam")) %>%
  filter(!is.na(TRT)) %>%
  mutate(Num = sexlink * Num) %>%
  group_by(CHROM, pos, TRT, REP, founder) %>%
  summarize(freq = mean(freq, na.rm = TRUE),
            Num  = mean(Num,  na.rm = TRUE), .groups = "drop") %>%
  arrange(CHROM, pos) %>%
  group_by(TRT, REP, founder) %>%
  mutate(freq = running_mean(freq, smooth_half)) %>%
  ungroup()

# ── Smooth reconstruction covariance matrices ─────────────────────────────────
# Unnest: one row per (window, pool), then expand each matrix to (fi, fj, v) rows
# Average within (window, TRT, REP, fi, fj) then smooth across windows
cat("Smoothing covariance matrices...\n")

err_unnested <- xx1 %>%
  select(CHROM, pos, sample, Err) %>%
  unnest(c(sample, Err)) %>%
  rename(pool = sample) %>%
  left_join(design.df %>% select(bam, TRT, REP), by = c("pool" = "bam")) %>%
  filter(!is.na(TRT))

# Vectorized covariance expansion — column-major order matches as.vector()
nF_cov  <- nrow(as.matrix(err_unnested$Err[[1]]))
nF2     <- nF_cov^2L
fi_tmpl <- rep(seq_len(nF_cov), nF_cov)               # row indices
fj_tmpl <- rep(seq_len(nF_cov), each = nF_cov)        # col indices
v_mat   <- do.call(rbind, lapply(err_unnested$Err,     # n_rows x nF2
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

# ── Save smoothed data for steps 2 and 3 ─────────────────────────────────────
founder_names <- sort(unique(freq_smoothed$founder))
nrepl         <- n_distinct(freq_smoothed$REP)

cat("Writing smoothed RDS:", fileout_rds, "\n")
saveRDS(
  list(freq          = freq_smoothed,
       err           = err_smoothed,
       founder_names = founder_names,
       nrepl         = nrepl),
  fileout_rds
)

# ── meansBySample from smoothed frequencies ───────────────────────────────────
cat("Writing meansBySample:", fileout_means, "\n")
freq_smoothed %>%
  select(chr = CHROM, pos, TRT, REP, founder, freq) %>%
  filter(!is.na(freq)) %>%
  write.table(fileout_means)

cat("Done.\n")
