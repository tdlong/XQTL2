###############################################################################
# haps2scan.stable.code.R
#
# Stabilized Wald scan with covariance matrix smoothing and eigenvalue
# regularization.  Drop-in replacement for haps2scan.Apr2025.code.R.
#
# Changes from Apr2025:
#   1. Pooled covariance matrices are smoothed across a ±COV_SMOOTH_KB half-window
#   2. Eigenvalues are floored at RIDGE_FRACTION * mean(eigenvalue)
#   3. Frequencies (means) are NOT smoothed
#   4. Pseu_log10p and Heritability computed from raw (unsmoothed) data
#
# Same inputs, same output format — all downstream scripts work unchanged.
###############################################################################

# ── Tuning parameters ────────────────────────────────────────────────────────
RIDGE_FRACTION <- 0.01   # Floor eigenvalues at 1% of mean eigenvalue

# COV_SMOOTH_KB and FREQ_SMOOTH_KB are set by the driver .R file in kb.
# 0 = no smoothing.  Window counts are derived from actual data spacing below.
if (!exists("COV_SMOOTH_KB"))  COV_SMOOTH_KB  <- 0L
if (!exists("FREQ_SMOOTH_KB")) FREQ_SMOOTH_KB <- 0L

# ── Helper: edge-aware running mean (vectorized, O(n)) ──────────────────────
running_mean <- function(x, half_window) {
  n <- length(x)
  if (n == 0) return(numeric(0))
  is_valid <- !is.na(x)
  x_clean  <- replace(x, !is_valid, 0)
  cs <- c(0, cumsum(x_clean))
  cn <- c(0, cumsum(as.numeric(is_valid)))
  lo <- pmax(1L, seq_len(n) - half_window)
  hi <- pmin(n,  seq_len(n) + half_window)
  totals <- cs[hi + 1L] - cs[lo]
  counts <- cn[hi + 1L] - cn[lo]
  ifelse(counts > 0, totals / counts, NA_real_)
}

# ── Helper: extract per-window data ─────────────────────────────────────────
# Same as doscan2 steps 1–4: unnest, join design, check founders,
# aggregate by (TRT, REP).
# Returns list(p1, p2, covar1, covar2, N1, N2, nrepl) or NULL.
extract_window_data <- function(df, chr, Nfounders) {
  sexlink <- if (chr == "chrX") 0.75 else 1

  df2 <- df %>%
    unnest(cols = c(sample, Groups, Haps, Err, Names)) %>%
    left_join(design.df, join_by(sample == bam)) %>%
    filter(!is.na(TRT))

  allFounders <- as.numeric(
    df2 %>% mutate(mm = max(unlist(Groups))) %>% summarize(max(mm))
  )
  if (allFounders != Nfounders) return(NULL)

  df3 <- df2 %>%
    select(-Groups) %>%
    group_by(TRT, REP) %>%
    summarise(
      Err_mean  = list(reduce(map(Err, ~as.matrix(.x)), `+`) / length(Err)),
      Haps_mean = list(reduce(map(Haps, ~as.vector(.x)), `+`) / length(Haps)),
      Names     = list(first(Names)),
      Num_mean  = sexlink * mean(Num)
    ) %>%
    rename(Haps = Haps_mean, Num = Num_mean, Err = Err_mean)

  p1 <- df3 %>% filter(TRT == "C") %>% pull(Haps) %>%
    as.data.frame() %>% as.matrix() %>% t()
  row.names(p1) <- NULL
  p2 <- df3 %>% filter(TRT == "Z") %>% pull(Haps) %>%
    as.data.frame() %>% as.matrix() %>% t()
  row.names(p2) <- NULL
  covar1 <- do.call(abind, c(df3 %>% filter(TRT == "C") %>% pull(Err), along = 3))
  covar2 <- do.call(abind, c(df3 %>% filter(TRT == "Z") %>% pull(Err), along = 3))
  nrepl  <- df3 %>% filter(TRT == "C") %>% nrow()
  N1     <- df3 %>% filter(TRT == "C") %>% pull(Num)
  N2     <- df3 %>% filter(TRT == "Z") %>% pull(Num)

  list(p1 = p1, p2 = p2, covar1 = covar1, covar2 = covar2,
       N1 = N1, N2 = N2, nrepl = nrepl)
}

# ── Helper: pool one window across replicates ────────────────────────────────
# Same as wald.test3 lines 37–59: compute N.eff per replicate,
# weight and sum covariance matrices, pool frequencies.
# Returns list(p1, p2, covar1, covar2) with pooled (vector/matrix) values.
pool_one_window <- function(wd) {
  if (is.null(wd)) return(NULL)

  nrepl  <- wd$nrepl
  p1     <- wd$p1;     p2     <- wd$p2
  covar1 <- wd$covar1; covar2 <- wd$covar2
  N1     <- wd$N1;     N2     <- wd$N2

  lp1    <- ncol(p1)
  N1.eff <- rep(NA_real_, nrepl)
  N2.eff <- rep(NA_real_, nrepl)
  cv1    <- array(NA_real_, c(lp1, lp1, nrepl))
  cv2    <- array(NA_real_, c(lp1, lp1, nrepl))

  for (i in 1:nrepl) {
    covmat1 <- mn.covmat((N1[i] * p1[i,] + N2[i] * p2[i,]) / (N1[i] + N2[i]),
                          2 * N1[i])
    covmat2 <- mn.covmat((N1[i] * p1[i,] + N2[i] * p2[i,]) / (N1[i] + N2[i]),
                          2 * N2[i])

    N1.eff[i] <- sum(diag(covmat1)) * 4 * N1[i]^2 /
      (sum(diag(covmat1)) * 2 * N1[i] + 2 * N1[i] * sum(diag(covar1[,,i])))
    N2.eff[i] <- sum(diag(covmat2)) * 4 * N2[i]^2 /
      (sum(diag(covmat2)) * 2 * N2[i] + 2 * N2[i] * sum(diag(covar2[,,i])))

    cv1[,,i] <- (covmat1 + covar1[,,i]) * N1.eff[i]^2
    cv2[,,i] <- (covmat2 + covar2[,,i]) * N2.eff[i]^2
  }

  p1_pooled     <- as.vector(N1.eff %*% p1 / sum(N1.eff))
  p2_pooled     <- as.vector(N2.eff %*% p2 / sum(N2.eff))
  covar1_pooled <- rowSums(cv1, dims = 2) / sum(N1.eff)^2
  covar2_pooled <- rowSums(cv2, dims = 2) / sum(N2.eff)^2

  list(p1 = p1_pooled, p2 = p2_pooled,
       covar1 = covar1_pooled, covar2 = covar2_pooled)
}

# ── Helper: Wald test from pooled data with eigenvalue regularization ────────
# Like wald.test3 lines 69–81 but floors small eigenvalues.
wald_from_pooled <- function(p1, p2, covar1, covar2) {
  df    <- length(p1) - 1
  covar <- covar1 + covar2
  eg    <- eigen(covar)
  ev    <- eg$vectors[, 1:df]
  eval  <- eg$values[1:df]

  # ── The key fix: regularize eigenvalues ──
  eval <- pmax(eval, RIDGE_FRACTION * mean(eval))

  trafo <- diag(1/sqrt(eval)) %*% t(ev)
  tstat <- sum((trafo %*% (p1 - p2))^2)
  pval  <- exp(pchisq(tstat, df, lower.tail = FALSE, log.p = TRUE))

  list(wald.test = tstat, p.value = pval,
       avg.var = average_variance(covar)$avg_var)
}


###############################################################################
# Main pipeline
###############################################################################

xx1 <- readRDS(filein)
Nfounders <- length(xx1$Groups[[1]][[1]])
ProportionSelect <- design.df %>%
  filter(TRT == "Z") %>% select(REP, Proportion) %>% arrange(REP)

# ── Step 1: Nest by window ───────────────────────────────────────────────────
cat("Step 1: Nesting by window...\n")
windows <- xx1 %>%
  group_by(CHROM, pos) %>%
  nest()
W <- nrow(windows)
cat(sprintf("  %d windows on %s\n", W, mychr))

# Derive window step size from data; convert kb smoothing distances to window counts
step_kb  <- median(diff(sort(unique(windows$pos)))) / 1000
COV_HALF  <- if (COV_SMOOTH_KB  > 0) round(COV_SMOOTH_KB  / step_kb) else 0L
FREQ_HALF <- if (FREQ_SMOOTH_KB > 0) round(FREQ_SMOOTH_KB / step_kb) else 0L
cat(sprintf("  Window step: %.1f kb | Cov smooth: %d kb | Freq smooth: %d kb\n",
            step_kb, COV_SMOOTH_KB, FREQ_SMOOTH_KB))

# ── Step 2: Extract per-window data ──────────────────────────────────────────
cat("Step 2: Extracting window data...\n")
window_data <- map2(windows$data, windows$CHROM,
                    extract_window_data, Nfounders = Nfounders)
valid <- !sapply(window_data, is.null)
cat(sprintf("  %d valid windows (%.1f%%)\n", sum(valid), 100 * mean(valid)))

# ── Step 3: Pool across replicates ───────────────────────────────────────────
cat("Step 3: Pooling replicates...\n")
pooled <- map(window_data, pool_one_window)

# ── Step 4: Stack pooled covariances into 3D arrays ─────────────────────────
cat(sprintf("Step 4: Smoothing covariance matrices (+/-%d kb)...\n",
            COV_SMOOTH_KB))
nF <- Nfounders

covar1_stack <- array(NA_real_, c(nF, nF, W))
covar2_stack <- array(NA_real_, c(nF, nF, W))
p1_vec       <- matrix(NA_real_, W, nF)
p2_vec       <- matrix(NA_real_, W, nF)

for (w in which(valid)) {
  covar1_stack[,,w] <- pooled[[w]]$covar1
  covar2_stack[,,w] <- pooled[[w]]$covar2
  p1_vec[w,]        <- pooled[[w]]$p1
  p2_vec[w,]        <- pooled[[w]]$p2
}

# Smooth each (i,j) element of the covariance matrices across windows
covar1_smooth <- array(NA_real_, c(nF, nF, W))
covar2_smooth <- array(NA_real_, c(nF, nF, W))

for (i in 1:nF) {
  for (j in 1:nF) {
    covar1_smooth[i,j,] <- running_mean(covar1_stack[i,j,], COV_HALF)
    covar2_smooth[i,j,] <- running_mean(covar2_stack[i,j,], COV_HALF)
  }
}

# ── Step 4b: Smooth haplotype frequency vectors (optional) ───────────────────
# A running mean of simplex vectors is simplex-preserving (convex combination).
# FREQ_SMOOTH_KB=0 skips this step entirely.
# When active, also fills gaps at originally-invalid windows: after smoothing,
# those windows have valid interpolated frequencies from neighbors, so we extend
# the Wald computation to all windows with non-NA smoothed data.
if (FREQ_HALF > 0L) {
  cat(sprintf("Step 4b: Smoothing frequency vectors (+/-%d kb)...\n", FREQ_SMOOTH_KB))
  for (f in seq_len(nF)) {
    p1_vec[, f] <- running_mean(p1_vec[, f], FREQ_HALF)
    p2_vec[, f] <- running_mean(p2_vec[, f], FREQ_HALF)
  }
  # Any window with non-NA smoothed frequencies and covariances can get a Wald score
  smoothed_valid <- !is.na(p1_vec[, 1]) & !is.na(p2_vec[, 1]) &
                    !is.na(covar1_smooth[1, 1, ]) & !is.na(covar2_smooth[1, 1, ])
  cat(sprintf("  Wald-eligible windows after freq smoothing: %d (was %d valid raw)\n",
              sum(smoothed_valid), sum(valid)))
} else {
  smoothed_valid <- valid
}

# ── Step 5: Wald test with smoothed covariances ─────────────────────────────
cat("Step 5: Computing stabilized Wald tests...\n")
Wald_log10p <- rep(NA_real_, W)
avg_var     <- rep(NA_real_, W)

for (w in which(smoothed_valid)) {
  wt <- wald_from_pooled(p1_vec[w,], p2_vec[w,],
                          covar1_smooth[,,w], covar2_smooth[,,w])
  Wald_log10p[w] <- -log10(wt$p.value)
  avg_var[w]     <- wt$avg.var
}

# ── Step 6: Pseu_log10p and Heritability from RAW data ───────────────────────
cat("Step 6: Computing pseudo-N test and heritability (raw data)...\n")
Pseu_log10p <- rep(NA_real_, W)
Falc_H2     <- rep(NA_real_, W)
Cutl_H2     <- rep(NA_real_, W)

for (w in which(valid)) {
  wd <- window_data[[w]]
  Pseu_log10p[w] <- pseudoN.test(wd$p1, wd$p2, wd$covar1, wd$covar2,
                                  wd$nrepl, wd$N1, wd$N2)
  af_cutoff <- 0.01
  temp <- Heritability(wd$p1, wd$p2, wd$nrepl, ProportionSelect, af_cutoff)
  Falc_H2[w] <- temp$Falconer_H2
  Cutl_H2[w] <- temp$Cutler_H2
}

# ── Step 7: Assemble scan output ─────────────────────────────────────────────
cat("Step 7: Assembling output...\n")
bb2 <- tibble(
  chr         = windows$CHROM,
  pos         = windows$pos,
  Wald_log10p = Wald_log10p,
  Pseu_log10p = Pseu_log10p,
  Falc_H2     = Falc_H2,
  Cutl_H2     = Cutl_H2,
  avg.var     = avg_var
)
bb3 <- add_genetic(bb2)

# ── Step 8: meansBySample (identical to current, uses raw frequencies) ───────
cat("Step 8: Computing meansBySample...\n")
keep <- valid & !is.na(Pseu_log10p)
valid_positions <- windows[keep, ] %>%
  ungroup() %>%
  select(CHROM, pos)

bb4 <- valid_positions %>%
  left_join(xx1, by = c("CHROM", "pos")) %>%
  select(-c(Err, Groups)) %>%
  unnest(c(sample, Haps, Names)) %>%
  unnest(c(Haps, Names)) %>%
  rename(chr = CHROM, pool = sample, freq = Haps, founder = Names) %>%
  left_join(design.df, by = c("pool" = "bam")) %>%
  select(c(chr, pos, TRT, REP, REPrep, freq, founder)) %>%
  filter(!is.na(TRT)) %>%
  group_by(chr, pos, TRT, REP, founder) %>%
  summarize(freq = mean(freq, na.rm = TRUE))

write.table(bb3, fileout)
write.table(bb4, fileout_meansBySample)
cat("Done.\n")
