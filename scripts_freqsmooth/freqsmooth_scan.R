#!/usr/bin/env Rscript
###############################################################################
# freqsmooth_scan.R  —  Step 2 of the freqsmooth pipeline
#
# Reads the smoothed RDS from smooth_haps.R. For each window, reconstructs the
# per-replicate frequency matrices and covariance arrays from the long-format
# smoothed data, pools replicates with effective-N weighting, and runs the
# stabilized Wald test (eigenvalue-regularized). Heritability is computed from
# the same smoothed frequencies as the Wald test.
#
# Usage:
#   Rscript scripts_freqsmooth/freqsmooth_scan.R \
#       --chr chrX --dir process/proj/SCAN_NAME \
#       --outdir SCAN_NAME --rfile helpfiles/proj/design.txt
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

RIDGE_FRACTION <- 0.01
AF_CUTOFF      <- 0.01

# ── Arguments ─────────────────────────────────────────────────────────────────
args   <- commandArgs(trailingOnly = TRUE)
parsed <- list()
i <- 1L
while (i <= length(args)) {
  switch(args[i],
    "--chr"    = { parsed$chr    <- args[i+1]; i <- i+2L },
    "--dir"    = { parsed$dir    <- args[i+1]; i <- i+2L },
    "--outdir" = { parsed$outdir <- args[i+1]; i <- i+2L },
    "--rfile"  = { parsed$rfile  <- args[i+1]; i <- i+2L },
    stop(paste("Unknown argument:", args[i]))
  )
}

mychr <- parsed$chr
design.df <- read.table(parsed$rfile, header = TRUE)
source("scripts/scan_functions.R")

filein  <- file.path(parsed$dir, paste0(parsed$outdir, ".smooth.", mychr, ".rds"))
fileout <- file.path(parsed$dir, paste0(parsed$outdir, ".scan.", mychr, ".txt"))

options(dplyr.summarise.inform = FALSE)

# ── Load smoothed data ────────────────────────────────────────────────────────
cat("Reading", filein, "\n")
sm            <- readRDS(filein)
freq_smoothed <- sm$freq       # (CHROM, pos, TRT, REP, founder, freq, Num)
err_smoothed  <- sm$err        # (CHROM, pos, TRT, REP, fi, fj, v)
founder_names <- sm$founder_names
nrepl         <- sm$nrepl
nF            <- length(founder_names)

ProportionSelect <- design.df %>%
  filter(TRT == "Z") %>% select(REP, Proportion) %>% arrange(REP)

# ── Pre-build arrays once — avoids repeated dplyr ops inside window loop ──────
cat(sprintf("Building arrays for %s...\n", mychr))

freq_chr <- freq_smoothed %>% filter(CHROM == mychr)
err_chr  <- err_smoothed  %>% filter(CHROM == mychr)

win_pos <- sort(unique(freq_chr$pos))
W       <- length(win_pos)

# Sample sizes: constant per REP across windows
N1 <- freq_chr %>% filter(TRT == "C") %>%
  distinct(REP, Num) %>% arrange(REP) %>% pull(Num)
N2 <- freq_chr %>% filter(TRT == "Z") %>%
  distinct(REP, Num) %>% arrange(REP) %>% pull(Num)

# Frequency arrays: W x nrepl x nF
# pivot_wider once for the full chromosome then reshape
p_wide <- freq_chr %>%
  select(pos, TRT, REP, founder, freq) %>%
  pivot_wider(names_from = founder, values_from = freq)

to_3d <- function(trt) {
  m <- p_wide %>% filter(TRT == trt) %>%
    arrange(pos, REP) %>%
    select(all_of(founder_names)) %>%
    as.matrix()                             # (W*nrepl) x nF, pos-major order
  aperm(array(t(m), c(nF, nrepl, W)), c(3, 2, 1))  # W x nrepl x nF
}
P1 <- to_3d("C")
P2 <- to_3d("Z")

# Covariance arrays: W x nF x nF x nrepl
# array() fills fi-major; sort order must be (pos, REP, fj, fi)
to_4d <- function(trt) {
  v <- err_chr %>% filter(TRT == trt) %>%
    arrange(pos, REP, fj, fi) %>% pull(v)
  aperm(array(v, c(nF, nF, nrepl, W)), c(4, 1, 2, 3))  # W x fi x fj x nrepl
}
E1 <- to_4d("C")
E2 <- to_4d("Z")

# ── Per-window scan: array indexing only, no dplyr inside loop ───────────────
cat(sprintf("Running Wald scan on %s (%d windows)...\n", mychr, W))
df <- nF - 1L

scan_list <- lapply(seq_len(W), function(w) {
  p1     <- matrix(P1[w,,], nrepl, nF)
  p2     <- matrix(P2[w,,], nrepl, nF)
  covar1 <- array(E1[w,,,], c(nF, nF, nrepl))
  covar2 <- array(E2[w,,,], c(nF, nF, nrepl))
  if (anyNA(p1) || anyNA(p2) || anyNA(covar1) || anyNA(covar2)) return(NULL)

  # Pool replicates
  N1e <- N2e <- numeric(nrepl)
  cv1 <- cv2 <- array(0, c(nF, nF, nrepl))
  for (r in seq_len(nrepl)) {
    pm  <- (N1[r]*p1[r,] + N2[r]*p2[r,]) / (N1[r]+N2[r])
    cm1 <- mn.covmat(pm, 2*N1[r]); cm2 <- mn.covmat(pm, 2*N2[r])
    N1e[r] <- sum(diag(cm1))*4*N1[r]^2 /
      (sum(diag(cm1))*2*N1[r] + 2*N1[r]*sum(diag(covar1[,,r])))
    N2e[r] <- sum(diag(cm2))*4*N2[r]^2 /
      (sum(diag(cm2))*2*N2[r] + 2*N2[r]*sum(diag(covar2[,,r])))
    cv1[,,r] <- (cm1 + covar1[,,r]) * N1e[r]^2
    cv2[,,r] <- (cm2 + covar2[,,r]) * N2e[r]^2
  }
  p1p <- as.vector(N1e %*% p1 / sum(N1e))
  p2p <- as.vector(N2e %*% p2 / sum(N2e))
  c1p <- rowSums(cv1, dims=2) / sum(N1e)^2
  c2p <- rowSums(cv2, dims=2) / sum(N2e)^2

  # Stabilized Wald test
  covar <- c1p + c2p
  eg    <- eigen(covar)
  eval  <- pmax(eg$values[1:df], RIDGE_FRACTION * mean(eg$values[1:df]))
  trafo <- diag(1/sqrt(eval)) %*% t(eg$vectors[,1:df])
  tstat <- sum((trafo %*% (p1p - p2p))^2)

  # Heritability
  h2 <- tryCatch(
    Heritability(p1, p2, nrepl, ProportionSelect, AF_CUTOFF),
    error = function(e) list(Falconer_H2 = NA_real_, Cutler_H2 = NA_real_)
  )

  list(chr         = mychr,
       pos         = win_pos[w],
       Wald_log10p = -log10(pchisq(tstat, df, lower.tail = FALSE)),
       Falc_H2     = h2$Falconer_H2,
       Cutl_H2     = h2$Cutler_H2)
})

scan_results <- bind_rows(Filter(Negate(is.null), scan_list))

# ── Add genetic position and write ───────────────────────────────────────────
cat("Writing", fileout, "\n")
scan_results %>% add_genetic() %>% write.table(fileout)
cat("Done.\n")
