###############################################################################
# h2_from_scan.R  --  bias-corrected Falconer h2 from a smoothed scan
#
# Recomputes per-window heritability from the smoothed RDS that smooth_haps.R
# wrote, without re-running the Wald scan.  Reports the raw estimator, the
# exact squaring bias, and the corrected value.
#
#   h2_raw  = 200/i^2 * sum_f (Z_f - C_f)^2                  / max(C_f, cutoff)
#   h2_bias = 200/i^2 * sum_f (vZ_f + vC_f)                  / max(C_f, cutoff)
#   h2_corr = 200/i^2 * sum_f ((Z_f-C_f)^2 - vZ_f - vC_f)    / max(C_f, cutoff)
#
# where i = dnorm(qnorm(1-P))/P and v is the variance of each arm's founder
# frequency estimate: multinomial sampling of that arm's own flies, plus the
# lsei reconstruction error.  Both are needed -- the multinomial term alone
# understates the bias, which is why this reads the RDS and not the
# meansBySample file (that file carries neither Num nor the reconstruction
# covariance).
#
# E[x^2] = x_true^2 + Var(x), so h2_raw is inflated by h2_bias whether or not
# there is any real signal.  The inflation scales with the variance, so it is
# larger wherever the frequencies are less well determined -- notably chrX in
# male pools, which carry half the chromosomes per fly (XQTL2 #40).  Comparing
# raw h2 between chrX and an autosome, or between sexes on chrX, compares two
# different offsets.  h2_corr is the comparable quantity.
#
# h2_corr is unbiased, so it scatters negative where the true h2 is zero; that
# is correct behaviour for an unbiased estimator of a non-negative quantity and
# not a defect.  h2_corr_pos is its positive part, for when a non-negative
# track is wanted.
#
# The chromosome and the pool sex do NOT have to be supplied: CHROM is in the
# data, Num in the RDS was already scaled for chrX dosage by smooth_haps.R, and
# the RDS records the --sex that produced it.  The design file is read for one
# thing only, Proportion, which sets the Falconer selection intensity.
#
# Usage:
#   Rscript scripts/h2_from_scan.R \
#       --dir process/proj --scan SCAN_NAME \
#       --rfile helpfiles/proj/design.txt \
#       [--chr chrX] [--out FILE]
#
# With no --chr, every <scan>.smooth.<chr>.rds in the scan directory is used.
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

# fill_gaps()/running_mean() are shared with smooth_haps.R.
script_dir <- dirname(normalizePath(sub("--file=", "",
                grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "scan_functions.R"))

AF_CUTOFF   <- 0.01   # matches hap_scan.R
CONFOUND_R  <- -0.99  # correlation below which a founder pair is unidentified
FILL_HALF   <- 50L    # flanking windows averaged for the gap-fill anchors


# ── Per-window variance component ────────────────────────────────────────────
# h2 is a variance, so estimate it as one.  On the scale u_f = (Z_f-C_f)/sqrt(C_f)
# the statistic is h2 = 200/i^2 * sum_f u_f^2, and each observed u_hat_f carries
# known noise s_f = (vZ_f+vC_f)/C_f.  Model u_f ~ N(0, tau2) and fit tau2 by ML
# from the nF founders AT THAT WINDOW:
#
#   u_hat_f ~ N(0, tau2 + s_f)
#
# tau2 is fit per window, not genome-wide.  That matters: with a single global
# prior, the ~95% of windows that are null drag tau2 to nothing and real QTL get
# shrunk by 5x.  Fitting locally has no such failure.
#
# Two properties this buys over subtracting the bias:
#   - Non-negativity is a boundary solution (tau2 = 0 is where the likelihood
#     actually peaks at a null window), not a clamp applied afterwards.
#   - Founders are weighted by 1/(tau2+s_f), so a badly determined founder
#     contributes less.  Plain subtraction weights every founder equally.
#
# Fisher-scoring fixed point, vectorised over windows.  At convergence
# sum(w^2 u^2) = sum(w), the ML score equation, with w = 1/(tau2+s).
fit_tau2 <- function(U, S, iter = 50L) {
  t2 <- pmax(rowMeans(U^2 - S), 0)                  # method-of-moments start
  for (k in seq_len(iter)) {
    w  <- 1 / (t2 + S)
    t2 <- pmax(rowSums(w^2 * (U^2 - S)) / rowSums(w^2), 0)
  }
  t2
}

# Posterior E[sum u^2 | u_hat] under the fitted tau2: shrink each founder by
# k_f = tau2/(tau2+s_f), then add back the posterior variance k_f*s_f.
vc_sum <- function(U, S, t2) {
  k <- t2 / (t2 + S)
  rowSums(k^2 * U^2 + k * S)
}

# ── Arguments ─────────────────────────────────────────────────────────────────
args   <- commandArgs(trailingOnly = TRUE)
parsed <- list()
i <- 1L
while (i <= length(args)) {
  switch(args[i],
    "--dir"   = { parsed$dir   <- args[i+1]; i <- i+2L },
    "--scan"  = { parsed$scan  <- args[i+1]; i <- i+2L },
    "--rfile" = { parsed$rfile <- args[i+1]; i <- i+2L },
    "--chr"   = { parsed$chr   <- args[i+1]; i <- i+2L },
    "--out"   = { parsed$out   <- args[i+1]; i <- i+2L },
    "--check" = { parsed$check <- TRUE;          i <- i+1L },
    stop(paste("Unknown argument:", args[i]))
  )
}
for (need in c("dir", "scan", "rfile"))
  if (is.null(parsed[[need]])) stop("missing required argument: --", need)

dirin <- file.path(parsed$dir, "Scans", parsed$scan)
if (!dir.exists(dirin)) stop("no such scan directory: ", dirin)

# ── R2 smoothing calibration ──────────────────────────────────────────────────
# mn.covmat() + covar describe a SINGLE window's estimate, but the frequencies
# being squared are smoothed over +/-smooth_kb, so their variance is smaller.
# hap_scan.R already corrects for this on the Wald side: it evaluates
# pchisq(tstat / R2_SMOOTH), which is the statistic recomputed against a
# covariance of R2 * Sigma_single.  So Var(smoothed) = R2 * Var(single), and
# the same factor belongs on the h2 bias.  Without it the bias over-subtracts
# by 1/R2 and drives the floor negative (XQTL2 #40).
#
# Same file and same default as hap_scan.R: absent means no correction.
r2_file <- file.path(dirin, paste0(parsed$scan, ".smooth_r2.txt"))
R2_SMOOTH <- if (file.exists(r2_file)) {
  r2 <- as.numeric(readLines(r2_file))
  cat(sprintf("Smoothing R2: %.4f (variance of the smoothed estimate is %.4f of a single window's)\n", r2, r2))
  r2
} else {
  cat("No ", basename(r2_file), " -- no smoothing correction (R2 = 1).\n",
      "  The bias will over-subtract; run smooth_r2_diag.sh to generate it.\n", sep = "")
  1.0
}

design.df <- read.table(parsed$rfile, header = TRUE)
if (!"Proportion" %in% names(design.df))
  stop("design file has no Proportion column: ", parsed$rfile)

# Proportion is a property of the selected arm, keyed by replicate.
# REP is a LABEL, not an index (XQTL2 #32), and read.table may type it
# differently than the RDS does.  Compare as character on both sides.
ProportionSelect <- design.df %>%
  filter(TRT == "Z") %>% select(REP, Proportion) %>% distinct() %>%
  mutate(REP = as.character(REP))

rds_files <- if (!is.null(parsed$chr)) {
  f <- file.path(dirin, paste0(parsed$scan, ".smooth.", parsed$chr, ".rds"))
  if (!file.exists(f)) stop("no such file: ", f)
  f
} else {
  f <- list.files(dirin, pattern = paste0("^", parsed$scan, "\\.smooth\\..*\\.rds$"),
                  full.names = TRUE)
  if (!length(f)) stop("no ", parsed$scan, ".smooth.*.rds in ", dirin)
  sort(f)
}

fileout <- if (!is.null(parsed$out)) parsed$out else
  file.path(dirin, paste0(parsed$scan, ".h2_falconer.txt"))

# ── Per-chromosome ────────────────────────────────────────────────────────────
one_chr <- function(f) {
  sm <- readRDS(f)
  founder_names <- sm$founder_names
  nF <- length(founder_names)

  # sex/xfactor are recorded only by smooth_haps.R from d721210 onward.
  sex     <- if (!is.null(sm$sex))     sm$sex     else NA_character_
  xfactor <- if (!is.null(sm$xfactor)) sm$xfactor else NA_real_

  freq <- sm$freq

  # Repair the covariance for unresolvable founders (XQTL2 #40).
  #
  # Where >1 founder shares a group, lsei identifies only the group sum, not the
  # individual founders.  It returns a perfectly anti-correlated block: huge,
  # equal marginal variances that cancel in their sum.  smooth_haps.R masks and
  # gap-fills the FREQUENCIES of those founders but, before this fix, left the
  # covariance untouched -- so the bias subtracted error bars up to 48,000x too
  # large for frequencies that had already been repaired.
  #
  # smooth_haps.R now masks the covariance at source.  Files written before that
  # still carry the raw block, so detect and repair it here too, which lets an
  # existing scan be corrected without re-running the smooth step.  A confounded
  # pair is identified by its correlation: exactly -1 up to rounding, against a
  # bulk that sits near 0.  On real chrX data this flags 6.2% of founder-windows
  # and they carry 99.3% of the total variance mass.
  vd <- sm$err %>% filter(fi == fj) %>% select(CHROM, pos, TRT, REP, fi, vdiag = v)
  confounded <- sm$err %>% filter(fi != fj) %>%
    left_join(vd,                     by = c("CHROM","pos","TRT","REP","fi")) %>%
    left_join(vd %>% rename(fj = fi, vdiag_j = vdiag),
                                      by = c("CHROM","pos","TRT","REP","fj")) %>%
    mutate(r = v / sqrt(vdiag * vdiag_j)) %>%
    group_by(CHROM, pos, TRT, REP, fi) %>%
    summarize(bad = any(r < CONFOUND_R, na.rm = TRUE), .groups = "drop")

  errd <- vd %>%
    left_join(confounded, by = c("CHROM","pos","TRT","REP","fi")) %>%
    mutate(vrec = if_else(coalesce(bad, FALSE), NA_real_, vdiag)) %>%
    arrange(CHROM, pos) %>%
    group_by(CHROM, TRT, REP, fi) %>%
    # gap-fill from flanking windows, as the frequencies already were
    mutate(vrec = fill_gaps(vrec, FILL_HALF)) %>%
    ungroup() %>%
    transmute(CHROM, pos, TRT, REP = as.character(REP),
              founder = founder_names[fi], vrec, repaired = coalesce(bad, FALSE))

  n_rep <- sum(errd$repaired)
  if (n_rep)
    cat(sprintf("  repaired %d / %d founder-windows (%.1f%%) with unresolvable covariance\n",
                n_rep, nrow(errd), 100 * n_rep / nrow(errd)))

  freq <- freq %>% mutate(REP = as.character(REP))

  d <- freq %>%
    inner_join(errd, by = c("CHROM", "pos", "TRT", "REP", "founder")) %>%
    # Multinomial term, matching mn.covmat() exactly.  That function rescales
    # its p to sum to 1 before taking p(1-p)/n, so do the same here: smoothing
    # applies pmax(freq, 3e-4) and gap-filling after the group rescaling, and
    # neither preserves the sum, so the founder frequencies at a window need
    # not add to 1.  Using the raw freq here would silently disagree with
    # hap_scan.R by the square of that discrepancy.  The Falconer denominator
    # below uses the UNnormalised freq, which is what Heritability() does.
    group_by(CHROM, pos, TRT, REP) %>%
    mutate(pnorm_ = freq / sum(freq)) %>%
    ungroup() %>%
    # Num is already chrX-scaled; 2*Num is the chromosome count.
    mutate(mult = pnorm_ * (1 - pnorm_) / (2 * Num),
           # R2 converts a single-window variance to the smoothed estimate's
           vtot = R2_SMOOTH * (mult + vrec)) %>%
    select(CHROM, pos, REP, founder, TRT, freq, vtot, mult, vrec, repaired) %>%
    pivot_wider(names_from = TRT, values_from = c(freq, vtot, mult, vrec, repaired))

  need <- c("freq_C", "freq_Z", "vtot_C", "vtot_Z")
  if (!all(need %in% names(d)))
    stop(basename(f), ": need both TRT arms, have ",
         paste(setdiff(need, names(d)), collapse = ", "), " missing")

  d %>%
    inner_join(ProportionSelect, by = "REP") %>%     # drops REPs with no Z row
    filter(!is.na(Proportion)) %>%
    mutate(denom = pmax(freq_C, AF_CUTOFF),
           u     = (freq_Z - freq_C) / sqrt(denom),
           sv    = (vtot_Z + vtot_C) / denom) %>%
    arrange(CHROM, pos, REP, founder) %>%
    group_by(CHROM, pos, REP, Proportion) %>%
    summarize(raw_sum  = sum(u^2),
              bias_sum = sum(sv),
              mult_sum = sum((mult_Z + mult_C)   / denom),
              rec_sum  = sum((vrec_Z + vrec_C)   / denom),
              n_floor  = sum(freq_C <= AF_CUTOFF),
              n_repair = sum(repaired_C | repaired_Z),
              nf       = n(),
              u_v      = list(u), s_v = list(sv),
              .groups  = "drop") %>%
    # variance component, one tau2 per (window, replicate)
    { dd <- .
      keep <- dd$nf == nF
      vcs  <- rep(NA_real_, nrow(dd))
      if (any(keep)) {
        U  <- matrix(unlist(dd$u_v[keep]), ncol = nF, byrow = TRUE)
        S  <- matrix(unlist(dd$s_v[keep]), ncol = nF, byrow = TRUE)
        t2 <- fit_tau2(U, S)
        vcs[keep] <- vc_sum(U, S, t2)
      }
      dd$vc_sum <- vcs
      dd %>% select(-u_v, -s_v, -nf)
    } %>%
    mutate(Falcon_i = dnorm(qnorm(1 - Proportion)) / Proportion,
           h2_raw   = 200 * raw_sum  / Falcon_i^2,
           h2_bias  = 200 * bias_sum / Falcon_i^2,
           h2_mult  = 200 * mult_sum / Falcon_i^2,
           h2_rec   = 200 * rec_sum  / Falcon_i^2,
           h2_vc    = 200 * vc_sum   / Falcon_i^2) %>%
    # the pipeline's h2 is the mean over replicates
    group_by(CHROM, pos) %>%
    summarize(nrepl    = n(),
              h2_raw   = mean(h2_raw),
              h2_bias  = mean(h2_bias),
              h2_mult  = mean(h2_mult),
              h2_rec   = mean(h2_rec),
              h2_vc    = mean(h2_vc),
              n_floor  = mean(n_floor),
              n_repair = mean(n_repair),
              .groups  = "drop") %>%
    mutate(h2_corr     = h2_raw - h2_bias,
           h2_corr_pos = pmax(0, h2_corr),
           sex         = sex,
           xfactor     = xfactor)
}

res <- bind_rows(lapply(rds_files, function(f) {
  cat("Reading", basename(f), "\n")
  one_chr(f)
}))

cat("\n")
# Which half of the variance dominates the bias, and how many founders sit at
# the AF_CUTOFF floor -- those enter the bias divided by 0.01, so a handful of
# badly determined low-frequency founders can carry the whole term.
cat("Bias breakdown (medians): mult = multinomial half, rec = lsei reconstruction half\n")
res %>% group_by(CHROM) %>%
  summarize(h2_raw = median(h2_raw), h2_bias = median(h2_bias),
            from_mult = median(h2_mult), from_recon = median(h2_rec),
            founders_at_floor = median(n_floor), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n")

# Cross-check against the pipeline's own columns.  h2_raw must reproduce
# Falc_H2 and h2_bias must reproduce Falc_H2_bias; if they do, any oddity in
# h2_corr is in the pipeline's bias, not in this script.
if (isTRUE(parsed$check)) {
  sf <- list.files(dirin, pattern = paste0("^", parsed$scan, "\\.scan\\..*\\.txt$"),
                   full.names = TRUE)
  if (!length(sf)) {
    cat("--check: no", paste0(parsed$scan, ".scan.<chr>.txt"), "in", dirin, "\n\n")
  } else {
    sc <- bind_rows(lapply(sf, function(x) read.table(x, header = TRUE)))
    if (!all(c("chr","pos","Falc_H2","Falc_H2_bias") %in% names(sc))) {
      cat("--check: scan file lacks Falc_H2/Falc_H2_bias; skipping\n\n")
    } else {
      cmp <- res %>% inner_join(sc %>% select(chr, pos, Falc_H2, Falc_H2_bias),
                                by = c("CHROM" = "chr", "pos" = "pos"))
      cat(sprintf("--check: %d windows matched against %s\n", nrow(cmp), basename(sf[1])))
      if (nrow(cmp)) {
        cat(sprintf("  max |h2_raw  - Falc_H2|      = %.3e\n",
                    max(abs(cmp$h2_raw  - cmp$Falc_H2),      na.rm = TRUE)))
        cat(sprintf("  max |h2_bias - Falc_H2_bias| = %.3e\n",
                    max(abs(cmp$h2_bias - cmp$Falc_H2_bias), na.rm = TRUE)))
        cat("  (both ~0 means this script agrees with the pipeline exactly,\n",
            "   so any large negative h2_corr is the pipeline's bias term.)\n", sep = "")
      }
      cat("\n")
    }
  }
}

res %>% group_by(CHROM, sex, xfactor) %>%
  summarize(windows = n(),
            raw     = median(h2_raw),
            bias    = median(h2_bias),
            corr    = median(h2_corr),
            .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nWriting", fileout, "\n")
write.table(res, fileout)
cat("Done.\n")
