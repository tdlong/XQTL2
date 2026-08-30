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

AF_CUTOFF <- 0.01     # matches hap_scan.R

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
  # Only the diagonal of the reconstruction covariance enters the Falconer
  # bias; fi indexes founder_names, the same order hap_scan.R assumes.
  errd <- sm$err %>% filter(fi == fj) %>%
    transmute(CHROM, pos, TRT, REP = as.character(REP),
              founder = founder_names[fi], vrec = v)

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
           vtot = mult + vrec) %>%
    select(CHROM, pos, REP, founder, TRT, freq, vtot, mult, vrec) %>%
    pivot_wider(names_from = TRT, values_from = c(freq, vtot, mult, vrec))

  need <- c("freq_C", "freq_Z", "vtot_C", "vtot_Z")
  if (!all(need %in% names(d)))
    stop(basename(f), ": need both TRT arms, have ",
         paste(setdiff(need, names(d)), collapse = ", "), " missing")

  d %>%
    inner_join(ProportionSelect, by = "REP") %>%     # drops REPs with no Z row
    filter(!is.na(Proportion)) %>%
    mutate(denom = pmax(freq_C, AF_CUTOFF)) %>%
    group_by(CHROM, pos, REP, Proportion) %>%
    summarize(raw_sum  = sum((freq_Z - freq_C)^2 / denom),
              bias_sum = sum((vtot_Z + vtot_C)   / denom),
              mult_sum = sum((mult_Z + mult_C)   / denom),
              rec_sum  = sum((vrec_Z + vrec_C)   / denom),
              n_floor  = sum(freq_C <= AF_CUTOFF),
              .groups  = "drop") %>%
    mutate(Falcon_i = dnorm(qnorm(1 - Proportion)) / Proportion,
           h2_raw   = 200 * raw_sum  / Falcon_i^2,
           h2_bias  = 200 * bias_sum / Falcon_i^2,
           h2_mult  = 200 * mult_sum / Falcon_i^2,
           h2_rec   = 200 * rec_sum  / Falcon_i^2) %>%
    # the pipeline's h2 is the mean over replicates
    group_by(CHROM, pos) %>%
    summarize(nrepl    = n(),
              h2_raw   = mean(h2_raw),
              h2_bias  = mean(h2_bias),
              h2_mult  = mean(h2_mult),
              h2_rec   = mean(h2_rec),
              n_floor  = mean(n_floor),
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
