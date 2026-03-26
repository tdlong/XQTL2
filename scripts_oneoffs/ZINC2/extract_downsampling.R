#!/usr/bin/env Rscript
###############################################################################
# extract_downsampling.R  (server-side, self-contained)
#
# Concatenates all downsampled pseudoscan files into a single compact TSV.
# Adds columns: sex (F/M) and N (7/9/11/13).
# N=15 is covered by the main ZINC2_F_freqs250 / ZINC2_M_freqs250 scans.
#
# Run from /dfs7/adl/tdlong/fly_pool/XQTL2:
#   Rscript scripts/extract_downsampling.R
#
# Output: process/ZINC2_v2/downsample_pseudoscans.tsv.gz
#   ~200k rows, ~3-5 MB gzipped — download this instead of the full tarball.
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

BASE <- "process/ZINC2_v2"

scans <- expand.grid(
  sex = c("F", "M"),
  N   = c(7, 9, 11, 13),
  stringsAsFactors = FALSE
)

dat <- lapply(seq_len(nrow(scans)), function(i) {
  sex  <- scans$sex[i]
  N    <- scans$N[i]
  Npad <- sprintf("%02d", N)
  name <- sprintf("ZINC2_%s_freqs250_N%s", sex, Npad)
  path <- file.path(BASE, name, paste0(name, ".pseudoscan.txt"))

  if (!file.exists(path)) {
    warning("Missing: ", path)
    return(NULL)
  }

  read.table(path, header = TRUE) %>%
    as_tibble() %>%
    mutate(sex = sex, N = N) %>%
    select(sex, N, chr, pos, Wald_log10p, Pseu_log10p, Falc_H2, Cutl_H2, avg.var, cM)
}) %>%
  bind_rows()

outfile <- file.path(BASE, "downsample_pseudoscans.tsv.gz")
write_tsv(dat, outfile)
cat(sprintf("Wrote %d rows to %s\n", nrow(dat), outfile))
