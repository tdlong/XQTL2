#!/usr/bin/env Rscript
###############################################################################
# catalog2snptable.R — build the SNP-scan founder table from the SNP catalog
#
# The SNP scan (snp_scan.R) imputes each pool's ALT frequency at a SNP as
#
#     f_ALT(pool, SNP) = h(pool, window) . s(SNP)
#
# where h is the smoothed founder haplotype frequency vector and s is the vector
# of per-founder ALT frequencies at that SNP.  h is derived from RefAlt, which is
# counted against the catalog — so s must come from the same catalog, or the two
# halves of that product describe different ascertainments of the same founders.
#
# This script produces s directly from the catalog:
#
#   positions  : catalog.tsv.gz  (the filtered catalog — exactly the sites the
#                samples were counted at, under whatever thresholds you cut with)
#   founder ALT: catalog.annot.tsv.gz AD_<founder> columns, alt/(ref+alt)
#   founders   : catalog.founders.txt, in catalog column order
#   cM         : helpfiles/flymap.r6.txt, same interpolation as the scans use
#
# Output has the schema snp_scan.R expects:
#
#   CHROM  POS  <founder1> ... <founderN>  cM
#
# Rebuild this whenever you rebuild or re-cut the catalog. A re-cut changes which
# positions are in catalog.tsv.gz, so the SNP table changes with it.
#
# Usage:
#   Rscript scripts/catalog2snptable.R \
#       --catdir process/<project>/Catalog \
#       --out    helpfiles/<pop>_SNPs.cM.txt.gz
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

args   <- commandArgs(trailingOnly = TRUE)
parsed <- list()
i <- 1L
while (i <= length(args)) {
  switch(args[i],
    "--catdir" = { parsed$catdir <- args[i+1]; i <- i+2L },
    "--out"    = { parsed$out    <- args[i+1]; i <- i+2L },
    stop(paste("Unknown argument:", args[i]))
  )
}
for (need in c("catdir","out"))
  if (is.null(parsed[[need]])) stop("Missing required argument: --", need)

script_dir <- dirname(normalizePath(sub("--file=", "",
                grep("--file=", commandArgs(FALSE), value = TRUE))))

annot_f <- file.path(parsed$catdir, "catalog.annot.tsv.gz")
filt_f  <- file.path(parsed$catdir, "catalog.tsv.gz")
fnd_f   <- file.path(parsed$catdir, "catalog.founders.txt")
for (f in c(annot_f, filt_f, fnd_f))
  if (!file.exists(f)) stop("Not found: ", f,
    "\nRun build_catalog.sh (and catalog_filter.sh) first — see README, 'Build the catalog'.")

zread <- function(f, ...) fread(cmd = paste("gzip -dc", shQuote(f)), ...)

founders <- readLines(fnd_f)
founders <- founders[nzchar(founders)]
cat(sprintf("Founders (%d, catalog order): %s\n",
            length(founders), paste(founders, collapse = ", ")))

# ── Filtered catalog: the positions the samples were actually counted at ──────
cat("Reading", filt_f, "...\n")
keep <- zread(filt_f, header = FALSE, select = 1:2,
              col.names = c("CHROM", "POS"))
cat(sprintf("  %d filtered-catalog positions\n", nrow(keep)))

# ── Annotated catalog: per-founder AD at every candidate site ────────────────
cat("Reading", annot_f, "...\n")
ad_cols <- paste0("AD_", founders)
# AD is "ref,alt". fread will otherwise read that comma as a DECIMAL separator --
# "19,0" becomes 19.0 and "0,33" becomes 0.33 -- silently destroying every count.
# Force character and parse the pair ourselves.
annot   <- zread(annot_f, select = c("CHROM", "POS", ad_cols),
                 colClasses = setNames(rep("character", length(ad_cols)), ad_cols))
cat(sprintf("  %d candidate sites annotated\n", nrow(annot)))

setkey(annot, CHROM, POS)
setkey(keep,  CHROM, POS)
snp <- annot[keep, nomatch = 0L]
cat(sprintf("  %d sites carried through to the SNP table\n", nrow(snp)))
if (nrow(snp) == 0L)
  stop("No filtered-catalog position matched the annotated catalog — are these from the same build?")

# ── AD "ref,alt" -> ALT frequency ────────────────────────────────────────────
# bcftools writes "." for a site with no reads in that founder; the catalog's
# --min-dp filter should have removed those already, so treat any that survive
# as missing rather than silently calling them REF.
cat("Converting AD to per-founder ALT frequency...\n")
for (j in seq_along(founders)) {
  v   <- as.character(snp[[ad_cols[j]]])
  # A well-formed AD is "ref,alt". Anything else -- "." for no call, or a
  # leftover multiallelic "ref,alt1,alt2" -- must not be read as REF-fixed, so
  # take the first two fields and let a missing second field become NA.
  ref <- suppressWarnings(as.numeric(sub(",.*$", "", v)))
  alt <- suppressWarnings(as.numeric(sub("^[^,]*,([^,]*).*$", "\\1", v)))
  alt[!grepl(",", v, fixed = TRUE)] <- NA_real_
  dp  <- ref + alt
  set(snp, j = founders[j], value = fifelse(is.na(dp) | dp == 0, NA_real_, alt / dp))
}
snp[, (ad_cols) := NULL]

n_bad <- sum(!complete.cases(snp[, ..founders]))
if (n_bad > 0) {
  cat(sprintf("  dropping %d site(s) with a founder lacking reads\n", n_bad))
  snp <- snp[complete.cases(snp[, ..founders])]
}

# ── Genetic map position, same interpolation the scans use ───────────────────
cat("Adding cM from flymap.r6.txt ...\n")
fm <- fread(file.path(script_dir, "..", "helpfiles", "flymap.r6.txt"),
            header = FALSE, col.names = c("chr", "pos", "cM"))
snp[, cM := NA_real_]
for (ch in unique(snp$CHROM)) {
  fmX <- fm[chr == ch]
  if (nrow(fmX) == 0L) { cat("  no map for", ch, "- cM left NA\n"); next }
  sm  <- ksmooth(fmX$pos, fmX$cM, kernel = "normal", bandwidth = 3e6)
  f_x <- splinefun(sm$x, sm$y)
  snp[CHROM == ch, cM := f_x(POS)]
}

setcolorder(snp, c("CHROM", "POS", founders, "cM"))
setorder(snp, CHROM, POS)

cat("Writing", parsed$out, "...\n")
fwrite(snp, parsed$out, sep = "\t", compress = "gzip")

cat(sprintf("\nDone. %d SNPs x %d founders.\n", nrow(snp), length(founders)))
cat(sprintf("Pass to snp_scan with:\n  --snp-table %s \\\n  --founders  %s\n",
            parsed$out, paste(founders, collapse = ",")))
