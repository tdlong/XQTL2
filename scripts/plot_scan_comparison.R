#!/usr/bin/env Rscript
# Three-way scan comparison for a single chromosome:
#   points  = unsmoothed (haps2scan legacy)
#   line 1  = smoothed, no correction (back-calculated from v4)
#   line 2  = smoothed + R² correction (v4 directly)
#
# Run from XQTL2-dev root:
#   Rscript pipeline/scripts/temp/plot_scan_comparison.R \
#       --unsmoothed process/ZINC2/ZINC2_F_unsmoothed/ZINC2_F_unsmoothed.scan.txt \
#       --smoothed   process/ZINC2/ZINC2_F_v4/ZINC2_F_v4.scan.txt \
#       --chr        chr2L \
#       --out        process/ZINC2/scan_comparison_chr2L.png

suppressPackageStartupMessages(library(tidyverse))

args   <- commandArgs(trailingOnly=TRUE)
parsed <- list(chr="chr2L", R2=0.937, df=7)
i <- 1L
while (i <= length(args)) {
  switch(args[i],
    "--unsmoothed" = { parsed$unsmoothed <- args[i+1]; i <- i+2L },
    "--smoothed"   = { parsed$smoothed   <- args[i+1]; i <- i+2L },
    "--chr"        = { parsed$chr <- args[i+1]; i <- i+2L },
    "--out"        = { parsed$out <- args[i+1]; i <- i+2L },
    "--R2"         = { parsed$R2  <- as.numeric(args[i+1]); i <- i+2L },
    stop(paste("Unknown argument:", args[i]))
  )
}

uns <- read.table(parsed$unsmoothed, header=TRUE) %>%
  filter(chr == parsed$chr, !is.na(Wald_log10p))

smo <- read.table(parsed$smoothed, header=TRUE) %>%
  filter(chr == parsed$chr, !is.na(Wald_log10p)) %>%
  mutate(
    tstat_corrected = qchisq(10^(-Wald_log10p), df=parsed$df, lower.tail=FALSE),
    tstat_raw       = tstat_corrected * parsed$R2,
    lp_uncorrected  = -log10(pchisq(tstat_raw, df=parsed$df, lower.tail=FALSE)),
    lp_corrected    = Wald_log10p
  )

xvar <- if ("cM" %in% names(uns)) "cM" else "pos"
xlab <- if (xvar == "cM") "Position (cM)" else "Position (Mb)"
xdiv <- if (xvar == "pos") 1e6 else 1

ylim <- c(0, max(uns$Wald_log10p, smo$lp_corrected, na.rm=TRUE) * 1.05)

png(parsed$out, width=1400, height=500, res=120)
par(mar=c(4,4,3,1))

plot(uns[[xvar]] / xdiv, uns$Wald_log10p,
     pch=16, cex=0.3, col=adjustcolor("steelblue", 0.4),
     ylim=ylim, xlab=xlab, ylab="-log10(p)",
     main=sprintf("%s  ZINC2 Female", parsed$chr))

lines(smo[[xvar]] / xdiv, smo$lp_uncorrected,
      col="darkorange", lwd=1.5, lty=2)

lines(smo[[xvar]] / xdiv, smo$lp_corrected,
      col="darkred", lwd=2)

legend("topright",
       legend=c("unsmoothed",
                "smoothed (no correction)",
                sprintf("smoothed + R²=%.3f correction", parsed$R2)),
       pch=c(16, NA, NA), lty=c(NA, 2, 1), lwd=c(NA, 1.5, 2),
       col=c(adjustcolor("steelblue", 0.6), "darkorange", "darkred"),
       pt.cex=0.8, bty="n", cex=0.9)

dev.off()
cat("Saved", parsed$out, "\n")
