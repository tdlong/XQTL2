library(tidyverse)
library(patchwork)
library(ggplot2)
library(RColorBrewer)

source("scripts/XQTL_plotting_functions.R")

scan_dir <- "process/XQTL_talk/AGING/AGE_SY20_M_freqs125"
name     <- "AGE_SY20_M_freqs125"

df1 <- as_tibble(read.table(file.path(scan_dir, paste0(name, ".pseudoscan.txt"))))
df2 <- as_tibble(read.table(file.path(scan_dir, paste0(name, ".meansBySample.txt"))))

# ── Zoom to the chr3L peak ────────────────────────────────────────────────────
# Broad window first; drop = 30 log10p units on each side to frame the QTL
out <- XQTL_zoom(df1, "chr3L", 8e6, 11e6, left_drop = 30, right_drop = 30)
cat(sprintf("Peak window: %s  %d – %d bp\n", out$chr, out$start, out$stop))

# ── Plot ──────────────────────────────────────────────────────────────────────
fig <- XQTL_change_average(df2, out$chr, out$start, out$stop)

# ── Save ──────────────────────────────────────────────────────────────────────
outfile <- file.path(scan_dir, paste0(name, ".chr3L_QTL.png"))
png(outfile, width = 7, height = 4, units = "in", res = 300)
print(fig)
dev.off()
cat("Saved:", outfile, "\n")
