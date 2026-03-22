#!/usr/bin/env Rscript
# plot_meansBySample_QC.R
#
# Quick QC plot: one founder's smoothed frequency vs position,
# control vs selected, for a single chromosome.
#
# Edit the four parameters below, then:
#   Rscript helpfiles/malathion_test/plot_meansBySample_QC.R

library(tidyverse)

MEANS_FILE <- "MALATHION_TEST_v2_smooth125.hap/MALATHION_TEST_v2_smooth125.meansBySample.txt"
FOUNDER    <- "A3"
CHR        <- "chr3L"
OUT_FILE   <- file.path(dirname(MEANS_FILE), paste0("qc_means_", FOUNDER, "_", CHR, ".png"))

df <- read.table(MEANS_FILE, header = TRUE) %>%
  filter(founder == FOUNDER, chr == CHR) %>%
  mutate(
    pos_mb  = pos / 1e6,
    Pool    = ifelse(TRT == "C", "Control", "Selected"),
    Rep     = factor(REP)
  )

if (nrow(df) == 0) stop("No data found — check FOUNDER and CHR values.")

p <- ggplot(df, aes(x = pos_mb, y = freq, colour = Rep, linetype = Pool)) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  scale_linetype_manual(values = c(Control = "dashed", Selected = "solid")) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    title    = paste(FOUNDER, "smoothed frequency —", CHR),
    x        = "Position (Mb)",
    y        = "Founder frequency",
    colour   = "Replicate",
    linetype = "Treatment"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

ggsave(OUT_FILE, p, width = 8, height = 3.5, dpi = 150)
cat("Saved:", OUT_FILE, "\n")
