#!/usr/bin/env Rscript
###############################################################################
# chrX_freq_check.R  —  Plot smoothed founder frequencies across chrX
#                        for ZINC2_F_v3, all 8 founders in facets.
#                        C=grey, Z=red; all reps overlaid as thin lines.
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

MEANS <- "process/ZINC2/ZINC2_F_v3/ZINC2_F_v3.meansBySample.chrX.txt"
OUT   <- "output/chrX_freq_check.png"

cat("Reading", MEANS, "\n")
df <- read.table(MEANS, header=TRUE) %>%
  as_tibble() %>%
  mutate(pos_mb = pos / 1e6,
         TRT    = factor(TRT, levels=c("C","Z")))

cat(sprintf("  %d rows | founders: %s\n", nrow(df),
            paste(sort(unique(df$founder)), collapse=", ")))

p <- ggplot(df, aes(x=pos_mb, y=freq, group=interaction(TRT, REP),
                    colour=TRT, alpha=TRT)) +
  geom_line(linewidth=0.3) +
  scale_colour_manual(values=c(C="grey60", Z="#CC3333"), name="Treatment") +
  scale_alpha_manual(values=c(C=0.5, Z=0.7), guide="none") +
  facet_wrap(~founder, ncol=2, scales="free_y") +
  labs(x="Position (Mb)", y="Smoothed frequency",
       title="ZINC2_F_v3: smoothed founder frequencies on chrX",
       subtitle="All replicates; C=grey, Z=red") +
  theme_bw(base_size=10) +
  theme(panel.grid.minor=element_blank(),
        strip.text=element_text(face="bold"))

ggsave(OUT, p, width=10, height=10, dpi=150)
cat("Written:", OUT, "\n")
