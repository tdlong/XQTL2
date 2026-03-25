#!/usr/bin/env Rscript
###############################################################################
# diag_freqsmooth.R  (server-side, self-contained)
#
# Generates a 3-panel chr2R diagnostic figure comparing frequency-smoothing
# bandwidths.  Run after the 4 freqsmooth scan jobs complete.
#
# Run from /dfs7/adl/tdlong/fly_pool/XQTL2:
#   Rscript scripts/diag_freqsmooth.R
#
# Output: process/ZINC2_v2/freqsmooth_diag.png
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

BASE <- "process/ZINC2_v2"
THRESHOLD <- 20

# euchromatin boundaries for chr2R
EU_START <- 5.398
EU_END   <- 24.685

# colours
COL_F <- "#CC3333"
COL_M <- "#3366CC"

read_scan <- function(path) {
  read.table(path, header = TRUE) %>%
    as_tibble() %>%
    mutate(pos_mb = pos / 1e6)
}

panels <- list(
  "No freq smooth\n(cov smooth only)" = list(
    F = file.path(BASE, "ZINC2_F_stable",  "ZINC2_F_stable.pseudoscan.chr2R.txt"),
    M = file.path(BASE, "ZINC2_M_stable",  "ZINC2_M_stable.pseudoscan.chr2R.txt")
  ),
  "+ freq smooth 100 kb" = list(
    F = file.path(BASE, "ZINC2_F_freqs100", "ZINC2_F_freqs100.pseudoscan.chr2R.txt"),
    M = file.path(BASE, "ZINC2_M_freqs100", "ZINC2_M_freqs100.pseudoscan.chr2R.txt")
  ),
  "+ freq smooth 250 kb" = list(
    F = file.path(BASE, "ZINC2_F_freqs250", "ZINC2_F_freqs250.pseudoscan.chr2R.txt"),
    M = file.path(BASE, "ZINC2_M_freqs250", "ZINC2_M_freqs250.pseudoscan.chr2R.txt")
  )
)

dat <- lapply(names(panels), function(label) {
  bind_rows(
    read_scan(panels[[label]]$F) %>% mutate(sex = "Female"),
    read_scan(panels[[label]]$M) %>% mutate(sex = "Male")
  ) %>% mutate(panel = label)
}) %>%
  bind_rows() %>%
  mutate(
    panel = factor(panel, levels = names(panels)),
    sex   = factor(sex,   levels = c("Female", "Male"))
  )

het_rects <- tibble(
  xmin = c(-Inf, EU_END), xmax = c(EU_START, Inf),
  ymin = -Inf,             ymax = Inf
)

p <- ggplot(dat, aes(x = pos_mb, y = Wald_log10p, colour = sex)) +
  geom_rect(data = het_rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey90", alpha = 0.6) +
  geom_point(size = 0.3, alpha = 0.6, shape = 16) +
  geom_hline(yintercept = THRESHOLD, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  scale_colour_manual(values = c(Female = COL_F, Male = COL_M), name = NULL) +
  scale_x_continuous(breaks = seq(5, 25, 5),
                     limits = c(4, 26),
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  labs(
    title = "chr2R Wald profiles — frequency smoothing diagnostic",
    x     = "Position (Mb)",
    y     = "-log10(p) Wald"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", hjust = 0),
    legend.position  = "top",
    panel.spacing    = unit(0.5, "lines")
  )

outfile <- file.path(BASE, "freqsmooth_diag.png")
ggsave(outfile, p, width = 7, height = 8, dpi = 150)
cat("Saved:", outfile, "\n")
