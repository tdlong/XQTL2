#!/usr/bin/env Rscript
###############################################################################
# chrX_snp_manhattan.R  —  chrX-only SNP Manhattan: haplotype scan lines
#                           with S/S SNPs overlaid (pA grey, pB-F red)
#
# Reads per-chr scan files directly (no concat needed).
# Thresholds: pA >= 7, pB >= 10 (genome-wide, no chrX adjustment).
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

PA_SCAN  <- "process/ZINC_Hanson/ZINC_Hanson_v3/ZINC_Hanson_v3.scan.chrX.txt"
PB_SCAN  <- "process/ZINC2/ZINC2_F_v3/ZINC2_F_v3.scan.chrX.txt"
SNP_TAB  <- "output/snp_shared_poly.tsv"
OUT      <- "output/chrX_snp_manhattan.png"

THRESH_PA <- 7
THRESH_PB <- 10

read_wt <- function(path) {
  hdr <- scan(path, what="", nlines=1, quiet=TRUE)
  df  <- read.table(path, header=FALSE, skip=1, col.names=c("rn", hdr))
  as_tibble(df[, -1])
}

cat("Reading haplotype scans...\n")
pA <- read_wt(PA_SCAN) %>% mutate(panel="pA (Hanson)", pos_mb=pos/1e6)
pB <- read_wt(PB_SCAN) %>% mutate(panel="pB females", pos_mb=pos/1e6)
cat(sprintf("  pA chrX: %d windows, max Wald=%.1f\n", nrow(pA), max(pA$Wald_log10p)))
cat(sprintf("  pB chrX: %d windows, max Wald=%.1f\n", nrow(pB), max(pB$Wald_log10p)))

scans <- bind_rows(pA, pB)

cat("Reading S/S SNPs...\n")
ss <- read_tsv(SNP_TAB, show_col_types=FALSE) %>%
  filter(quadrant == "S/S", chr == "chrX") %>%
  mutate(pos_mb = pos / 1e6)
cat(sprintf("  S/S SNPs on chrX: %d\n", nrow(ss)))

ss_long <- ss %>%
  select(pos_mb, Wald_log10p_pA, Wald_log10p_pB) %>%
  pivot_longer(c(Wald_log10p_pA, Wald_log10p_pB),
               names_to="panel", values_to="Wald_log10p") %>%
  mutate(panel = recode(panel,
                        Wald_log10p_pA="pA (Hanson)",
                        Wald_log10p_pB="pB females"))

panel_cols <- c("pA (Hanson)"="grey50", "pB females"="#E31A1C")

p <- ggplot() +
  geom_line(data=scans,
            aes(x=pos_mb, y=Wald_log10p, colour=panel),
            linewidth=0.5, alpha=0.85) +
  geom_hline(yintercept=THRESH_PA, colour="grey50",  linetype="dotted", linewidth=0.4) +
  geom_hline(yintercept=THRESH_PB, colour="#E31A1C", linetype="dotted", linewidth=0.4) +
  geom_point(data=ss_long,
             aes(x=pos_mb, y=Wald_log10p, colour=panel),
             size=0.5, alpha=0.4, show.legend=FALSE) +
  scale_colour_manual(values=panel_cols) +
  scale_x_continuous(breaks=seq(0, 24, 2)) +
  labs(x="Position (Mb)",
       y=expression(-log[10](p)~"(haplotype Wald)"),
       colour=NULL,
       title="chrX — haplotype scans with S/S SNPs overlaid",
       subtitle=sprintf("pA threshold=%.0f  pB threshold=%.0f  |  S/S SNPs: %d",
                        THRESH_PA, THRESH_PB, nrow(ss))) +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.position="top")

ggsave(OUT, p, width=10, height=4, dpi=150)
cat("Written:", OUT, "\n")
