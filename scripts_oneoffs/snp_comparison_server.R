library(tidyverse)

# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
# Submitted by rerun_v4.sh after ZINC_Hanson_v4 and ZINC2_F_v4 figs complete.

PA_PREFIX <- "ZINC_Hanson_v3"
PA_DIR    <- "process/ZINC_Hanson/ZINC_Hanson_v3"
PB_PREFIX <- "ZINC2_F_v3"
PB_DIR    <- "process/ZINC2/ZINC2_F_v3"
PA_SCAN   <- file.path(PA_DIR, paste0(PA_PREFIX, ".snp_scan.txt"))
PB_SCAN   <- file.path(PB_DIR, paste0(PB_PREFIX, ".snp_scan.txt"))
OUT_TAB   <- "output/snp_shared_poly.tsv"

POLY_LO   <- 0.05
POLY_HI   <- 0.95
THRESH_PA <- 7
THRESH_PB <- 10
CHRS      <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

# Pipeline writes with write.table(row.names=TRUE)
read_wt <- function(path) {
    hdr <- scan(path, what="", nlines=1, quiet=TRUE)
    df  <- read.table(path, header=FALSE, skip=1, col.names=c("rn", hdr))
    as_tibble(df[, -1])
}

# ── 1. Polymorphism filter — per-chr to avoid loading all at once ────────────
poly_from_scan <- function(prefix, dir, label) {
    cat(sprintf("  %s:\n", label))
    map_dfr(CHRS, function(chr) {
        f <- file.path(dir, sprintf("%s.snp_meansBySample.%s.txt", prefix, chr))
        cat(sprintf("    %s\n", chr))
        read_wt(f) %>%
            group_by(chr, pos, TRT) %>%
            summarize(mean_freq = mean(F_alt, na.rm=TRUE), .groups="drop") %>%
            pivot_wider(names_from=TRT, values_from=mean_freq,
                        names_prefix="mean_") %>%
            filter(mean_C >= POLY_LO, mean_C <= POLY_HI,
                   mean_Z >= POLY_LO, mean_Z <= POLY_HI) %>%
            select(chr, pos)
    })
}

cat("Computing polymorphic SNP sets...\n")
pA_poly <- poly_from_scan(PA_PREFIX, PA_DIR, "pA")
pB_poly <- poly_from_scan(PB_PREFIX, PB_DIR, "pB")

cat(sprintf("pA polymorphic: %d  pB polymorphic: %d\n",
            nrow(pA_poly), nrow(pB_poly)))

shared_poly <- inner_join(pA_poly, pB_poly, by=c("chr","pos"))
cat(sprintf("Shared polymorphic: %d\n", nrow(shared_poly)))

# ── 2. Join with Wald scores and classify quadrants ─────────────────────────
cat("Reading SNP scans...\n")
pA_scan <- read_wt(PA_SCAN) %>% select(chr, pos, Wald_log10p)
pB_scan <- read_wt(PB_SCAN) %>% select(chr, pos, Wald_log10p)

shared <- inner_join(pA_scan, shared_poly, by=c("chr","pos")) %>%
    inner_join(pB_scan, by=c("chr","pos"), suffix=c("_pA","_pB")) %>%
    mutate(
        sig_pA   = Wald_log10p_pA >= THRESH_PA,
        sig_pB   = Wald_log10p_pB >= THRESH_PB,
        quadrant = case_when(
             sig_pA &  sig_pB ~ "S/S",
             sig_pA & !sig_pB ~ "S/NS",
            !sig_pA &  sig_pB ~ "NS/S",
            TRUE               ~ "NS/NS"
        )
    )
cat(sprintf("SNPs in final table: %d\n", nrow(shared)))

# ── 3. Save ──────────────────────────────────────────────────────────────────
dir.create("output", showWarnings=FALSE)
write_tsv(shared, OUT_TAB)
cat("Written:", OUT_TAB, "\n")

shared %>% count(quadrant) %>% arrange(quadrant) %>% print()
