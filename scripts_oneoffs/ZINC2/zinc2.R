SCAN_FILES <- c(
  "process/XQTL_talk/ZINC2/ZINC2_M_freqs250/ZINC2_M_freqs250.pseudoscan.txt",
  "process/XQTL_talk/ZINC2/ZINC2_F_freqs250/ZINC2_F_freqs250.pseudoscan.txt"
)
SCAN_LABELS  <- c("Males", "Females")
SCAN_COLOURS <- c("#1F78B4", "#E31A1C")
THRESHOLD    <- 10
OUT_FILE     <- "figures/slides/zinc2.png"
FORMAT       <- "powerpoint"
OUT_HEIGHT_IN <- 5.5
PEAKS        <- NULL
GENES        <- NULL

source("plot_pseudoscan.R")
