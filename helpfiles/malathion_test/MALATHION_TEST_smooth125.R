SCAN_FILES   <- c("process/malathion_test/MALATHION_TEST_smooth125/MALATHION_TEST_smooth125.scan.txt")
SCAN_COLOURS <- c("#1B9E77")
SCAN_LABELS  <- NULL
THRESHOLD    <- 20
OUT_FILE     <- "process/malathion_test/MALATHION_TEST_smooth125/MALATHION_TEST_smooth125.png"
FORMAT       <- "powerpoint"
PEAKS        <- NULL
GENES        <- data.frame(
    name   = c("Ace",   "Cyp6g1", "Mdr65"),
    chr    = c("chr3R", "chr2R",  "chr3L"),
    pos_mb = c(9.07,    12.19,    6.24)
)

source("scripts/plot_pseudoscan.R")
