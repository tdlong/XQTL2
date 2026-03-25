SCAN_FILES <- c(
  "process/XQTL_talk/MALATHION/MALATHION_freqs250/MALATHION_freqs250.pseudoscan.txt"
)
SCAN_LABELS  <- c("Malathion selected vs control")
SCAN_COLOURS <- c("black")
THRESHOLD    <- 10
OUT_FILE     <- "figures/slides/malathion.png"
FORMAT       <- "powerpoint"
OUT_HEIGHT_IN <- 6.0
PEAKS        <- NULL
GENES <- data.frame(
  name   = c("Ace",    "Cyp6g1", "Mdr65"),
  chr    = c("chr3R",  "chr2R",  "chr3L"),
  pos_bp = c(9069721,  12185667, 6240500)
)

source("plot_pseudoscan.R")
