library(tidyverse)
data_path <- "process/ZINC2/ZINC2_F_v3"
name      <- "ZINC2_F_v3"
mfiles <- dir(data_path, pattern=paste0(name, "\.snp_meansBySample\..*\.txt"), full.names=TRUE)
mfiles <- grep("chr", mfiles, value=TRUE)
cat("Merging", length(mfiles), "snp_meansBySample files\n")
mdf <- mfiles %>% map(~ as_tibble(read.table(.x, header=TRUE))) %>% bind_rows()
moutfile <- file.path(data_path, paste0(name, ".snp_meansBySample.txt"))
write.table(mdf, moutfile, quote=FALSE)
cat("Written:", moutfile, "\n")
cat(sprintf("Total rows: %d\n", nrow(mdf)))
