library(tidyverse)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(slider)

source("scripts/XQTL_plotting_functions.R")

args = commandArgs(trailingOnly=TRUE)
data_path = as.character(args[1])
name = basename(data_path)

# merge scan
files = dir(data_path, pattern = paste0(name,".pseudoscan.*.txt")) # get file names
files = grep("chr",files,value=TRUE)
df1 = files %>%
  # read in all the files, appending the path before the filename
  map(~ as_tibble(read.table(file.path(data_path, .)))) %>% 
  reduce(rbind) %>%
  # Apply sliding window smoothing to summary statistics (5 before, current, 5 after)
  group_by(chr) %>%
  arrange(chr, pos) %>%
  mutate(
    Wald_log10p = slide_dbl(Wald_log10p, mean, .before = 5, .after = 5),
    Pseu_log10p = slide_dbl(Pseu_log10p, mean, .before = 5, .after = 5),
    Falc_H2 = slide_dbl(Falc_H2, mean, .before = 5, .after = 5),
    Cutl_H2 = slide_dbl(Cutl_H2, mean, .before = 5, .after = 5),
    avg.var = slide_dbl(avg.var, mean, .before = 5, .after = 5)
  ) %>%
  ungroup()

write.table(df1,paste0(data_path, "/",name,".pseudoscan.txt"))

# merge means
files = dir(data_path, pattern = paste0(name,".meansBySample.*.txt")) # get file names
files = grep("chr",files,value=TRUE)

df2 = files %>%
  # read in all the files, appending the path before the filename
  map(~ as_tibble(read.table(file.path(data_path, .)))) %>% 
  reduce(rbind)

write.table(df2,paste0(data_path,"/",name,".meansBySample.txt"))

p1 = XQTL_Manhattan_5panel(df1, cM = FALSE)
p2 = XQTL_Manhattan_5panel(df1, cM = TRUE)
p3 = XQTL_Manhattan(df1, cM = FALSE)
p4 = XQTL_Manhattan(df1, cM = TRUE)

png(paste0(data_path,"/",name,".5panel.cM.png"), width=8, height=8, units="in", res=600)
p2
dev.off()

png(paste0(data_path,"/",name,".5panel.Mb.png"), width=8, height=8, units="in", res=600)
p1
dev.off()

png(paste0(data_path,"/",name,".Manhattan.png"), width=8, height=8, units="in", res=600)
p3 / p4
dev.off()

