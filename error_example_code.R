library(tidyverse)
library(limSolve)
library(abind)

source("scripts/scan_functions.R")

names_in_bam=c("Con", "FHV")
samples=c("Con_1_1", "FHV_1_1")
Numflies = data.frame(pool=samples,Num=c(100,100))
ProportionSelect = data.frame(REP=c(1),Proportion=c(0.10))
TreatmentMapping = data.frame(longTRT=c("Con","FHV"),TRT=c("C","Z"))

recodeTable = as_tibble(data.frame(sample=names_in_bam,pool=samples))

xx1 = readRDS("R.haps.chr2L.out.rds")
Nfounders=length(xx1$Groups[[1]][[1]])

bb1 = xx1 %>%
        group_by(CHROM,pos) %>%
        nest() %>%
        mutate(out = map2(data, CHROM, doscan, Nfounders=Nfounders))

