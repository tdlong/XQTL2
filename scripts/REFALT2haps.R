library(tidyverse)
library(limSolve)
library(data.table)
library(dtplyr)
# library(tidyfast)
args = commandArgs(trailingOnly=TRUE)
mychr=as.character(args[1])
myparfile = as.character(args[2])
mydir = as.character(args[3])

source(myparfile)
# Stage layout: ref/alt in Calls/, haplotypes in Haps/ (see reorganize_project.sh).
filein=paste0(mydir,"/Calls/RefAlt.",mychr,".txt")
rdsfile=paste0(mydir,"/Haps/R.haps.",mychr,".rds")
fileout=paste0(mydir,"/Haps/R.haps.",mychr,".out.rds")
# Guide a project that predates the Calls/Haps/Scans layout to the migration script.
if (!file.exists(filein) && file.exists(paste0(mydir,"/RefAlt.",mychr,".txt")))
  stop("This project uses the old flat layout (RefAlt in ", mydir, "/). Run:\n",
       "  bash pipeline/scripts/reorganize_project.sh ", mydir, "\nthen rerun.")
dir.create(paste0(mydir,"/Haps"), showWarnings=FALSE)
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE))))
source(file.path(script_dir, "REFALT2haps.code.R"))

