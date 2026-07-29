library(tidyverse)
library(limSolve)
library(abind)
args = commandArgs(trailingOnly=TRUE)
mychr=as.character(args[1])
myRfile = as.character(args[2])
mydir = as.character(args[3])
myoutdir = as.character(args[4])
design.df = read.table(myRfile)
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE))))
# This script lives in scripts/legacy/; the shared library stays in scripts/.
source(file.path(script_dir, "..", "scan_functions.R"))
# Stage layout: haplotypes in Haps/, scans under Scans/ (see reorganize_project.sh).
filein=paste0(mydir,"/Haps/R.haps.",mychr,".out.rds")
# Guide a project that predates the Calls/Haps/Scans layout to the migration script.
if (!file.exists(filein) && file.exists(paste0(mydir,"/R.haps.",mychr,".out.rds")))
  stop("This project uses the old flat layout (R.haps in ", mydir, "/). Run:\n",
       "  bash pipeline/scripts/reorganize_project.sh ", mydir, "\nthen rerun.")
mydirout=paste0(mydir,"/Scans/",myoutdir)
dir.create(mydirout, showWarnings=FALSE, recursive=TRUE)
fileout=paste0(mydirout,"/",myoutdir,".scan.",mychr,".txt")
fileout_meansBySample=paste0(mydirout,"/",myoutdir,".meansBySample.",mychr,".txt")
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
source(file.path(script_dir, "haps2scan.Apr2025.code.R"))

