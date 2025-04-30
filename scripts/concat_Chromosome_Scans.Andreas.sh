#!/bin/bash
module load R/4.2.2
mydir=$1
name=$(basename $mydir)
Rscript scripts/concat_Chromosome_Scans.Andreas.R $mydir
tar -czvf $mydir/${name}.tar.gz -C $mydir ${name}.pseudoscan.txt ${name}.meansBySample.txt ${name}.5panel.cM.png ${name}.5panel.Mb.png ${name}.Manhattan.png

