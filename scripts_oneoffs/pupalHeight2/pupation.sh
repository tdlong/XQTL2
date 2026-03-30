#!/bin/bash
# Run from repo root: bash scripts_oneoffs/pupalHeight2/pupation.sh
Rscript scripts/plot_pseudoscan.R \
  --scan    process/XQTL_talk/PUPATION/PUPATION_TB_freqs250/PUPATION_TB_freqs250.pseudoscan.txt \
  --label   "Top vs Bottom" \
  --colour  black \
  --threshold 10 \
  --out     figures/slides/pupation.png \
  --format  powerpoint \
  --height  5.8
